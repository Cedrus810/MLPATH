"""One controlled experiment: probe, free NVE response, unbiased quench."""
from dataclasses import dataclass
from pathlib import Path
import traceback
import numpy as np
from ase import units
from ase.io import write
from ase.io.trajectory import Trajectory
from ase.md.verlet import VelocityVerlet
from ase.optimize import FIRE
from .internal import named_values
from .io import append_json, snapshot
from .perturbations import (apply as apply_probe, pair_axis, rotatable_torsions,
                            step_torsion)
from .reliability import GateRejected, GuardedCalculator, PairPulse
from .state import classify, encode, same_basin, validate_atoms


@dataclass
class Outcome:
    record: dict
    endpoint: object = None


class PathRecorder:
    """Per-frame path data for the stages that are not molecular dynamics.

    The trial already writes a trajectory and an observation row per sample, so the probe
    and the free response are on the record. The climb to a saddle, the descents off it and
    the boundary continuation wrote nothing at all -- a walk that ended in the wrong basin
    left behind its verdict and no way to see why. Every one of those stages now records the
    same fields on the same named coordinates, one file per stage, so "where did it go and
    what was the reaction coordinate doing" is answerable from the run directory.

    Nothing here triggers a calculation. Energy and force are recorded when the caller
    already has them or when the attached calculator has them cached; a recorder that
    evaluated the potential itself would change the cost of the walk it is observing.
    """

    def __init__(self, directory, config, stage):
        self.directory = Path(directory)
        self.config = config
        self.stage = stage
        self.frames = 0
        # Fields the caller wants on every row from here on -- which side of a saddle, which
        # relay hop. Set by the caller because the walker does not know it is one of several.
        self.context = {}
        self.directory.mkdir(parents=True, exist_ok=True)

    @staticmethod
    def _cached(atoms, name):
        results = getattr(atoms.calc, "results", None) or {}
        return results.get(name)

    def observe(self, phase, step, atoms, energy=None, forces=None, extra=None):
        energy = energy if energy is not None else self._cached(atoms, "energy")
        if forces is None:
            forces = self._cached(atoms, "forces")
        fmax = (float(np.linalg.norm(np.asarray(forces), axis=1).max())
                if forces is not None else None)
        row = {"frame": self.frames, "stage": self.stage, "phase": phase, "step": step,
               "potential_eV": None if energy is None else float(energy),
               "max_force_eV_A": fmax,
               "tracked": named_values(atoms.positions, self.config.tracked_coordinates),
               "edges": [list(e) for e in sorted(
                   encode(atoms, self.config.bond_scale,
                          active=self.config.active_atoms).edges)]}
        row.update(self.context)
        if extra:
            row.update(extra)
        append_json(self.directory / f"{self.stage}.jsonl", row)
        frame = atoms.copy()
        frame.calc = None
        frame.info = {"stage": self.stage, "phase": phase, "step": step}
        write(str(self.directory / f"{self.stage}.extxyz"), frame, append=self.frames > 0)
        self.frames += 1
        return row


def thermal_momenta(atoms, temperature, seed):
    rng = np.random.default_rng(seed)
    masses = atoms.get_masses()
    p = rng.normal(size=(len(atoms), 3)) * np.sqrt(masses[:, None] * units.kB * temperature)
    p -= masses[:, None] * p.sum(axis=0) / masses.sum()
    positions = atoms.positions - atoms.get_center_of_mass()
    inertia = sum(m * (np.dot(r, r) * np.eye(3) - np.outer(r, r))
                  for m, r in zip(masses, positions))
    angular = np.cross(positions, p).sum(axis=0)
    omega = np.linalg.pinv(inertia) @ angular
    p -= masses[:, None] * np.cross(omega, positions)
    atoms.set_momenta(p)


def _trivial_modes(positions, masses):
    """Mass-weighted translations and rotations, orthonormalized.

    These span the directions along which a converged structure is free by symmetry;
    projecting them out is what lets the lowest remaining eigenvalue be compared
    against zero without a threshold that has to absorb rigid-body noise.
    """
    root = np.sqrt(np.repeat(masses, 3))
    centered = positions - np.average(positions, axis=0, weights=masses)
    modes = []
    for axis in range(3):
        translation = np.zeros_like(positions)
        translation[:, axis] = 1.0
        modes.append(translation.ravel() * root)
        rotation = np.cross(np.eye(3)[axis], centered)
        modes.append(rotation.ravel() * root)
    basis = []
    for mode in modes:
        for kept in basis:
            mode = mode - np.dot(mode, kept) * kept
        norm = np.linalg.norm(mode)
        if norm > 1e-8:
            basis.append(mode / norm)
    return np.array(basis)


def hessian_vector_product(atoms, guard, direction, step):
    """-dF along a direction, which is Hv: the curvature the structure shows when pushed.

    Expanding about a reference geometry, F(R + dR) = F(R) - H dR + O(dR^2), so a small
    displacement returns one Hessian-vector product for two force evaluations. The full
    Hessian is 3N of these; an iterative lowest-eigenvalue solver needs a few dozen
    regardless of size; and a single one gives the stiffness of one direction. It is the
    same measurement the search performs at finite amplitude -- the Hessian is the
    infinitesimal limit of perturbation response, not a separate object.
    """
    direction = np.asarray(direction, dtype=float).reshape(-1, 3)
    reference = atoms.positions.copy()
    gradients = []
    for sign in (1, -1):
        atoms.set_positions(reference + sign * step * direction)
        gradients.append(-np.asarray(guard.get_forces(atoms)))
    atoms.set_positions(reference)
    return (gradients[0] - gradients[1]) / (2 * step)


def _projected_direction(positions, masses, direction):
    """Mass-weighted unit direction with the rigid-body components removed.

    Curvature along a translation is exactly zero and along a rotation is proportional to
    the residual gradient, so leaving either in a probe direction only dilutes whatever
    signal the probe was meant to carry.
    """
    weighted = (np.sqrt(masses[:, None])
                * np.asarray(direction, dtype=float).reshape(-1, 3)).ravel()
    basis = _trivial_modes(positions, masses)
    weighted = weighted - basis.T @ (basis @ weighted)
    norm = np.linalg.norm(weighted)
    if norm < 1e-10:
        return None
    return (weighted / norm).reshape(-1, 3)


def _curvature_along(atoms, guard, direction, step):
    """v^T H v in mass-weighted coordinates for one already-mass-weighted unit direction.

    Two force evaluations, independent of system size. Sign is the whole content: negative
    means the reference is not a minimum along this direction.
    """
    masses = atoms.get_masses()
    product = hessian_vector_product(atoms, guard, direction / np.sqrt(masses[:, None]), step)
    return float(np.sum(direction * product / np.sqrt(masses[:, None])))


def _mass_weighted_direction(masses, direction):
    """Unit direction in mass-weighted coordinates, from a Cartesian displacement.

    Mass-weighted coordinates are q = sqrt(m) x, so a Cartesian displacement maps to
    sqrt(m) times itself, not to itself over sqrt(m). Getting this backwards is invisible
    on a homonuclear pair, where it only rescales, and picks the wrong direction entirely
    as soon as the masses differ.
    """
    weighted = np.sqrt(masses[:, None]) * np.asarray(direction, dtype=float).reshape(-1, 3)
    norm = np.linalg.norm(weighted)
    if norm < 1e-12:
        raise ValueError("Direction has no mass-weighted component")
    return weighted / norm


def direction_stiffness(atoms, factory, config, direction):
    """v^T H v in mass-weighted coordinates, for two force evaluations.

    Small and positive means a soft mode, where another basin may be close; large means
    a stiff one that will need a large amplitude; negative means the reference is not a
    minimum along it. The direction is a Cartesian displacement; the value is in
    eV/(Angstrom^2 amu), so it is comparable with a hessian_spectrum eigenvalue.
    """
    masses = atoms.get_masses()
    weighted = _mass_weighted_direction(masses, direction)
    probe = atoms.copy()
    guard = GuardedCalculator(factory(), config, None)
    probe.calc = guard
    product = hessian_vector_product(probe, guard, weighted / np.sqrt(masses[:, None]),
                                     config.minimum_check_step_A)
    return float(np.sum(weighted * product / np.sqrt(masses[:, None])))


def _rigid_body_expectation(positions, masses, gradient, mode):
    """How much of a rigid-body Rayleigh quotient the residual gradient already explains.

    E(R(s)x) = E(x) holds for every s, so differentiating twice gives

        xdot^T H xdot = omega^2 (g . x_perp),

    which means a rotational curvature is proportional to the gradient and vanishes only
    at a stationary point. Read raw, it says nothing about the quality of the calculation:
    a structure converged to the search's own force tolerance shows a rotational quotient
    of order |g| / r, which on ethanol is 1e-3 -- the eigenvalue tolerance itself. What is
    left after subtracting this is the part that is genuinely error, and that residual is
    the quantity worth gating on. Translations have no such term at all: the identity is
    exact at every geometry because the energy depends on coordinate differences.

    Returns (expectation, kind). Measured against a central-difference Hessian at
    non-stationary geometries the prediction is good to four digits, the remainder being
    the finite difference's own chord-for-arc bias.
    """
    root = np.sqrt(np.repeat(masses, 3))
    velocity = (mode / root).reshape(-1, 3)
    scale = max(np.linalg.norm(velocity), 1e-30)
    if np.linalg.norm(velocity - velocity.mean(axis=0)) < 1e-8 * scale:
        return 0.0, "translation"
    centered = positions - np.average(positions, axis=0, weights=masses)
    inertia = sum(np.dot(r, r) * np.eye(3) - np.outer(r, r) for r in centered)
    angular = sum(np.cross(r, v) for r, v in zip(centered, velocity))
    omega = np.linalg.pinv(inertia) @ angular
    speed = float(np.linalg.norm(omega))
    if speed < 1e-12:
        return 0.0, "translation"
    axis = omega / speed
    perpendicular = centered - np.outer(centered @ axis, axis)
    expectation = speed ** 2 * float(np.sum(gradient * perpendicular)) / float(np.sum(mode ** 2))
    return expectation, "rotation"


def trivial_mode_floor(mass_weighted, basis, positions, masses, gradient):
    """Per-mode error floor of a curvature calculation, needing no reference data.

    Translations and rotations carry no curvature on an invariant potential, so whatever
    is measured along them beyond what the gradient explains is error: model
    non-invariance plus whatever the differentiation added. This is the only quality
    measure available without an external reference, and it costs no force evaluations
    beyond the gradient, which the caller already has. A residual that is not far below
    the eigenvalue tolerance means every decision at that scale is noise.
    """
    report = []
    for mode in basis:
        expectation, kind = _rigid_body_expectation(positions, masses, gradient, mode)
        measured = float(mode @ mass_weighted @ mode)
        report.append({"kind": kind, "measured": measured,
                       "gradient_expectation": expectation,
                       "residual": measured - expectation})
    return report


def floor_residual(report):
    """The worst unexplained rigid-body curvature: the number the gate compares."""
    return max(abs(entry["residual"]) for entry in report)


def _hessian_providers(physical):
    """The backend's own analytic second derivative, if every member exposes one."""
    members = getattr(physical, "members", None)
    candidates = list(members) if members else [physical]
    if candidates and all(hasattr(member, "get_hessian") for member in candidates):
        return candidates
    return []


def _analytic_hessian(candidate, physical, config):
    """d2E/dx2 straight from the backend, verified against one finite difference.

    A model whose forces are already an autograd gradient carries the second derivative
    for free, so differencing forces to rebuild it is numerical differentiation of an
    analytic derivative: slower and only accurate to O(step^2). Two traps make this a
    checked path rather than a swap. The unit factors a calculator applies to energy and
    forces are not necessarily applied to its Hessian, and a wrong one is silent whenever
    the factors happen to be 1; and nothing guarantees the array layout. So the result is
    proven by execution -- one Hessian-vector product along a fixed direction must agree
    with the finite-difference measurement -- rather than trusted because the method
    exists. Returns (hessian, provenance) or (None, reason).
    """
    providers = _hessian_providers(physical)
    if not providers:
        return None, "backend exposes no analytic Hessian"
    size = 3 * len(candidate)
    matrices = []
    for provider in providers:
        scale = (float(getattr(provider, "energy_units_to_eV", 1.0))
                 / float(getattr(provider, "length_units_to_A", 1.0)) ** 2)
        try:
            raw = np.asarray(provider.get_hessian(candidate), dtype=float)
        except Exception as exc:                        # noqa: BLE001 - fail closed below
            return None, f"analytic Hessian raised {type(exc).__name__}: {exc}"
        if raw.size != size * size:
            return None, f"analytic Hessian has {raw.size} entries, expected {size * size}"
        matrices.append(scale * raw.reshape(size, size))
    hessian = np.mean(matrices, axis=0)
    hessian = 0.5 * (hessian + hessian.T)

    direction = np.zeros(size)
    direction[::2] = 1.0
    direction /= np.linalg.norm(direction)
    guard = GuardedCalculator(physical, config, None)
    probe = candidate.copy()
    probe.calc = guard
    measured = hessian_vector_product(probe, guard, direction.reshape(-1, 3),
                                      config.minimum_check_step_A).ravel()
    predicted = hessian @ direction
    reference = max(np.abs(measured).max(), 1e-12)
    deviation = float(np.abs(predicted - measured).max() / reference)
    if deviation > config.analytic_hessian_check:
        return None, (f"analytic Hessian disagrees with a finite difference by "
                      f"{deviation:.3e} relative, above {config.analytic_hessian_check}")
    return hessian, {"source": "analytic", "members": len(providers),
                     "finite_difference_deviation": deviation,
                     "meaning": "backend second derivative, cross-checked against one HVP"}


def _mass_weighted_hessian(candidate, factory, config, source="auto"):
    """Mass-weighted Hessian plus the rigid-body basis, from whichever path is available.

    Two sources, and which one was used is recorded rather than assumed: "analytic" is
    the backend's own second derivative and "fd" is 6N force evaluations of central
    differences. "auto" prefers analytic and falls back with the reason kept.
    """
    positions = candidate.positions
    masses = candidate.get_masses()
    count = len(candidate)
    physical = factory()
    hessian, provenance = None, None
    if source in ("auto", "analytic"):
        hessian, provenance = _analytic_hessian(candidate, physical, config)
        if hessian is None:
            if source == "analytic":
                raise GateRejected("analytic_hessian_unavailable", str(provenance), candidate)
            fallback_reason = provenance
    if hessian is None:
        guard = GuardedCalculator(physical, config, candidate.get_potential_energy())
        probe = candidate.copy()
        probe.calc = guard
        step = config.minimum_check_step_A
        hessian = np.zeros((3 * count, 3 * count))
        for column in range(3 * count):
            direction = np.zeros(3 * count)
            direction[column] = 1.0
            hessian[column] = hessian_vector_product(probe, guard, direction, step).ravel()
        hessian = 0.5 * (hessian + hessian.T)
        provenance = {"source": "fd", "step_A": step,
                      "force_evaluations": 6 * count,
                      "meaning": "central differences; truncation error is O(step^2)"}
        if source == "auto":
            provenance["analytic_unavailable"] = fallback_reason
    root = np.sqrt(np.repeat(masses, 3))
    mass_weighted = hessian / np.outer(root, root)
    basis = _trivial_modes(positions, masses)
    return mass_weighted, basis, provenance


def hessian_spectrum(candidate, factory, config, return_vectors=False, source="auto"):
    """Mass-weighted Hessian spectrum, rigid-body modes projected out.

    Eigenvalues are in eV/(Angstrom^2 amu); the sign structure is what classifies the
    stationary point, and negative eigenvalues are also reported as imaginary wavenumbers
    because that is how they are read. This is the only regime in which the numbers are
    genuinely eigenvalues: at finite probe amplitude the response is neither a linear nor
    a symmetric map, so there is no spectrum to speak of.
    """
    mass_weighted, basis, _ = _mass_weighted_hessian(candidate, factory, config, source)
    projector = np.eye(mass_weighted.shape[0]) - basis.T @ basis
    projected = projector @ mass_weighted @ projector
    if not return_vectors:
        return np.sort(np.linalg.eigvalsh(projected))
    eigenvalues, vectors = np.linalg.eigh(projected)
    order = np.argsort(eigenvalues)
    return eigenvalues[order], vectors[:, order]


def curvature_spectrum(candidate, factory, config, source="auto"):
    """hessian_spectrum plus its provenance and its reference-free error floor.

    The floor needs the gradient to separate what a non-stationary geometry explains from
    what is error, so this costs one force evaluation more than the spectrum alone -- and
    usually none, because the caller has just quenched and the value is cached.
    """
    mass_weighted, basis, provenance = _mass_weighted_hessian(candidate, factory, config,
                                                              source)
    projector = np.eye(mass_weighted.shape[0]) - basis.T @ basis
    eigenvalues, vectors = np.linalg.eigh(projector @ mass_weighted @ projector)
    order = np.argsort(eigenvalues)
    floor = trivial_mode_floor(mass_weighted, basis, candidate.positions,
                               candidate.get_masses(), -np.asarray(candidate.get_forces()))
    return eigenvalues[order], vectors[:, order], floor, provenance


def _wavenumbers(eigenvalues):
    """Report negative curvature the way it is normally quoted, as icm."""
    scale = units._hbar * 1e10 / np.sqrt(units._e * units._amu)
    return [float(-scale * np.sqrt(-value) / units.invcm)
            for value in eigenvalues if value < 0]


def response_curvature(atoms, guard, direction, config, reference_energy=None):
    """Response curvature along one direction, regressed over a band of amplitudes.

    The symmetric second difference of the energy at amplitude a is not an approximation
    to the curvature that improves as a shrinks; it is the curvature of the potential
    low-passed at wavenumber ~2/a, because in the Fourier variable it multiplies each
    component by sinc^2(ka/2). So a is a resolution, and a single small a is a choice of
    bandwidth made silently. Sampling a band and fitting

        kappa(a) = kappa_0 + c a^2 + O(a^4)

    reports the extrapolated curvature, the leading anharmonicity that the harmonic
    number throws away, and a standard error -- which is the part that matters, because a
    curvature whose magnitude is not several times its own uncertainty cannot classify
    anything. The energy form is used rather than the force form because sinc^2 is
    non-negative: the force form's sinc(ka) inverts the sign of every component with a
    wavelength between a and 2a.

    Amplitudes and the returned curvature are mass-weighted, matching hessian_spectrum:
    A*sqrt(amu) and eV/(Angstrom^2 amu). Costs 2K evaluations for K amplitudes; only the
    energies are used, but each one carries the forces too, because the reliability gate
    needs them to decide whether the probe is physical at all.
    """
    masses = atoms.get_masses()
    weighted = _mass_weighted_direction(masses, direction)
    displacement = weighted / np.sqrt(masses[:, None])
    reference = atoms.positions.copy()
    if reference_energy is None:
        reference_energy = float(guard.get_potential_energy(atoms))
    low, high = config.response_band_A
    amplitudes = np.geomspace(low, high, config.response_samples)
    curvatures = []
    try:
        for amplitude in amplitudes:
            energies = []
            for sign in (1, -1):
                atoms.set_positions(reference + sign * amplitude * displacement)
                energies.append(float(guard.get_potential_energy(atoms)))
            curvatures.append(
                (energies[0] + energies[1] - 2 * reference_energy) / amplitude ** 2)
    finally:
        # A rejected probe must not leave the caller holding a displaced structure.
        atoms.set_positions(reference)
    curvatures = np.array(curvatures)

    design = np.vstack([np.ones_like(amplitudes), amplitudes ** 2]).T
    coefficients, *_ = np.linalg.lstsq(design, curvatures, rcond=None)
    residual = curvatures - design @ coefficients
    degrees = max(1, len(amplitudes) - design.shape[1])
    variance = float((residual ** 2).sum() / degrees)
    covariance = variance * np.linalg.inv(design.T @ design)
    kappa = float(coefficients[0])
    sigma = float(np.sqrt(max(covariance[0, 0], 0.0)))
    signs = np.sign(curvatures)
    return {"kappa": kappa, "sigma": sigma, "anharmonicity": float(coefficients[1]),
            "band_A_sqrt_amu": [float(low), float(high)],
            "samples": [[float(a), float(k)] for a, k in zip(amplitudes, curvatures)],
            "sign_stable": bool(np.all(signs == signs[0])),
            # Deliberately not called significant. sigma is the ordinary-least-squares
            # residual scale of a deterministic fit -- it measures how well kappa_0 + c a^2
            # describes the samples, not a calibrated statistical uncertainty, and it cannot
            # see a systematic error such as probing the wrong direction. Measured on
            # malonaldehyde it read 3.7e-06 while the value was wrong by 0.18.
            "exceeds_residual_scale": bool(abs(kappa)
                                           > config.response_significance * sigma),
            "sigma_meaning": ("OLS residual scale of the a^2 fit, not a calibrated "
                              "uncertainty; blind to a wrong direction"),
            "residual_rms": float(np.sqrt((residual ** 2).mean())),
            "units": "eV/(Angstrom^2 amu); amplitude A*sqrt(amu)",
            "meaning": "band-regressed response curvature; NOT an eigenvalue"}


def response_gradient(atoms, guard, direction, amplitude):
    """The gradient of response_curvature with respect to the direction, for free.

    Differentiating the energy form with respect to v gives exactly twice the odd part of
    the gradient, [g(x+av) - g(x-av)] / a, which the same pair of evaluations already
    produced. So minimizing the response curvature over the unit sphere costs two force
    evaluations per step -- the price of one Hessian-vector product -- with no matrix and
    no spectrum, and reduces to the usual dimer rotation as the amplitude goes to zero.
    Mass-weighted, like response_curvature.
    """
    masses = atoms.get_masses()
    weighted = _mass_weighted_direction(masses, direction)
    cartesian = hessian_vector_product(atoms, guard, weighted / np.sqrt(masses[:, None]),
                                       amplitude)
    return 2.0 * cartesian / np.sqrt(masses[:, None])


def _response_and_gradient(atoms, guard, direction, amplitude, reference_energy):
    """kappa_a^E(v) and its gradient with respect to v, from one pair of evaluations.

    Differentiating the energy form with respect to the direction gives exactly twice the
    odd part of the gradient,

        d/dv [E(x+av) + E(x-av) - 2E(x)] / a^2  =  [g(x+av) - g(x-av)] / a  =  2 D_a(v),

    and the guard returns forces with every energy. So the objective and its gradient come
    out of the same two evaluations, which is what makes rotating a direction downhill in
    response curvature cost what one Hessian-vector product costs -- with no matrix, and
    without needing a spectrum that at finite amplitude does not exist.

    The direction is a mass-weighted unit vector; the gradient is returned in the same
    coordinates. Two force evaluations.
    """
    masses = atoms.get_masses()
    root = np.sqrt(masses[:, None])
    reference = atoms.positions.copy()
    energies, gradients = [], []
    try:
        for sign in (1, -1):
            atoms.set_positions(reference + sign * amplitude * direction / root)
            energies.append(float(guard.get_potential_energy(atoms)))
            gradients.append(-np.asarray(guard.get_forces(atoms)) / root)
    finally:
        atoms.set_positions(reference)
    curvature = (energies[0] + energies[1] - 2 * reference_energy) / amplitude ** 2
    return curvature, (gradients[0] - gradients[1]) / amplitude


def lowest_response_direction(atoms, guard, config, initial, amplitude=None,
                              reference_energy=None):
    """Rotate a direction downhill in response curvature, on the unit sphere.

    This is the dimer idea moved off the Hessian. At finite amplitude there is no operator
    to diagonalize -- the response is neither linear nor symmetric in the direction -- so
    the bottom direction is defined variationally, as a minimizer of kappa_a^E over the
    sphere, and found by descent rather than by an eigensolver. As the amplitude goes to
    zero it reduces to minimizing the Rayleigh quotient, which is what a dimer does.

    Rigid-body components are projected out of the direction and of every step: curvature
    along a translation is zero and along a rotation is proportional to the residual
    gradient, so leaving them in would let the descent drift into directions that carry no
    chemistry. Two force evaluations per step. Convergence is the fraction of the gradient
    still tangential to the sphere, not an absolute norm, because that norm has units and
    no fixed tolerance can serve a soft coordinate and a stiff one at once.

    A converged result is one direction, not an index: it says nothing about how many other
    negative directions exist. Counting them needs the orthogonal complement, and at finite
    amplitude that count depends on the order it is taken in.
    """
    masses = atoms.get_masses()
    amplitude = float(config.response_band_A[0] if amplitude is None else amplitude)
    if reference_energy is None:
        reference_energy = float(guard.get_potential_energy(atoms))
    direction = _projected_direction(atoms.positions, masses, initial)
    if direction is None:
        raise ValueError("Initial direction has no internal component")
    start = direction.copy()

    curvature, gradient = _response_and_gradient(atoms, guard, direction, amplitude,
                                                reference_energy)
    evaluations, step, history = 2, float(config.response_descent_step), []
    converged, stalled = False, False
    for iteration in range(config.response_descent_max_steps):
        tangent = _projected_direction(
            atoms.positions, masses,
            (gradient - np.sum(gradient * direction) * direction) / np.sqrt(masses[:, None]))
        residual = 0.0 if tangent is None else float(
            np.linalg.norm(gradient - np.sum(gradient * direction) * direction))
        # Relative, because the tangent norm carries units: what says the direction is a
        # critical point of the response on the sphere is that hardly any of the gradient
        # is tangential, not that some number in eV/(A^2 amu) is small.
        fraction = residual / max(float(np.linalg.norm(gradient)), 1e-30)
        history.append({"iteration": iteration, "curvature": curvature,
                        "tangent_norm": residual, "tangent_fraction": fraction,
                        "step": step})
        if tangent is None or fraction < config.response_descent_tol:
            converged = True
            break
        accepted = False
        while step > 1e-6:
            trial = _projected_direction(atoms.positions, masses,
                                         (direction - step * tangent)
                                         / np.sqrt(masses[:, None]))
            if trial is None:
                step *= 0.5
                continue
            value, trial_gradient = _response_and_gradient(atoms, guard, trial, amplitude,
                                                           reference_energy)
            evaluations += 2
            if value < curvature:
                direction, curvature, gradient = trial, value, trial_gradient
                step *= 1.5
                accepted = True
                break
            step *= 0.5
        if not accepted:
            # No trial step lowered the objective. That is the line search giving up, not
            # the direction being a critical point, and calling it convergence would
            # report a stall as a result.
            stalled = True
            break
    overlap = float(np.clip(abs(np.sum(direction * start)), 0.0, 1.0))
    return {"direction": direction.tolist(), "curvature": curvature,
            "amplitude_A_sqrt_amu": amplitude, "converged": converged, "stalled": stalled,
            "iterations": len(history), "force_evaluations": evaluations,
            "rotation_from_initial_deg": float(np.degrees(np.arccos(overlap))),
            "final_tangent_norm": history[-1]["tangent_norm"] if history else None,
            "final_tangent_fraction": history[-1]["tangent_fraction"] if history else None,
            "history": history,
            "units": "eV/(Angstrom^2 amu); amplitude A*sqrt(amu)",
            "meaning": ("variational bottom response direction at one amplitude; "
                        "NOT an eigenvector and its value is NOT an eigenvalue")}


def response_persistence(atoms, guard, config, initial, amplitudes=None):
    """Track the bottom response direction across amplitudes.

    A direction that is lowest at one resolution need not be lowest at another, and the
    ethanol syn barrier is a measured case where the sign itself changes with amplitude. So
    "this is the unstable direction" is a statement about a band, not a point, and what
    makes it trustworthy is that the direction moves continuously and the sign does not
    flip inside the band.

    No verdict is returned. The tolerance on how far the direction may rotate cannot be
    set from a molecule whose every torsion is soft, so the numbers are reported and the
    threshold is left to whoever has a benchmark that exercises it.
    """
    if amplitudes is None:
        low, high = config.response_band_A
        amplitudes = np.geomspace(low, high, config.response_samples)
    reference_energy = float(guard.get_potential_energy(atoms))
    results, previous = [], None
    for amplitude in amplitudes:
        found = lowest_response_direction(atoms, guard, config, initial, float(amplitude),
                                          reference_energy)
        direction = np.asarray(found["direction"])
        rotation = None
        if previous is not None:
            overlap = float(np.clip(abs(np.sum(direction * previous)), 0.0, 1.0))
            rotation = float(np.degrees(np.arccos(overlap)))
        results.append({"amplitude_A_sqrt_amu": float(amplitude),
                        "curvature": found["curvature"],
                        "rotation_from_previous_deg": rotation,
                        "force_evaluations": found["force_evaluations"],
                        "converged": found["converged"]})
        previous = direction
    signs = [np.sign(entry["curvature"]) for entry in results]
    return {"band": results,
            "sign_stable": bool(all(s == signs[0] for s in signs)),
            "total_rotation_deg": float(sum(e["rotation_from_previous_deg"] or 0.0
                                            for e in results)),
            "force_evaluations": sum(e["force_evaluations"] for e in results),
            "meaning": "scale-resolved bottom direction; no persistence verdict is implied"}


def follow_min_mode(candidate, factory, config, direction=None, seed=0,
                    ceiling=None, reference=None, recorder=None):
    """Walk to a first-order saddle from a seed near the dividing surface.

    Finding the unstable direction is not finding the saddle. `lowest_response_direction`
    answers where the negative curvature points; getting there still needs the centre to
    move. Min-mode following does that by inverting the force along the lowest mode,

        F_TS = F - 2 (F . v) v,

    so the walk descends in every direction except v and ascends on the ridge. Implemented
    here as the Newton step that expression approximates, taken in the eigenbasis of the
    analytic Hessian: a Newton step converges to the nearest stationary point whatever its
    index, which is the same inversion done exactly and with the curvature as its scale.
    That choice is only available because the analytic Hessian turned out to be clean to
    1e-15 and cheaper than the finite-difference one at these sizes.

    Rigid-body directions are excluded rather than inverted, and the step is capped by a
    trust radius because a Newton step is only as good as the quadratic model.

    Seeding matters more than the optimizer. A quench stalling on a saddle is a
    soft-coordinate coincidence, not an algorithm -- measured on malonaldehyde, sixty
    probes found the product twice and the saddle never, while the saddle sits at
    lambda_1 = -37.3, a ridge 120 times sharper than the ethanol barriers where the
    coincidence did happen. What is PRRS-native is the seed: an amplitude near the dividing
    surface, which the transient crossings now preserve.

    `direction` is the coordinate to follow, in Cartesian components. Supplying it turns
    this from "find the softest saddle" into "find the saddle on this coordinate", which are
    different questions and were being conflated.

    `ceiling` is an energy the walk must stay under. Climbing is unbounded otherwise and
    will happily reach a saddle that has nothing to do with the reaction that prompted it:
    seeded from a malonaldehyde probe frame it reached an O-H homolysis 2.34 eV up, whose
    endpoints are a radical pair and outside the model's domain entirely.

    A step that overshoots the ceiling is backtracked, not fatal. Treating one overshoot as
    the end of the search confuses a step-length problem with a domain problem, and it did:
    the transient seed that had reached the proton-transfer saddle was refused because a
    fixed ascent step carried it over the line on the way. A seed that starts above the
    ceiling is a different verdict again and is skipped so the next seed gets its turn.

    Returns (structure, diagnostics); the structure is None unless the walk reached an
    index-one stationary point.
    """
    atoms = candidate.copy()
    guard = GuardedCalculator(factory(), config, None)
    atoms.calc = guard
    masses = atoms.get_masses()
    root = np.sqrt(masses[:, None])
    tolerance = config.minimum_check_eigenvalue_tol
    history, eigenvalues, vectors = [], None, None
    tracked, selected_order, selected_overlap = None, 0, None
    target, target_overlap = None, None
    if direction is not None:
        # d_q = M^(1/2) d_x normalized, rigid-body part removed. A tangent that survives
        # neither step is not a coordinate to follow, and saying so is the point: silently
        # falling back to the softest mode is how a proton transfer becomes a methyl rotor.
        projected = _projected_direction(atoms.positions, masses, direction)
        if projected is None:
            return None, {"reason": "missing_target_direction", "steps": 0,
                          "meaning": ("the supplied tangent has no internal component, so "
                                      "there is no coordinate to follow; the global lowest "
                                      "mode is a different question and is not substituted")}
        target = projected.ravel()
    # max_rise_eV starts as None, not -inf. It used to start at -inf and get mapped back to
    # None only on the successful return, so all ten early returns carried a -inf straight
    # into the record -- and `atomic_json` writes with allow_nan=False, so the write raised.
    # Worse, the write that raised was publish(), which is also what the failure handler
    # calls, so the run could not even record that it had failed: three seeds of the T=300
    # scan died with "Out of range float values are not JSON compliant: -inf" and left
    # status frozen at "running", which reads exactly like being killed by a signal.
    bounds = {"initial_rise_eV": None, "max_rise_eV": None, "backtracks": 0,
              "ceiling_hit_step": None, "target_supplied": direction is not None,
              "target_overlap": None, "selected_mode_index": None}
    try:
      for step in range(config.min_mode_max_steps):
          if eigenvalues is None or step % config.min_mode_refresh == 0:
              eigenvalues, vectors, floor, provenance = curvature_spectrum(
                  atoms, factory, config, source=config.hessian_source)
              # Which mode to climb. Taking the lowest eigenvalue every time is not
              # min-mode following toward a given reaction -- it is a search for the softest
              # saddle anywhere, and on 3-oxobutanal, whose methyl and hydroxyl rotors are
              # the softest things in the molecule, it returned eight rotor barriers at -65
              # to -80 cm^-1 while the proton transfer sits at -3186. When a target direction
              # is supplied the mode is instead the eigenvector with the largest
              # mass-weighted overlap with the direction being followed, carried forward with
              # a consistent sign, so the walk stays on the coordinate it was asked about.
              # Without a target the question really is "which is the softest saddle", and
              # the lowest mode is the answer; tracking would drift onto whatever mode it
              # started next to. With a target the question is "where is the saddle on this
              # coordinate", and only overlap answers that. The caller picks by supplying
              # `direction` or not, because the two are different questions.
              internal = [i for i, value in enumerate(eigenvalues)
                          if abs(value) > tolerance]
              if not internal:
                  return None, {"reason": "no_internal_mode", "steps": step, **bounds}
              if target is None:
                  # The lowest NON-TRIVIAL mode. Index zero is whichever rigid-body residual
                  # sorts first, and rigid-body directions are skipped below, so selecting it
                  # would leave no mode climbing and turn the walk into a plain minimiser.
                  selected_order, selected_overlap, tracked = internal[0], None, None
              else:
                  # First update: the internal mode most like the tangent, argmax |v_i . d_q|.
                  # Later updates: the one most like the branch already being followed, with the
                  # sign carried so tracking is continuous across an eigenvalue crossing -- the
                  # order index reshuffles as the structure moves and is not a mode's identity.
                  anchor = tracked if tracked is not None else target
                  overlaps = np.abs(np.array([vectors[:, i] @ anchor for i in internal]))
                  order = internal[int(np.argmax(overlaps))]
                  overlap = float(vectors[:, order] @ anchor)
                  tracked = vectors[:, order] * (1.0 if overlap >= 0 else -1.0)
                  selected_order, selected_overlap = order, abs(overlap)
                  if target_overlap is None:
                      target_overlap = selected_overlap
                  if selected_overlap < config.min_mode_overlap:
                      return None, {"reason": "lost_target_mode", "steps": step,
                                    "mode_overlap": selected_overlap,
                                    "selected_mode_index": order,
                                    "target_overlap": target_overlap, **bounds,
                                    "meaning": ("no eigenvector still resembles the "
                                                "coordinate being followed, so continuing "
                                                "would climb a different reaction")}
          # Carried on every return path, not only the successful one: which mode was
          # followed and how well it matched the tangent is exactly what has to be
          # readable when a walk fails or finds the wrong saddle.
          bounds["selected_mode_index"] = selected_order
          bounds["target_overlap"] = target_overlap
          energy = float(guard.get_potential_energy(atoms))
          if reference is not None:
              rise = energy - reference
              bounds["max_rise_eV"] = (rise if bounds["max_rise_eV"] is None
                                       else max(bounds["max_rise_eV"], rise))
              if step == 0:
                  bounds["initial_rise_eV"] = rise
          if ceiling is not None and step == 0 and energy > ceiling:
              # Not a failed climb: a seed that is already outside the region. The next seed
              # deserves its turn rather than inheriting this one's verdict.
              return None, {"reason": "seed_above_ceiling", "steps": 0,
                            "energy_eV": energy, "ceiling_eV": ceiling, **bounds}
          forces = np.asarray(guard.get_forces(atoms))
          fmax = float(np.linalg.norm(forces, axis=1).max())
          history.append({"step": step, "fmax_eV_A": fmax,
                          "lowest_eigenvalue": float(eigenvalues[0]),
                          "followed_eigenvalue": float(eigenvalues[selected_order]),
                          "followed_order": selected_order,
                          "mode_overlap": selected_overlap,
                          "climbing": bool(eigenvalues[selected_order] >= -tolerance)})
          if recorder is not None:
              recorder.observe("climb", step, atoms, energy=energy, forces=forces,
                               extra={"followed_eigenvalue":
                                          float(eigenvalues[selected_order]),
                                      "followed_order": selected_order,
                                      "mode_overlap": selected_overlap,
                                      "climbing": bool(eigenvalues[selected_order]
                                                       >= -tolerance)})
          if fmax < config.min_mode_fmax_eV_A:
              break
          weighted = (forces / root).ravel()
          displacement = np.zeros_like(weighted)
          for index, (value, mode) in enumerate(zip(eigenvalues, vectors.T)):
              if abs(value) < tolerance:
                  continue                      # rigid-body and flat directions carry no step
              along = float(weighted @ mode)
              if index == selected_order:
                  # The softest non-trivial mode is the one to climb. Where its curvature is
                  # already negative a Newton step ascends by itself; where it is positive a
                  # Newton step would descend, so the walk pushes against the force by a fixed
                  # amount instead. That distinction is the whole of min-mode following, and
                  # without it the walk is just a minimizer that happens to start high: on the
                  # sharp-ridge potential the probe frames all sit past the inflection, where
                  # a pure Newton step returns the product basin every time.
                  displacement += (along / value if value < -tolerance
                                   else -np.sign(along) * config.min_mode_step_A) * mode
                  continue
              displacement += (along / value) * mode
          cartesian = displacement.reshape(-1, 3) / root
          extent = float(np.linalg.norm(cartesian, axis=1).max())
          if not np.isfinite(extent) or extent < 1e-12:
              return None, {"reason": "degenerate_step", "steps": step,
                          "history": history[-5:], **bounds}
          if extent > config.min_mode_step_A:
              cartesian *= config.min_mode_step_A / extent
          origin = atoms.positions.copy()
          accepted = False
          for _ in range(config.min_mode_backtracks + 1):
              atoms.set_positions(origin + cartesian)
              if ceiling is None:
                  accepted = True
                  break
              trial = float(guard.get_potential_energy(atoms))
              if trial <= ceiling:
                  accepted = True
                  break
              if bounds["ceiling_hit_step"] is None:
                  bounds["ceiling_hit_step"] = step
              bounds["backtracks"] += 1
              cartesian = cartesian * 0.5
              if float(np.linalg.norm(cartesian, axis=1).max()) < 1e-6:
                  break
          if not accepted:
              atoms.set_positions(origin)
              return None, {"reason": "ceiling_blocked", "steps": step,
                            "ceiling_eV": ceiling, **bounds,
                            "meaning": ("no step short enough to stay under the ceiling; the "
                                        "ridge this walk is on lies outside the region")}
          eigenvalues = None if (step + 1) % config.min_mode_refresh == 0 else eigenvalues
      else:
        return None, {"reason": "not_converged", "steps": config.min_mode_max_steps,
                      "history": ([history[0]] + history[-4:]) if history else [], **bounds}
    except GateRejected as exc:
        # The walk drove the structure somewhere the reliability gate refuses. That is a
        # stated outcome, not an exception for the search loop to handle: aiming at a stiff
        # mode pushes atoms into each other, and the guard is right to stop it.
        return None, {"reason": f"gate_{exc.code}", "steps": len(history), **bounds}

    eigenvalues, vectors, floor, provenance = curvature_spectrum(
        atoms, factory, config, source=config.hessian_source)
    negative = [float(v) for v in eigenvalues if v < -tolerance]
    worst_floor = floor_residual(floor)
    if worst_floor > config.trivial_floor_fraction * tolerance:
        # The same gate the minimum check applies. Without it a rigid-body residual can be
        # read as the one negative mode and a numerical artefact becomes an index-one
        # saddle.
        return None, {"reason": "trivial_mode_floor_above_tolerance", "steps": len(history),
                      "trivial_mode_floor_worst_residual": worst_floor, **bounds}
    final = float(atoms.get_potential_energy())
    if ceiling is not None and final > ceiling:
        # Backtracking bounds the path; the result has to be inside the region too, or the
        # bound could be walked around by taking smaller steps.
        return None, {"reason": "above_ceiling_at_saddle", "steps": len(history),
                      "energy_eV": final, "ceiling_eV": ceiling, **bounds}
    diagnostics = {"steps": len(history), "fmax_eV_A": history[-1]["fmax_eV_A"],
                   "followed_order": selected_order,
                   "final_mode_overlap": selected_overlap,
                   "saddle_order": len(negative), "negative_eigenvalues": negative,
                   "imaginary_wavenumbers_icm": _wavenumbers(negative),
                   "trivial_mode_floor_worst_residual": floor_residual(floor),
                   "provenance": provenance,
                   "history": ([history[0]] + history[-4:]) if history else [],
                   **bounds,
                   "meaning": "Newton walk in the analytic Hessian eigenbasis; index confirmed"}
    if len(negative) != 1:
        diagnostics["reason"] = "wrong_index"
        return None, diagnostics
    mode = vectors[:, 0].reshape(-1, 3) / root
    mode = mode / np.linalg.norm(mode)
    diagnostics["unstable_mode"] = mode.tolist()
    # A climbed saddle needs the same scale-resolved reading as one a quench stalled on,
    # or the significance test has nothing to work with and the acceptance criterion that
    # asks for it cannot be evaluated at all.
    if config.response_band_enabled:
        try:
            diagnostics["response"] = response_curvature(atoms, guard, mode, config)
        except GateRejected as exc:
            diagnostics["response"] = {"reason": f"gate_{exc.code}"}
    return snapshot(atoms, guard), diagnostics


def torsion_response(atoms, guard, torsion, step_rad, reference_energy=None):
    """Generalized force and stiffness of one torsion, in eV/rad and eV/rad^2.

    Measured by differencing the energy along the coordinate's own exact path, so there
    is no chain rule and no ambiguity about how the other coordinates were held: for a
    bridge bond the path is a rigid rotation, under which nothing else moves at all.
    Costs two energy evaluations when the reference energy is already known.
    """
    energies = []
    for delta in (-step_rad, 0.0, step_rad):
        if delta == 0.0 and reference_energy is not None:
            energies.append(float(reference_energy))
            continue
        probe, _ = step_torsion(atoms, torsion, delta, in_place=False)
        energies.append(float(guard.get_potential_energy(probe)))
    gradient_ = (energies[2] - energies[0]) / (2 * step_rad)
    curvature = (energies[2] - 2 * energies[1] + energies[0]) / step_rad ** 2
    return float(gradient_), float(curvature)


def polish_soft_modes(atoms, guard, config):
    """Minimize each soft torsion along its own path, leaving stiff coordinates alone.

    Three tiers, because one rule cannot cover them:

      stiff (k above the maximum)   the global force tolerance already pins it
      soft, positively curved       require |g/k| < tolerance, and Newton-step to it
      free (|k| below floor, g too) nothing to converge; every value is equivalent, so
                                    it is recorded rather than demanded
      flat but biased               |k| below the floor while g is not: the local second
                                    derivative is too small to be a Newton scale, but the
                                    coordinate is not stationary and work is available

    Negative curvature is left alone as well: it means the structure is not a minimum
    along that coordinate, which is the minimum check's decision to make, not this one.

    The fourth tier exists because the third one was wrong, and measurably so. On
    3-oxobutanal the run's product minimum was reported 9.97 meV above the same minimum
    reached with a tight tolerance, and the whole difference sat in one methyl torsion that
    this function had labelled "free": curvature 0.0041 eV/rad^2, below the 0.01 floor,
    while its gradient was -0.014 eV/rad and the coordinate was 0.71 rad from its minimum.
    A flat coordinate has no gradient; a small local curvature only means the Newton scale
    is unusable near an inflection of the rotor barrier. Because the "free" verdict also
    feeds `free_bonds_of`, the same misreading quotiented that coordinate out of microstate
    comparison. Flatness is now tested with both quantities, and the gradient bound is
    derived rather than chosen: at the smallest curvature this function is willing to trust,
    a gradient below `floor * tolerance` cannot move the coordinate past the tolerance.
    """
    from .chemistry import canonical_labels
    graph = encode(atoms, config.bond_scale, active=config.active_atoms)
    labels = canonical_labels(atoms.numbers, graph.edges, graph.index)
    genuine, rotors = rotatable_torsions(graph, labels)
    # Whether the geometry moved is a separate fact from whether every coordinate is now
    # within tolerance, and the caller needs both: a polish that moved the structure has
    # invalidated the force tolerance it was standing on, however happy its own criterion is.
    report = {"tolerance_rad": config.torsion_tolerance_rad, "modes": [],
              "converged": True, "moved": False, "energy_evaluations": 0}
    for torsion in genuine + rotors:
        entry = {"indices": list(torsion.indices), "symmetry_order": torsion.symmetry_order}
        polished = False
        reference = float(guard.get_potential_energy(atoms))
        gradient_, curvature = torsion_response(atoms, guard, torsion,
                                                config.soft_mode_step_rad, reference)
        report["energy_evaluations"] += 2
        entry.update(curvature_eV_rad2=curvature, initial_gradient_eV_rad=gradient_)
        if curvature > config.soft_mode_curvature_max_eV_rad2:
            entry["tier"] = "stiff"
        elif curvature < 0 and abs(curvature) >= config.soft_mode_curvature_floor_eV_rad2:
            entry["tier"] = "negative_curvature"
            entry["note"] = "not a minimum along this coordinate; left to the minimum check"
        elif abs(curvature) < config.soft_mode_curvature_floor_eV_rad2:
            bound = (config.soft_mode_curvature_floor_eV_rad2
                     * config.torsion_tolerance_rad)
            entry["free_gradient_bound_eV_rad"] = bound
            if abs(gradient_) <= bound:
                entry["tier"] = "free"
                entry["note"] = "flat coordinate; its value carries no information"
            else:
                entry["tier"] = "flat_biased"
                entry["note"] = ("curvature below the floor but the coordinate is not "
                                 "stationary; scanned along its own path instead of "
                                 "Newton-stepping on an unusable curvature")
                if config.soft_mode_scan_points >= 3:
                    # The coordinate's path is exact -- a rigid rotation about the bond --
                    # so a scan over its symmetry-reduced range is well defined and does not
                    # depend on a local model that has just been shown not to hold. One
                    # range, because values a symmetry period apart are the same structure.
                    span = 2 * np.pi / max(int(torsion.symmetry_order), 1)
                    offsets = np.linspace(-span / 2, span / 2,
                                          config.soft_mode_scan_points + 1)[:-1]
                    scanned = []
                    for delta in offsets:
                        probe, _ = step_torsion(atoms, torsion, float(delta),
                                                in_place=False)
                        scanned.append(float(guard.get_potential_energy(probe)))
                    report["energy_evaluations"] += len(offsets)
                    best = int(np.argmin(scanned))
                    entry["scan_span_rad"] = float(span)
                    entry["scan_points"] = int(len(offsets))
                    entry["scan_offset_rad"] = float(offsets[best])
                    entry["scan_energy_gain_eV"] = float(reference - scanned[best])
                    if abs(offsets[best]) > 1e-12:
                        step_torsion(atoms, torsion, float(offsets[best]))
                        report["moved"] = True
                        gradient_, curvature = torsion_response(
                            atoms, guard, torsion, config.soft_mode_step_rad)
                        report["energy_evaluations"] += 3
                        entry["scanned_gradient_eV_rad"] = gradient_
                        entry["scanned_curvature_eV_rad2"] = curvature
                polished = True
        else:
            entry["tier"] = "soft"
            polished = True
        if polished:
            # A torsion is periodic, so no step ever needs to exceed half of its own
            # symmetry-reduced period: every value is reachable inside that, and a longer
            # step only wraps past a symmetry copy of where it started. That is where the
            # cap comes from -- the coordinate's periodicity, not a tuned trust radius.
            cap = np.pi / max(int(torsion.symmetry_order), 1)
            for _ in range(config.soft_polish_steps):
                # Below the floor the Newton scale is not trustworthy either, so the step is
                # taken against the floor rather than against a curvature that would produce
                # a meaningless offset -- measured, -0.014/0.0041 asked for 3.4 radians and
                # -0.15/0.01 on the test rotor asked for 15.
                scale = max(curvature, config.soft_mode_curvature_floor_eV_rad2)
                offset = float(np.clip(-gradient_ / scale, -cap, cap))
                entry["last_offset_rad"] = offset
                if abs(offset) < config.torsion_tolerance_rad:
                    break
                step_torsion(atoms, torsion, offset)
                report["moved"] = True
                report["energy_evaluations"] += 1
                gradient_, curvature = torsion_response(
                    atoms, guard, torsion, config.soft_mode_step_rad)
                report["energy_evaluations"] += 3
                if curvature <= 0:
                    entry["note"] = "curvature changed sign while polishing; stopped"
                    break
            entry["final_gradient_eV_rad"] = gradient_
            entry["final_curvature_eV_rad2"] = curvature
            scale = max(curvature, config.soft_mode_curvature_floor_eV_rad2)
            entry["implied_offset_rad"] = (abs(gradient_ / scale) if scale > 0 else None)
            entry["within_tolerance"] = bool(
                scale > 0 and abs(gradient_ / scale) < config.torsion_tolerance_rad)
            if not entry["within_tolerance"]:
                report["converged"] = False
        report["modes"].append(entry)
    return report


def quench(atoms, guard, config, observe=None):
    """Force-tolerance minimization, then a mode-resolved polish of the soft coordinates.

    The global tolerance stays where it is. Polishing moves the structure slightly, so
    the stiff coordinates are re-settled and the soft ones re-checked, for a bounded
    number of rounds.
    """
    total_steps = 0
    diagnostics = {"fmax_eV_A": config.quench_fmax_eV_A, "rounds": []}

    def minimize():
        nonlocal total_steps
        optimizer = FIRE(atoms, logfile=None)
        if observe is not None:
            optimizer.attach(lambda: observe(total_steps + optimizer.nsteps), interval=1)
        converged = optimizer.run(fmax=config.quench_fmax_eV_A, steps=config.quench_steps)
        total_steps += optimizer.nsteps
        return bool(converged)

    if not minimize():
        diagnostics.update(converged=False, reason="fmax_not_reached", steps=total_steps)
        return False, diagnostics

    for _ in range(config.soft_polish_rounds):
        report = polish_soft_modes(atoms, guard, config)
        diagnostics["rounds"].append(report)
        # Convergence needs the polish to be satisfied AND to have left the geometry alone.
        # Returning on "satisfied" alone reported a converged quench on a structure the
        # polish had just displaced: measured on the 3-oxobutanal product, a 0.35 rad methyl
        # scan step returned converged=True at fmax = 0.179 against a 0.02 tolerance. The old
        # polish only ever took steps at the tolerance scale, so the same hole was there and
        # small; it is not a new one, and it is not tolerable now that a step can be large.
        if report["converged"] and not report["moved"]:
            diagnostics.update(converged=True, steps=total_steps)
            return True, diagnostics
        if not minimize():
            diagnostics.update(converged=False, reason="fmax_not_reached_after_polish",
                               steps=total_steps)
            return False, diagnostics
    report = polish_soft_modes(atoms, guard, config)
    diagnostics["rounds"].append(report)
    settled = bool(report["converged"] and not report["moved"])
    diagnostics.update(converged=settled, steps=total_steps,
                       reason=None if settled else
                              ("soft_mode_moved_on_the_last_round" if report["converged"]
                               else "soft_mode_not_converged"))
    return settled, diagnostics


def follow_unstable_mode(saddle, mode, factory, config, sign, recorder=None):
    """Step off a first-order saddle along its unstable mode and quench.

    One negative Hessian eigenvalue means exactly one direction leads downhill on both
    sides, and its eigenvector is already in hand from the minimum check. Two quenches
    therefore name the two minima the saddle connects, which is the difference between
    "a perturbation reached B" and "A and B are joined through this saddle".

    This is a steepest-descent endpoint, not an integrated reaction path: it identifies
    the endpoints, and says nothing about the path between them or about any barrier.

    The step stays small on purpose. A large displacement can cross more than one dividing
    surface and land in a basin that is not adjacent to the saddle, and then the endpoint it
    names is not the one this saddle connects. Where the descent stops on a further saddle
    the caller relays from there rather than pushing harder from here.
    """
    mode = np.asarray(mode, dtype=float).reshape(-1, 3)
    extent = np.linalg.norm(mode, axis=1).max()
    if extent < 1e-12:
        raise ValueError("Unstable mode has no displacement")
    atoms = saddle.copy()
    guard = GuardedCalculator(factory(), config, saddle.get_potential_energy())
    atoms.calc = guard
    atoms.set_positions(saddle.positions + sign * config.irc_step_A / extent * mode)
    atoms.set_momenta(np.zeros((len(atoms), 3)))
    phase = "descent_plus" if sign > 0 else "descent_minus"

    def watch(step):
        if recorder is not None:
            recorder.observe(phase, step, atoms)

    try:
        watch(0)
        converged, diagnostics = quench(atoms, guard, config, watch)
    except GateRejected as exc:
        return None, {"reason": f"gate_{exc.code}"}
    if not converged:
        return None, {"reason": diagnostics.get("reason") or "quench_nonconverged"}
    atoms.set_momenta(np.zeros((len(atoms), 3)))
    endpoint = snapshot(atoms, guard)
    # The quench's verdict on which coordinates are flat travels with the geometry it was
    # measured on, the way relax_source already does it.
    endpoint.info["quench"] = diagnostics
    return endpoint, diagnostics


def _arrived(decision):
    """Whether a rung's named verdict counts as reaching somewhere other than the source.

    Decided here rather than by the naming callback, because the callback is run
    speculatively -- on both descents of a stalled rung -- and must not be able to end the
    walk or hand anything back. An explicit "stop" is still honoured so a caller can be
    stricter than this default, never looser.
    """
    if not decision:
        return False
    if "stop" in decision:
        return bool(decision["stop"])
    return decision.get("verdict") not in (None, "returned_to_source",
                                           "minimum_unclassified")


def _projected_direction_from(displacement, masses):
    """Mass-weighted unit direction with rigid-body translation removed.

    Rotation is left in: the two structures are one small step apart, so any rotation
    between them is part of how the walk actually moved rather than a symmetry to quotient.
    """
    weighted = (np.sqrt(masses[:, None])
                * np.asarray(displacement, dtype=float).reshape(-1, 3))
    weighted = weighted - weighted.mean(axis=0)
    norm = np.linalg.norm(weighted)
    return None if norm < 1e-12 else weighted / norm


def align_mode_with_travel(mode, travel, masses, sign):
    """Orient a relay saddle's unstable mode so the walk keeps going the way it came.

    An eigenvector's sign is arbitrary -- `numpy.linalg.eigh` picks one, and which one is
    not a property of the physics. The relay used to reapply the same numeric `sign` to
    whatever `confirm_minimum` handed back, so after a hop the descent could reverse and
    walk back toward the saddle it had just left. The two-sided descent would then report
    both sides landing in the same basin for a reason that has nothing to do with the
    potential.

    The direction the walk was travelling is the displacement that brought it here. The
    next step will be along `sign * mode`, so require

        (sign * mode) . travel > 0

    and flip the mode when it is not. The overlap is measured in mass-weighted coordinates,
    the same metric the min-mode walk uses for target tracking, so the two agree about what
    "the same direction" means.

    Returns (mode, overlap); overlap is None when the travel has no internal component, in
    which case nothing is claimed and the mode is returned untouched.
    """
    basis_travel = _projected_direction_from(travel, masses)
    if basis_travel is None:
        return np.asarray(mode, dtype=float), None
    weighted = np.sqrt(masses[:, None]) * np.asarray(mode, dtype=float).reshape(-1, 3)
    norm = np.linalg.norm(weighted)
    if norm < 1e-12:
        return np.asarray(mode, dtype=float), None
    overlap = float(np.sum(basis_travel * (weighted / norm)))
    mode = np.asarray(mode, dtype=float)
    return (mode if sign * overlap >= 0 else -mode), overlap


def descend_saddle(saddle, mode, factory, config, sign, seed=0, recorder=None):
    """One side of a saddle, followed downhill until it reaches a minimum or gives up.

    Where the descent stops on a further first-order saddle it continues from THAT saddle
    along its own unstable mode rather than pushing harder from the first one: a larger
    displacement can cross several dividing surfaces and name a basin the saddle does not
    connect. On 3-oxobutanal the soft torsions put several ridges within one step of each
    other, so relaying is the difference between naming an endpoint and reporting failure.

    Every hop is recorded, because a chain of them is weaker connectivity evidence than a
    single step.

    Returns (endpoint, checked, relays); endpoint is None if no minimum was reached.
    """
    current, current_mode, relays = saddle, mode, []
    endpoint, checked = None, {}
    for hop in range(config.irc_relay_hops + 1):
        endpoint, diagnostics = follow_unstable_mode(current, current_mode, factory, config,
                                                    sign, recorder=recorder)
        if endpoint is None:
            return None, {"reason": (diagnostics or {}).get("reason")}, relays
        settled, checked = confirm_minimum(endpoint, factory, config, seed + hop)
        if settled:
            return endpoint, checked, relays
        if checked.get("saddle_order") != 1:
            return None, checked, relays
        relayed = checked.get("unstable_mode")
        if relayed is None:
            return None, checked, relays
        # The new eigenvector's sign is arbitrary; orient it along the way the walk came.
        travel = endpoint.positions - current.positions
        relayed, overlap = align_mode_with_travel(relayed, travel,
                                                  endpoint.get_masses(), sign)
        relays.append({"hop": hop, "saddle_order": checked.get("saddle_order"),
                       "icm": checked.get("imaginary_wavenumbers_icm"),
                       "travel_overlap": overlap,
                       "overlap_meaning": ("mass-weighted overlap between the relayed mode "
                                           "and the displacement that reached this saddle; "
                                           "null means the travel had no internal "
                                           "component and no orientation was claimed")})
        current, current_mode = endpoint, relayed
    return None, checked, relays


def arrival_through_saddle(saddle, mode, factory, config, verdict, recorder=None,
                           seed=0):
    """Does a first-order saddle's neighbourhood belong to one chemical state?

    A quench stopping on a saddle does not by itself say the geometry failed to arrive. A
    first-order saddle lies on the boundary between exactly two basins, so if both of its
    descents name the same state then every neighbourhood of it except the ridge itself
    belongs to that state -- and a structure that quenched onto the ridge came from one side
    or the other of it.

    That is the situation the 3-oxobutanal benchmark was stuck in. Every rung of every
    boundary continuation in round four quenched onto a rotor saddle 11 meV above the
    product minimum, at -65 to -79 cm^-1, while the proton had demonstrably crossed
    (q_PT +0.70 -> -0.69). Nine trials did the same. Nothing was ever a confirmed minimum,
    so nothing was ever admitted, and the directed edge stayed missing while the evidence
    for it sat in the record unread.

    Both sides, small steps, one calculator and one set of settings. Disagreement is the
    interesting case and is not smoothed over: it means the geometry really is on a reaction
    ridge, and the caller keeps its "stopped on a saddle" verdict.

    Returns (agreed, sides). `agreed` is the side whose named verdict both descents produced,
    or None; `sides` is the per-side record either way.
    """
    sides = []
    for sign in (1, -1):
        landed, checked, relays = descend_saddle(saddle, np.asarray(mode), factory, config,
                                                sign, seed=seed + 7919 * sign,
                                                recorder=recorder)
        if landed is None:
            sides.append({"sign": sign, "reason": checked.get("reason"),
                          "saddle_order": checked.get("saddle_order"),
                          "relays": relays})
            continue
        named = verdict(landed, True, checked)
        sides.append({"sign": sign, "confirmed_minimum": True,
                      "saddle_order": checked.get("saddle_order"),
                      "relays": relays,
                      "verdict": (named or {}).get("verdict"),
                      "broken": (named or {}).get("broken"),
                      "formed": (named or {}).get("formed"),
                      "energy_eV": float(landed.get_potential_energy()),
                      "named": named, "structure": landed})
    complete = [side for side in sides if side.get("named")]
    agreed = None
    if len(complete) == 2 and all(complete[0]["named"][field] == complete[1]["named"][field]
                                 for field in ("verdict", "broken", "formed")):
        agreed = complete[0]
    return agreed, sides


def continue_across_boundary(frame, tangent, factory, config, verdict=None,
                             recorder=None):
    """Keep delivering the probe past the frame where the topology switched.

    A first-order saddle with a two-sided descent establishes that two minima are joined
    through it. It does not establish that the probe which was delivered gets from the
    source to the product, and those are two different claims: on 3-oxobutanal round three
    the proton-transfer channel had three saddles supporting it and not one trial had
    observed A reaching B. The delivered path stopped at the dividing surface and the
    record could only say that it stopped.

    So this continues the same motion. The tangent is the probe's own signed local
    direction -- the difference across the frame where the topology switched, or the
    coordinate's gradient when no adjacent frame was kept -- and the walk advances along it
    in rungs whose largest single-atom displacement is capped. Each rung is quenched
    independently from its displaced geometry, not continued from the previous rung's
    endpoint, so a rung's verdict is about that rung's position and nothing else.

    Four verdicts, and they are the point of the exercise:

      returned to the source    the crossing had not yet entered the product's basin
      reached another state     a directed observation of the delivered path arriving
      stopped on a saddle       the delivered path runs along the dividing surface
      left the domain           the rung is outside what the potential describes, so the
                                step or the amplitude was too large

    This is emphatically not the two-sided descent wearing a different hat. The descent
    starts at a saddle and goes downhill along an eigenvector; this starts at a probe frame
    and goes forward along the probe's tangent. Presenting the first as an observation of
    the delivered path is the fabrication this function exists to avoid.

    `verdict` is a callback taking (endpoint, confirmed, diagnostics) and returning a dict
    with at least a "verdict" key; the caller owns state identity, so the naming of basins
    does not happen here. It must be pure: it is run speculatively on both descents of a
    stalled rung, including on rungs that turn out not to have arrived, so a callback that
    recorded anything would record those too. An explicit "stop" is honoured; otherwise any
    verdict other than "returned_to_source" ends the walk.

    Returns (endpoint, report). The endpoint is None unless a rung arrived somewhere other
    than the source; nothing is admitted, recorded or claimed by this function.
    """
    report = {"step_A": config.boundary_continuation_step_A,
              "max_steps": config.boundary_continuation_max_steps,
              "rungs": [], "outcome": None,
              "meaning": ("continuation of the delivered probe past the topology switch; "
                          "each rung quenched independently, no path or barrier implied")}
    masses = frame.get_masses()
    projected = _projected_direction(frame.positions, masses, tangent)
    if projected is None:
        report["outcome"] = "tangent_has_no_internal_component"
        return None, report
    # Back to Cartesian and scaled so the cap is on the largest single-atom displacement.
    # A cap on the norm would let a twelve-atom molecule move every atom by the full step.
    direction = projected / np.sqrt(masses[:, None])
    extent = float(np.linalg.norm(direction, axis=1).max())
    if not np.isfinite(extent) or extent < 1e-12:
        report["outcome"] = "degenerate_direction"
        return None, report
    direction = direction / extent
    origin = frame.positions.copy()
    for rung in range(1, config.boundary_continuation_max_steps + 1):
        travelled = rung * config.boundary_continuation_step_A
        entry = {"rung": rung, "displacement_A": travelled}
        advanced = frame.copy()
        guard = GuardedCalculator(factory(), config, None)
        advanced.calc = guard
        advanced.set_positions(origin + travelled * direction)
        advanced.set_momenta(np.zeros((len(advanced), 3)))
        try:
            entry["energy_eV"] = float(guard.get_potential_energy(advanced))
            if recorder is not None:
                recorder.observe("continuation", rung, advanced,
                                 energy=entry["energy_eV"])
            quenched = advanced.copy()
            quenched.calc = guard
            converged, diagnostics = quench(quenched, guard, config,
                                            (lambda step, r=rung:
                                             recorder.observe(f"continuation_quench{r}",
                                                              step, quenched))
                                            if recorder is not None else None)
        except GateRejected as exc:
            entry.update(verdict="left_domain", gate=exc.code,
                         note=("the potential's reliability gate refused this rung, so the "
                               "step or the amplitude carried the structure out of the "
                               "region the model describes"))
            report["rungs"].append(entry)
            report["outcome"] = "left_domain"
            return None, report
        entry["quench_converged"] = bool(converged)
        if not converged:
            entry.update(verdict="quench_failed",
                         reason=diagnostics.get("reason") or "quench_nonconverged")
            report["rungs"].append(entry)
            continue
        quenched.set_momenta(np.zeros((len(quenched), 3)))
        endpoint = snapshot(quenched, guard)
        endpoint.info["quench"] = diagnostics
        entry["endpoint_energy_eV"] = float(endpoint.get_potential_energy())
        confirmed, checked = confirm_minimum(endpoint, factory, config,
                                            config.seed + 104729 * rung)
        entry["saddle_order"] = checked.get("saddle_order")
        if not confirmed:
            entry.update(verdict="stopped_on_saddle" if checked.get("saddle_order") == 1
                                 else "not_a_minimum",
                         imaginary_wavenumbers_icm=checked.get("imaginary_wavenumbers_icm"))
            mode = checked.get("unstable_mode")
            agreed = None
            if (config.boundary_continuation_descend_saddles and verdict is not None
                    and checked.get("saddle_order") == 1 and mode is not None):
                agreed, sides = arrival_through_saddle(
                    endpoint, mode, factory, config, verdict, recorder=recorder,
                    seed=config.seed + 15485863 * rung)
                entry["descents"] = [{k: v for k, v in side.items()
                                      if k not in ("named", "structure")}
                                     for side in sides]
                entry["descent_meaning"] = (
                    "both sides of the saddle this rung quenched onto; agreement means the "
                    "rung is inside that state's basin, disagreement means it is on a ridge")
            if agreed is None:
                entry["note"] = ("this rung quenched onto a stationary point that is not a "
                                 "minimum and its two descents did not agree on one state, "
                                 "so it is not evidence of arrival")
                report["rungs"].append(entry)
                continue
            decision = dict(agreed["named"])
            entry.update({k: v for k, v in decision.items() if k != "stop"})
            entry["arrival"] = "via_saddle"
            entry["arrival_meaning"] = (
                "the rung quenched onto a saddle whose both descents name this state, which "
                "is weaker evidence than quenching into the minimum directly")
            entry["saddle_icm"] = checked.get("imaginary_wavenumbers_icm")
            report["rungs"].append(entry)
            if _arrived(decision):
                report["outcome"] = entry["verdict"]
                return agreed["structure"], report
            continue
        decision = verdict(endpoint, confirmed, checked) if verdict is not None else {}
        entry.update({k: v for k, v in decision.items() if k != "stop"})
        entry.setdefault("verdict", "minimum_unclassified")
        report["rungs"].append(entry)
        if _arrived(decision):
            report["outcome"] = entry["verdict"]
            return endpoint, report
    report["outcome"] = "exhausted"
    return None, report


def gate_collateral(delivered, config, atoms):
    """Refuse a probe that damaged coordinates it was not aiming at.

    Reaching the requested amplitude proves only that the named coordinate moved. The
    gated quantities are the ones that mean the structure was torn: a stretched covalent
    bond, an atom clash, a changed fragment count. Coupled angles and torsions are
    recorded but not gated, because a real bend probe has to move the angles that share
    its vertex.
    """
    if not config.collateral_gate:
        return
    collateral = delivered.get("collateral")
    if collateral is None:
        return
    if collateral["max_bond_change_A"] > config.collateral_bond_A:
        raise GateRejected(
            "collateral_bond",
            f"probe changed bond {collateral['worst_bond']} by "
            f"{collateral['max_bond_change_A']:.3f} A, above "
            f"{config.collateral_bond_A} A", atoms)
    if collateral["closest_distance_ratio"] < config.min_distance_ratio:
        raise GateRejected(
            "collateral_clash",
            f"probe brought {collateral['closest_pair']} to "
            f"{collateral['closest_distance_ratio']:.3f} of their covalent sum", atoms)
    if collateral["fragments_after"] != collateral["fragments_before"]:
        raise GateRejected(
            "collateral_fragmentation",
            f"probe changed the fragment count from {collateral['fragments_before']} "
            f"to {collateral['fragments_after']}", atoms)


def confirm_minimum(candidate, factory, config, seed):
    """Confirm a converged quench really sits in a local minimum.

    A vanishing gradient alone admits maxima and saddles, so a barrier top would be
    recorded as a product basin. Two methods are available and the choice is recorded:

    "hessian" diagonalizes a central-difference mass-weighted Hessian with rigid-body
    modes projected out and requires every remaining eigenvalue to be nonnegative. It
    costs 6N force evaluations and is decisive.

    "probe" applies seeded mass-weighted displacements and requenches, which is cheap
    but only samples directions. Measured on ethanol torsional barrier tops it caught
    one of two equivalent saddles, so it is a fallback for systems too large to
    diagonalize, not a default.

    Either way the outcome fails closed and the evidence is kept.
    """
    if config.minimum_check == "none":
        return True, {"method": "none",
                      "meaning": "minimum confirmation disabled; basins are unverified"}
    if config.minimum_check == "hessian":
        try:
            eigenvalues, vectors, floor, provenance = curvature_spectrum(
                candidate, factory, config, source=config.hessian_source)
        except GateRejected as exc:
            return False, {"method": "hessian", "reason": f"gate_{exc.code}",
                           "saddle_order": None}
        negative = [float(v) for v in eigenvalues
                    if v < -config.minimum_check_eigenvalue_tol]
        worst_floor = floor_residual(floor)
        diagnostics = {
            "method": "hessian",
            "provenance": provenance,
            "eigenvalue_tol": config.minimum_check_eigenvalue_tol,
            "lowest_eigenvalues": [float(v) for v in eigenvalues[:6]],
            "negative_eigenvalues": negative,
            "saddle_order": len(negative),
            "trivial_mode_floor": floor,
            "trivial_mode_floor_worst_residual": worst_floor,
            "imaginary_wavenumbers_icm": _wavenumbers(negative),
            "units": "eV/(Angstrom^2 amu)",
            "meaning": "rigid-body modes projected out; negative curvature means saddle",
        }
        # The rigid-body modes bound the error of this whole calculation once the part
        # the gradient explains is removed, so a residual that reaches the decision
        # threshold makes the sign structure unreadable. Failing closed here keeps a noisy
        # Hessian from silently minting saddles.
        if worst_floor > config.trivial_floor_fraction * config.minimum_check_eigenvalue_tol:
            diagnostics["reason"] = "trivial_mode_floor_above_tolerance"
            diagnostics["saddle_order"] = None
            return False, diagnostics
        if len(negative) == 1:
            # Exactly one negative mode is a first-order saddle: a transition state
            # candidate, not refuse. Its eigenvector is the reaction coordinate and is
            # what a dimer or NEB refinement would otherwise have to rediscover.
            masses = candidate.get_masses()
            mode = vectors[:, 0].reshape(-1, 3) / np.sqrt(masses[:, None])
            diagnostics["unstable_mode"] = (mode / np.linalg.norm(mode)).tolist()
            if config.response_band_enabled:
                # With the reference energy supplied, a band probe that blows up is
                # rejected by the same gate every other evaluation goes through.
                guard = GuardedCalculator(factory(), config,
                                          candidate.get_potential_energy())
                probe = candidate.copy()
                probe.calc = guard
                try:
                    diagnostics["response"] = response_curvature(
                        probe, guard, np.asarray(diagnostics["unstable_mode"]), config)
                except GateRejected as exc:
                    diagnostics["response"] = {"reason": f"gate_{exc.code}"}
        return not negative, diagnostics
    return _probe_minimum(candidate, factory, config, seed)


def _probe_minimum(candidate, factory, config, seed):
    """Sampled-direction curvature probe. Separate a local minimum from any other
    stationary point without building a Hessian.

    A converged quench only proves the gradient vanished, so a barrier top or saddle is
    accepted as readily as a basin. Each direction is measured rather than followed:
    two force evaluations give v^T H v, whose sign says directly whether the structure is
    a minimum along it. That replaces displace-and-requench, which cost a full quench per
    direction and returned only "came back" or "did not" -- on the ethanol torsional
    barrier tops it found one of the two equivalent saddles, because sliding away is a
    weaker signal than a negative number.

    Directions come from the molecule's own rotatable torsions first, since that is where
    the soft and unstable modes live, and are then filled out with seeded random ones.
    Rigid-body components are projected out of every direction.

    This is still a sampled probe, not a spectrum: it cannot see an unstable mode
    orthogonal to every direction tried, and it gives no eigenvalue. Use it where the
    system is too large to diagonalize. Returns (confirmed, diagnostics).
    """
    from .chemistry import canonical_labels
    masses = candidate.get_masses()
    tolerance = config.minimum_check_eigenvalue_tol
    guard = GuardedCalculator(factory(), config, candidate.get_potential_energy())
    probe = candidate.copy()
    probe.calc = guard

    directions = []
    graph = encode(candidate, config.bond_scale, active=config.active_atoms)
    labels = canonical_labels(candidate.numbers, graph.edges, graph.index)
    genuine, rotors = rotatable_torsions(graph, labels)
    delta = config.soft_mode_step_rad
    for torsion in genuine + rotors:
        plus, _ = step_torsion(candidate, torsion, delta, in_place=False)
        minus, _ = step_torsion(candidate, torsion, -delta, in_place=False)
        tangent = _projected_direction(candidate.positions, masses,
                                       (plus.positions - minus.positions) / (2 * delta))
        if tangent is not None:
            directions.append(({"source": "torsion",
                                "indices": list(torsion.indices)}, tangent))

    rng = np.random.default_rng(seed)
    attempts = 0
    while len(directions) < len(genuine) + len(rotors) + config.minimum_check_probes:
        attempts += 1
        if attempts > 100:
            raise ValueError("Could not draw an independent minimum-confirmation direction")
        drawn = _projected_direction(candidate.positions, masses,
                                     rng.normal(size=(len(candidate), 3)))
        if drawn is None:
            continue
        directions.append(({"source": "random"}, drawn))

    probes = []
    for entry, direction in directions:
        try:
            curvature = _curvature_along(probe, guard, direction,
                                        config.minimum_check_step_A)
        except GateRejected as exc:
            entry.update(curvature=None, negative=None, reason=f"gate_{exc.code}")
            probes.append(entry)
            continue
        entry.update(curvature_eV_A2_amu=curvature, negative=bool(curvature < -tolerance))
        probes.append(entry)

    measured = [entry for entry in probes if entry.get("curvature_eV_A2_amu") is not None]
    negative = [entry for entry in measured if entry["negative"]]
    confirmed = bool(measured) and not negative
    return confirmed, {
        "method": "probe",
        "probes": probes,
        "step_A": config.minimum_check_step_A,
        "eigenvalue_tol": tolerance,
        "negative_directions": len(negative),
        "lowest_curvature": min((e["curvature_eV_A2_amu"] for e in measured), default=None),
        "units": "eV/(Angstrom^2 amu)",
        "meaning": ("sampled curvature probe; a negative value is decisive evidence of a "
                    "saddle, an all-positive result is not proof of a minimum"),
    }


def relax_source(atoms, factory, config):
    validate_atoms(atoms, config.active_atoms)
    relaxed = atoms.copy()
    relaxed.set_momenta(np.zeros((len(relaxed), 3)))
    guard = GuardedCalculator(factory(), config)
    relaxed.calc = guard
    converged, diagnostics = quench(relaxed, guard, config)
    if not converged:
        raise ValueError(f"Initial structure failed to quench: {diagnostics}")
    relaxed.set_momenta(np.zeros((len(relaxed), 3)))
    source = snapshot(relaxed, guard)
    source.info["quench"] = diagnostics
    confirmed, diagnostics = confirm_minimum(source, factory, config, config.seed)
    if not confirmed:
        raise ValueError(f"Source structure is not a confirmed local minimum: {diagnostics}")
    return source


def run_trial(source, probe, factory, config, directory, trial_id):
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=False)
    atoms = source.copy()
    physical = GuardedCalculator(factory(), config, source.get_potential_energy())
    atoms.calc = physical
    record = {
        "attempt_id": trial_id, "probe": probe.to_dict(), "status": "running",
        "classification": None, "fingerprint": None, "failure": None,
        "trajectory": f"trials/{trial_id}/response.traj",
        "observations": f"trials/{trial_id}/observations.jsonl",
        "endpoint": None, "uncertainty_status": "unavailable",
        "max_uncertainty_eV_A": None, "max_nve_drift_eV_atom": 0.0,
        "pulse_work_eV": 0.0, "tail_persistence_steps": 0,
    }
    previous_edges = None
    free_graphs = []
    elapsed_fs = 0.0
    sample_count = 0

    with Trajectory(str(directory / "response.traj"), "w") as trajectory:
        def observe(phase, step, force_save=False, free_baseline=None):
            nonlocal previous_edges, sample_count
            frame = snapshot(atoms, physical)
            graph = encode(atoms, config.bond_scale, previous_edges, config.bond_hysteresis,
                           active=config.active_atoms)
            events = {"formed": [], "broken": []}
            if previous_edges is not None:
                events = {"formed": [list(e) for e in sorted(graph.edges - previous_edges)],
                          "broken": [list(e) for e in sorted(previous_edges - graph.edges)]}
            changed = any(events.values())
            previous_edges = graph.edges
            energy = frame.get_potential_energy()
            kinetic = atoms.get_kinetic_energy()
            if not np.isfinite(kinetic):
                raise GateRejected("nonfinite", "Non-finite kinetic energy", atoms)
            if free_baseline is not None:
                drift = abs(energy + kinetic - free_baseline) / len(atoms)
                record["max_nve_drift_eV_atom"] = max(record["max_nve_drift_eV_atom"], drift)
                if drift > config.max_nve_drift_eV_atom:
                    raise GateRejected("nve_drift", "Unbiased NVE energy drift exceeds limit", atoms)
            uncertainty = physical.results.get("force_uncertainty_eV_A")
            if uncertainty is not None:
                record["uncertainty_status"] = "available"
                record["max_uncertainty_eV_A"] = max(record["max_uncertainty_eV_A"] or 0, uncertainty)
            if phase == "free":
                free_graphs.append(encode(atoms, config.bond_scale,
                                          active=config.active_atoms).edges)
            if force_save or changed or step % config.sample_interval == 0:
                i, j = np.triu_indices(len(graph.index), k=1)
                com = atoms.get_center_of_mass()
                momenta = atoms.get_momenta()
                bias = (atoms.calc.results.get("bias_energy", 0.0)
                        if isinstance(atoms.calc, PairPulse) else 0.0)
                data = {
                    "frame": sample_count, "phase": phase, "step": step,
                    "time_fs": elapsed_fs, "potential_eV": energy, "kinetic_eV": kinetic,
                    "physical_total_eV": energy + kinetic, "bias_energy_eV": bias,
                    "extended_total_eV": energy + kinetic + bias,
                    "max_force_eV_A": float(np.linalg.norm(frame.get_forces(), axis=1).max()),
                    "uncertainty_eV_A": uncertainty, "coordination": graph.coordination.tolist(),
                    "pair_distances_A": graph.distances[i, j].tolist(),
                    "edges": [list(e) for e in sorted(graph.edges)], "events": events,
                    # Named coordinates, because an unlabelled upper-triangle vector of
                    # every pair distance held the information and answered no question.
                    "tracked": named_values(atoms.positions, config.tracked_coordinates),
                    "active_subset": graph.index.tolist(),
                    "com_A": com.tolist(), "total_momentum_ase": momenta.sum(axis=0).tolist(),
                    "angular_momentum_ase": np.cross(atoms.positions - com, momenta).sum(axis=0).tolist(),
                }
                frame.info.update({"phase": phase, "time_fs": elapsed_fs})
                trajectory.write(frame)
                append_json(directory / "observations.jsonl", data)
                sample_count += 1
            return energy + kinetic

        try:
            atoms.set_momenta(np.zeros((len(atoms), 3)))
            observe("source", 0, True)
            if probe.mode == "displace":
                # Geometry probes act on a quenched, zero-momentum structure.
                record["probe_delivered"] = apply_probe(atoms, probe, config.bond_scale,
                                                       config.active_atoms)
                gate_collateral(record["probe_delivered"], config, atoms)
            thermal_momenta(atoms, config.temperature_K, probe.seed)
            before_kick = atoms.get_kinetic_energy()
            if probe.mode == "kick":
                record["probe_delivered"] = apply_probe(atoms, probe)
            record["injected_kinetic_eV"] = atoms.get_kinetic_energy() - before_kick
            observe("perturbed", 0, True)
            record["geometry_energy_change_eV"] = (
                atoms.get_potential_energy() - source.get_potential_energy())
            if probe.family == "pulse":
                _, r0 = pair_axis(atoms, probe.pair)
                atoms.calc = PairPulse(physical, probe.pair, probe.sign, probe.amplitude, r0)
                dynamics = VelocityVerlet(atoms, timestep=config.timestep_fs * units.fs, logfile=None)
                pulse_initial = physical.get_potential_energy(atoms) + atoms.get_kinetic_energy()
                for step in range(1, config.pulse_steps + 1):
                    dynamics.run(1)
                    elapsed_fs += config.timestep_fs
                    observe("forced", step, step == config.pulse_steps)
                _, r1 = pair_axis(atoms, probe.pair)
                work = probe.sign * probe.amplitude * (r1 - r0)
                record["pulse_work_eV"] = work
                record["pulse_work_residual_eV"] = (
                    physical.get_potential_energy(atoms) + atoms.get_kinetic_energy() - pulse_initial - work)
                # New calculator/integrator: no stale pulse forces at the switch.
                atoms.calc = physical
            free_baseline = observe("free_start", 0, True)
            dynamics = VelocityVerlet(atoms, timestep=config.timestep_fs * units.fs, logfile=None)
            for step in range(1, config.response_steps + 1):
                dynamics.run(1)
                elapsed_fs += config.timestep_fs
                observe("free", step, step == config.response_steps, free_baseline)
            raw = snapshot(atoms, physical)
            write(str(directory / "raw_endpoint.extxyz"), raw)
            record["raw_endpoint"] = f"trials/{trial_id}/raw_endpoint.extxyz"
            atoms.set_momenta(np.zeros((len(atoms), 3)))
            converged, quench_diagnostics = quench(
                atoms, physical, config, lambda step: observe("quench", step))
            atoms.set_momenta(np.zeros((len(atoms), 3)))
            observe("quenched", quench_diagnostics["steps"], True)
            record["quench_converged"] = bool(converged)
            record["quench"] = quench_diagnostics
            if not converged:
                record.update(status="quench_failed",
                              failure={"code": "quench_nonconverged",
                                       "message": quench_diagnostics.get("reason")
                                                  or "quench did not converge"})
                return Outcome(record)
            endpoint = snapshot(atoms, physical)
            # The quench's own verdict on which coordinates are flat travels with the
            # geometry it was measured on, the way relax_source already does it.
            endpoint.info["quench"] = quench_diagnostics
            write(str(directory / "endpoint.extxyz"), endpoint)
            record["endpoint"] = f"trials/{trial_id}/endpoint.extxyz"
            label, fingerprint = classify(source, endpoint, config)
            endpoint_edges = encode(endpoint, config.bond_scale,
                                    active=config.active_atoms).edges
            for edges in reversed(free_graphs):
                if edges != endpoint_edges:
                    break
                record["tail_persistence_steps"] += 1
            record.update(status="completed", classification=label, fingerprint=fingerprint)
            if label in ("reactive", "dissociation") and record["tail_persistence_steps"] < 3:
                # A topology that did not hold through the free tail is a reason to
                # distrust the reaction claim, not a reason to throw the quenched endpoint
                # away. That endpoint is a perfectly good stationary point, and the
                # amplitude that produced it is the one nearest the dividing surface --
                # which is exactly where a saddle is reachable. Measured on malonaldehyde,
                # the 0.35 A probe crossed, was discarded whole, and the run found no
                # transition state in sixty trials.
                record.update(status="unconfirmed",
                              failure={"code": "transient_topology",
                                       "message": "New graph did not persist through the last three free MD steps"})
                record["crossing"] = {
                    "observed": True, "tail_persistence_steps": record["tail_persistence_steps"],
                    "meaning": ("the probe crossed but the product did not hold; the "
                                "endpoint is still returned for the minimum check and the "
                                "amplitude still brackets the dividing surface")}
                return Outcome(record, endpoint)
            return Outcome(record, endpoint)
        except GateRejected as exc:
            record.update(status="rejected", classification="unphysical",
                          failure={"code": exc.code, "message": str(exc)})
            if np.isfinite(exc.atoms.positions).all():
                write(str(directory / "rejected.extxyz"), exc.atoms)
                record["rejected_structure"] = f"trials/{trial_id}/rejected.extxyz"
                if exc.code == "ood":
                    record["qm_candidate"] = record["rejected_structure"]
            return Outcome(record)
        except Exception as exc:
            # Failures remain auditable and never enter the reaction network.
            (directory / "error.txt").write_text(traceback.format_exc(), encoding="utf-8")
            record.update(status="failed",
                          failure={"code": "trial_error", "type": type(exc).__name__, "message": str(exc)})
            return Outcome(record)

