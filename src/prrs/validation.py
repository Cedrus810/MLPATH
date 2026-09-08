"""Independent cross-checks of the curvature layer. Nothing here is used by the search.

The curvature numbers the search reports are, at present, self-consistent: one
implementation projects the rigid-body modes out, mass-weights, and diagonalizes, and
nothing else has ever computed the same quantity. A second implementation is the only
thing that can find a mistake in that convention, which is the same reason the OpenMM
backend exists as a numerical cross-check rather than as a faster path.

Two checks, and they are not interchangeable:

  the same Hessian through someone else's spectral machinery
      catches a mass-weighting, ordering or unit-conversion mistake, and costs nothing

  someone else's Hessian as well
      catches the above and any mistake in the differentiation, at 6N force evaluations

This module lives outside the search on purpose: it belongs to the analysis layer, like
`openmm_backend`, and importing it must never be necessary to run a search.
"""

from pathlib import Path
from tempfile import TemporaryDirectory
import numpy as np
from ase import units
from .reliability import GuardedCalculator
from .runner import _mass_weighted_hessian, _trivial_modes

# hbar * omega in eV for a mass-weighted eigenvalue in eV/(Angstrom^2 amu).
ENERGY_PER_ROOT_EIGENVALUE = units._hbar * 1e10 / np.sqrt(units._e * units._amu)


def eigenvalues_from_energies(energies):
    """Undo the vibrational-energy convention: lambda = (hbar omega / scale)^2, signed.

    A negative eigenvalue is reported by every code as an imaginary energy, so the sign
    has to be read off the imaginary part rather than inferred.
    """
    values = []
    for energy in np.asarray(energies):
        magnitude = (abs(energy) / ENERGY_PER_ROOT_EIGENVALUE) ** 2
        imaginary = abs(np.imag(energy)) > abs(np.real(energy))
        values.append(-magnitude if imaginary else magnitude)
    return np.sort(np.array(values))


def internal_eigenvalues(eigenvalues, trivial_count):
    """Drop the rigid-body modes by magnitude, whichever convention produced them.

    Our own spectrum has them projected to zero; a code that does not project leaves them
    as small numbers of either sign. Both are handled by discarding the trivial_count
    smallest in absolute value, which is only safe because a genuine mode of a molecule is
    orders of magnitude larger -- 0.2 against 1e-15 on ethanol.
    """
    order = np.argsort(np.abs(eigenvalues))
    return np.sort(np.asarray(eigenvalues)[order[trivial_count:]])


def cartesian_hessian(candidate, factory, config, source="auto"):
    """Our Hessian in eV/Angstrom^2, undoing the mass weighting exactly."""
    mass_weighted, basis, provenance = _mass_weighted_hessian(
        candidate, factory, config, source
    )
    root = np.sqrt(np.repeat(candidate.get_masses(), 3))
    return mass_weighted * np.outer(root, root), basis, provenance


def ase_energies_from_our_hessian(candidate, factory, config, source="auto"):
    """Our Hessian, ASE's mass weighting and unit conversion. No force evaluations.

    ASE does not project the rigid-body modes out, so its six smallest values are the
    error floor rather than zeros; `internal_eigenvalues` removes them from both sides.
    """
    from ase.vibrations import VibrationsData

    hessian, basis, provenance = cartesian_hessian(candidate, factory, config, source)
    count = len(candidate)
    data = VibrationsData(candidate, hessian.reshape(count, 3, count, 3))
    return data.get_energies(), len(basis), provenance


def ase_vibrations_energies(candidate, factory, config, delta=None):
    """ASE's own central-difference Hessian and spectrum. 6N force evaluations.

    Its step is its own, so any disagreement at the size of a truncation error is exactly
    that and not a defect: two finite differences at different steps are two different
    approximations of the same limit.
    """
    from ase.vibrations import Vibrations

    atoms = candidate.copy()
    guard = GuardedCalculator(factory(), config, None)
    atoms.calc = guard
    step = config.minimum_check_step_A if delta is None else float(delta)
    with TemporaryDirectory() as directory:
        vibrations = Vibrations(atoms, name=str(Path(directory) / "vib"), delta=step, nfree=2)
        vibrations.run()
        energies = vibrations.get_energies()
    trivial = len(_trivial_modes(candidate.positions, candidate.get_masses()))
    return (
        energies,
        trivial,
        {
            "source": "ase.vibrations",
            "delta_A": step,
            "force_evaluations": 6 * len(candidate),
            "meaning": "independent central differences, unprojected",
        },
    )
