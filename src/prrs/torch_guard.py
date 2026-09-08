"""The one wrapper around ``torch.load``, installed once and never taken off.

Two separate needs used to wrap ``torch.load`` independently, and the combination was the
defect -- not either one alone.

1. ``weights_only``. torch 2.6 flipped its default to True. e3nn's ``constants.pt``
   contains a pickled ``slice``, and a MACE checkpoint is a full pickled model, so both
   fail under the new default. mace's answer is to set
   ``TORCH_FORCE_NO_WEIGHTS_ONLY_LOAD=1`` itself, on import, and
   ``torch/serialization.py`` then warns on every load where the callsite left the
   argument unset. Unsetting the variable does not work: mace puts it back the next time
   it is imported. Passing ``weights_only`` explicitly at the callsite does, because that
   is the exact condition torch checks -- and it is also narrower than the variable, since
   only the checkpoints a caller has declared get the full-pickle exemption and everything
   else keeps torch's safe default.

2. The device move. openmm-ml's local-``modelPath`` MACE branch loads with ``map_location``
   and never calls ``.to(device)``, so constants inside the scripted submodules stay on the
   host and every CUDA evaluation fails.

Wrapping twice is what made this dangerous. The device-move wrapper saved whatever
``torch.load`` was on entry and restored it unconditionally on exit, so installing the
``weights_only`` wrapper while that context was open made the context's ``finally`` delete
it -- silently, with no error, the warning simply returning and nothing pointing at why.
Order decided whether the guard survived, and order is not something a caller should have
to reason about.

So there is one wrapper. It is installed at most once per process, it is never restored,
and the device move is a flag on it rather than a second wrapper. Nothing here has an
unwind path, which is the property that makes the failure mode impossible instead of
merely unlikely.
"""

from contextlib import contextmanager
from pathlib import Path

# The wrapper this module installed, or None. Compared by identity, so a caller can ask
# whether something else has since replaced it.
_INSTALLED = None
# Resolved paths a caller has declared loadable as a full pickle. Everything else gets
# torch's safe default, so a dependency that genuinely needs weights_only=False raises an
# UnpicklingError naming itself instead of being quietly permitted.
_TRUSTED = set()
# Depth, not a boolean: nesting must not let an inner exit switch the move off.
_MOVE_DEPTH = 0


def _resolved(target):
    """The absolute path of a torch.load target, or None if it is not a path at all."""
    if isinstance(target, (str, Path)):
        try:
            return str(Path(target).resolve())
        except OSError:
            return None
    return None


def trust(path):
    """Declare one checkpoint loadable as a full pickle.

    Scoped to the file, not to the process: a MACE checkpoint genuinely cannot be
    reconstructed by weights-only loading, but that is a statement about that file and
    should not become a blanket exemption for every pickle the process later opens.
    """
    resolved = _resolved(path)
    if resolved is None:
        raise ValueError(f"cannot trust a non-path torch.load target: {path!r}")
    _TRUSTED.add(resolved)
    return resolved


def weights_only_for(target):
    """Whether this target should load under torch's safe default.

    Split out from the wrapper so the decision is testable without torch present.
    """
    resolved = _resolved(target)
    return not (resolved is not None and resolved in _TRUSTED)


def install(trusted=()):
    """Wrap torch.load once. Idempotent, and there is deliberately no uninstall.

    `trusted` is declared at install time rather than by a later `trust` call so that a
    caller which must have its list in force before anything is loaded can hand it over
    in one step. Installing and then trusting leaves a window in between, which matters
    for a guard whose whole job is to be in place before the first import.

    Registering ``slice`` first means e3nn's constants.pt is a real weights-only load
    rather than an exemption. This has to run before e3nn or mace are imported, because
    that import is when the file is read; a caller that imports them first gets the
    UnpicklingError and no guard can retroactively help.
    """
    import torch

    torch.serialization.add_safe_globals([slice])
    global _INSTALLED, _TRUSTED
    if _INSTALLED is not None:
        for path in trusted:
            trust(path)
        return _INSTALLED
    original = torch.load

    # Is another guard of this kind already underneath? Two of them can legitimately
    # coexist -- the p3 bundle installs its own before importing anything -- and then the
    # outer one must not answer a question the inner one is better placed to answer.
    inner_guard = bool(getattr(original, "_prrs_wrapped", False))
    shared = getattr(original, "_prrs_trusted", None)
    if shared is not None:
        # One registry between them, so `trust` from either side is seen by both.
        _TRUSTED = shared

    def load(*args, **kwargs):
        """Set ``weights_only`` explicitly, which is the condition torch warns on.

        With an inner guard present the rule narrows to: override only where this guard
        has a positive reason to, never to impose the restrictive default. Imposing it
        was a real failure -- an outer guard that had not been told to trust the
        checkpoint set weights_only=True, the inner guard saw the argument already
        present and stood down as designed, and the model failed to unpickle even though
        the inner guard's own list would have allowed it. The restrictive default is only
        safe to apply when nothing below can know better.
        """
        target = kwargs.get("f", args[0] if args else None)
        if "weights_only" not in kwargs:
            trusted = not weights_only_for(target)
            if trusted:
                kwargs["weights_only"] = False
            elif not inner_guard:
                kwargs["weights_only"] = True
        loaded = original(*args, **kwargs)
        if _MOVE_DEPTH:
            device = kwargs.get("map_location")
            if device is not None and hasattr(loaded, "to"):
                return loaded.to(device)
        return loaded

    load._prrs_wrapped = True
    load._prrs_trusted = _TRUSTED
    load._prrs_defers_to_inner = inner_guard
    torch.load = load
    _INSTALLED = load
    for path in trusted:
        trust(path)
    return load


def outermost():
    """Is this module's wrapper the one ``torch.load`` currently names?

    It reports an ordering fact, and the name says so, because two different situations
    both make it False and only one of them is a problem:

        REMOVED  something restored torch.load over the top of this guard, so the guard
                 is gone from the chain entirely -- the failure this module exists to
                 make impossible.
        NESTED   something wrapped this guard, which still runs underneath and still
                 does its job.

    Telling those apart cheaply is not possible: it would mean walking the closure of an
    arbitrary wrapper. So False means "go and look", not "broken".
    """
    import torch

    return _INSTALLED is not None and torch.load is _INSTALLED


# The earlier name, which claimed more than the check can deliver. Kept so callers and
# recorded manifests do not break.
intact = outermost


@contextmanager
def moved_to_device():
    """Make loads inside this block honour their ``map_location``.

    A flag on the existing wrapper, never a second wrapper, so leaving the block cannot
    remove anything. It stays scoped rather than becoming permanent because
    test_backend_parity pins the unpatched upstream behaviour: openmm-ml failing on CUDA
    without the move is what tells us the workaround is still needed, and a global move
    would erase that signal.
    """
    global _MOVE_DEPTH
    install()
    _MOVE_DEPTH += 1
    try:
        yield
    finally:
        _MOVE_DEPTH -= 1


def state():
    """What the wrapper is doing right now, for a run manifest."""
    return {
        "wrapper_installed": _INSTALLED is not None,
        "wrapper_outermost": outermost() if _INSTALLED is not None else False,
        "defers_to_inner_guard": bool(getattr(_INSTALLED, "_prrs_defers_to_inner", False)),
        "trusted_checkpoints": sorted(_TRUSTED),
        "device_move_active": bool(_MOVE_DEPTH),
    }
