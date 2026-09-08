"""The single torch.load wrapper, tested against a fake torch so it always runs.

These tests take a fake `torch` rather than the real one on purpose. The behaviour at
issue is the wrapping mechanism -- who is installed, what survives a context exit, what
argument reaches the underlying load -- and none of that needs a tensor library. Testing
it against the real torch would put the whole file behind an importorskip, and a
module-level importorskip collects zero tests while pytest still reports green. That
failure mode is what let the backend-parity suite go unrun without anyone noticing, and
it is the same class of defect as the one this module exists to remove.
"""

import sys
import types

import pytest


@pytest.fixture
def fake_torch(monkeypatch):
    """A stand-in exposing exactly the surface torch_guard touches."""
    calls = []

    def load(*args, **kwargs):
        calls.append({"args": args, "kwargs": dict(kwargs)})
        return kwargs.get("_result", _Loadable())

    torch = types.ModuleType("torch")
    torch.load = load
    torch.serialization = types.SimpleNamespace(
        safe_globals=[], add_safe_globals=lambda g: torch.serialization.safe_globals.extend(g)
    )
    torch.calls = calls
    torch.original_load = load
    monkeypatch.setitem(sys.modules, "torch", torch)

    from prrs import torch_guard

    monkeypatch.setattr(torch_guard, "_INSTALLED", None)
    monkeypatch.setattr(torch_guard, "_TRUSTED", set())
    monkeypatch.setattr(torch_guard, "_MOVE_DEPTH", 0)
    return torch


class _Loadable:
    """Something with .to(), like a model; records where it was asked to go."""

    def __init__(self):
        self.moved_to = None

    def to(self, device):
        self.moved_to = device
        return self


def test_install_is_idempotent_and_registers_slice(fake_torch):
    from prrs import torch_guard

    first = torch_guard.install()
    second = torch_guard.install()
    assert first is second
    assert fake_torch.load is first
    assert slice in fake_torch.serialization.safe_globals
    assert torch_guard.intact()


def test_weights_only_is_always_passed_explicitly(fake_torch, tmp_path):
    """torch warns only when the callsite leaves the argument unset, so setting it is the
    fix at the cause; the env var mace re-asserts on import is then irrelevant."""
    from prrs import torch_guard

    torch_guard.install()
    checkpoint = tmp_path / "model.pt"
    checkpoint.write_bytes(b"x")

    fake_torch.load(str(checkpoint))
    assert fake_torch.calls[-1]["kwargs"]["weights_only"] is True

    torch_guard.trust(checkpoint)
    fake_torch.load(str(checkpoint))
    assert fake_torch.calls[-1]["kwargs"]["weights_only"] is False

    # An untrusted neighbour keeps the safe default even in the same directory.
    other = tmp_path / "other.pt"
    other.write_bytes(b"x")
    fake_torch.load(str(other))
    assert fake_torch.calls[-1]["kwargs"]["weights_only"] is True


def test_an_explicit_caller_argument_is_never_overridden(fake_torch, tmp_path):
    from prrs import torch_guard

    torch_guard.install()
    fake_torch.load(str(tmp_path / "m.pt"), weights_only=False)
    assert fake_torch.calls[-1]["kwargs"]["weights_only"] is False


def test_device_move_is_scoped_and_leaves_nothing_behind(fake_torch):
    """The move must stay inside createSystem: test_backend_parity pins that openmm-ml
    fails on CUDA *without* it, and a permanent move would erase that signal."""
    from prrs import torch_guard

    torch_guard.install()
    installed = fake_torch.load

    outside = fake_torch.load("m.pt", map_location="cuda")
    assert outside.moved_to is None

    with torch_guard.moved_to_device():
        inside = fake_torch.load("m.pt", map_location="cuda")
    assert inside.moved_to == "cuda"

    after = fake_torch.load("m.pt", map_location="cuda")
    assert after.moved_to is None
    # The wrapper is the same object it was before: nothing was installed or restored.
    assert fake_torch.load is installed
    assert torch_guard.intact()


def test_nested_scopes_do_not_switch_the_move_off_early(fake_torch):
    from prrs import torch_guard

    torch_guard.install()
    with torch_guard.moved_to_device():
        with torch_guard.moved_to_device():
            pass
        still_inside = fake_torch.load("m.pt", map_location="cuda")
    assert still_inside.moved_to == "cuda"


def test_the_scope_does_not_delete_a_wrapper_installed_while_it_was_open(fake_torch):
    """The defect this module removes.

    The old context manager saved torch.load on entry and restored it unconditionally on
    exit, so anything installed while it was open vanished when it closed -- with no
    error. Here the equivalent sequence must leave the guard in place, because the scope
    is a flag rather than a wrapper and has nothing to restore.
    """
    from prrs import torch_guard

    torch_guard.install()
    installed = fake_torch.load

    with torch_guard.moved_to_device():

        def someone_elses_wrapper(*args, **kwargs):
            return installed(*args, **kwargs)

        fake_torch.load = someone_elses_wrapper

    assert fake_torch.load is someone_elses_wrapper  # not clobbered back
    assert not torch_guard.intact()  # and the guard says so


def test_intact_reports_a_replaced_wrapper(fake_torch):
    from prrs import torch_guard

    torch_guard.install()
    assert torch_guard.intact()
    fake_torch.load = lambda *a, **k: None
    assert not torch_guard.intact()


def test_state_is_recordable_in_a_manifest(fake_torch, tmp_path):
    import json
    from prrs import torch_guard

    checkpoint = tmp_path / "m.pt"
    checkpoint.write_bytes(b"x")
    torch_guard.install()
    torch_guard.trust(checkpoint)
    state = torch_guard.state()
    assert state["wrapper_installed"] is True
    assert state["device_move_active"] is False
    assert str(checkpoint.resolve()) in state["trusted_checkpoints"]
    json.dumps(state)  # must not raise


def test_trust_refuses_a_non_path_target(fake_torch):
    from prrs import torch_guard

    with pytest.raises(ValueError):
        torch_guard.trust(object())


def _bundle_style_guard(torch, trusted_paths):
    """A second guard of the same family, the shape scripts/p3_node_bundle uses.

    Its trust list is its own closure, not shared, which is the situation that made
    stacking dangerous: the outer guard cannot read it and must not overrule it.
    """
    from pathlib import Path

    original = torch.load
    trusted = {str(Path(p).resolve()) for p in trusted_paths}

    def load(*args, **kwargs):
        target = kwargs.get("f", args[0] if args else None)
        if "weights_only" not in kwargs:
            resolved = str(Path(target).resolve()) if isinstance(target, (str, Path)) else None
            kwargs["weights_only"] = resolved not in trusted if resolved else True
        return original(*args, **kwargs)

    load._prrs_wrapped = True
    torch.load = load
    return load


def test_an_outer_guard_never_overrules_an_inner_one_into_failing(fake_torch, tmp_path):
    """The defect a second session measured on the real stack.

    An outer guard that had not been told to trust the checkpoint set weights_only=True;
    the inner guard saw the argument already present and stood down as designed; the
    model then failed to unpickle even though the inner guard's own list allowed it. The
    outer guard may only override where it has a positive reason to.
    """
    from prrs import torch_guard

    checkpoint = tmp_path / "model.pt"
    checkpoint.write_bytes(b"x")

    _bundle_style_guard(fake_torch, [checkpoint])  # inner knows the checkpoint
    torch_guard.install()  # outer does NOT

    fake_torch.calls.clear()
    fake_torch.load(str(checkpoint))
    # The inner guard's answer survives: a full pickle, as it intended.
    assert fake_torch.calls[-1]["kwargs"]["weights_only"] is False


def test_an_outer_guard_still_grants_its_own_trust_through_an_inner_one(fake_torch, tmp_path):
    from prrs import torch_guard

    checkpoint = tmp_path / "outer_only.pt"
    checkpoint.write_bytes(b"x")

    _bundle_style_guard(fake_torch, [])  # inner knows nothing
    torch_guard.install()
    torch_guard.trust(checkpoint)  # outer does

    fake_torch.calls.clear()
    fake_torch.load(str(checkpoint))
    assert fake_torch.calls[-1]["kwargs"]["weights_only"] is False


def test_the_restrictive_default_still_applies_when_nothing_is_underneath(fake_torch, tmp_path):
    """Deferring is only correct when something below can know better.

    Alone, the guard must still set the argument on every call -- that is what silences
    torch, and leaving it unset would bring the warning straight back.
    """
    from prrs import torch_guard

    torch_guard.install()
    fake_torch.calls.clear()
    fake_torch.load(str(tmp_path / "unknown.pt"))
    assert fake_torch.calls[-1]["kwargs"]["weights_only"] is True


def test_a_shared_trust_registry_is_adopted_when_one_is_offered(fake_torch, tmp_path):
    """Two guards that both expose a registry keep one between them, so a trust granted
    on either side is seen by both rather than each holding half the answer."""
    from prrs import torch_guard

    checkpoint = tmp_path / "m.pt"
    checkpoint.write_bytes(b"x")
    shared = set()

    def inner(*args, **kwargs):
        return fake_torch.original_load(*args, **kwargs)

    inner._prrs_wrapped = True
    inner._prrs_trusted = shared
    fake_torch.load = inner

    torch_guard.install()
    torch_guard.trust(checkpoint)
    assert str(checkpoint.resolve()) in shared


def test_outermost_reports_ordering_not_health(fake_torch):
    """False has two causes and only one is a problem, so the name says ordering."""
    from prrs import torch_guard

    torch_guard.install()
    ours = fake_torch.load
    assert torch_guard.outermost()
    assert torch_guard.intact is torch_guard.outermost  # old name still works

    def someone_wrapped_us(*args, **kwargs):
        return ours(*args, **kwargs)

    fake_torch.load = someone_wrapped_us
    # NESTED: still in the chain and still working, but no longer outermost.
    assert not torch_guard.outermost()
    assert torch_guard.state()["wrapper_outermost"] is False


def test_trust_can_be_declared_at_install_time(fake_torch, tmp_path):
    """Installing and then trusting leaves a window; a guard that must be in force before
    the first import needs to hand over its list in one step."""
    from prrs import torch_guard

    first, second = tmp_path / "a.pt", tmp_path / "b.pt"
    for path in (first, second):
        path.write_bytes(b"x")
    torch_guard.install(trusted=[first, second])
    fake_torch.calls.clear()
    fake_torch.load(str(first))
    assert fake_torch.calls[-1]["kwargs"]["weights_only"] is False
    # A second install call is still idempotent but may add to the list.
    third = tmp_path / "c.pt"
    third.write_bytes(b"x")
    torch_guard.install(trusted=[third])
    fake_torch.load(str(third))
    assert fake_torch.calls[-1]["kwargs"]["weights_only"] is False
