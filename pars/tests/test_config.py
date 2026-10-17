from __future__ import annotations

import dataclasses

import pytest

from pars.config import ParsConfig


def test_round_trip_to_dict_from_dict():
    config = ParsConfig()
    restored = ParsConfig.from_dict(config.to_dict())
    assert restored == config


def test_unknown_key_rejected():
    with pytest.raises(ValueError, match="unknown ParsConfig keys"):
        ParsConfig.from_dict({"margin_gamma": 0.2, "typo_option": 1})


def test_out_of_range_values_rejected():
    with pytest.raises(ValueError, match="margin_gamma"):
        ParsConfig(margin_gamma=-0.1)
    with pytest.raises(ValueError, match="refinement_images"):
        ParsConfig(refinement_images=2)
    with pytest.raises(ValueError, match="uncertainty_quantile"):
        ParsConfig(uncertainty_quantile=1.5)
    with pytest.raises(ValueError, match="dominance_mode"):
        ParsConfig(dominance_mode="something_else")


@pytest.mark.parametrize("scale", [-1.0, 0.0, float("nan"), float("inf")])
def test_a_fallback_scale_must_be_a_positive_number(scale):
    """ALG-14: the fallback divides an objective whose pool IQR is zero; a negative one
    made a worse raw value (1) normalise below a better one (0)."""
    with pytest.raises(ValueError, match="objective_scales"):
        ParsConfig.from_dict({"objective_scales": [["energy_excursion", scale]]})


def test_a_positive_fallback_scale_keeps_the_order():
    from pars.frontier.normalize import fit_normalization

    config = ParsConfig.from_dict({"objective_scales": [["energy_excursion", 2.0]]})
    raw = {f"c{i}": {"energy_excursion": 0.0} for i in range(8)}
    raw["worse"] = {"energy_excursion": 1.0}
    normalization = fit_normalization(raw, config)
    assert normalization.scale["energy_excursion"] == 2.0
    assert normalization.center["energy_excursion"] == 0.0


def test_every_field_has_a_documented_rationale_when_nonobvious():
    from pars.config import RATIONALE

    documented = {f.name for f in dataclasses.fields(ParsConfig)}
    undocumented = set(RATIONALE) - documented
    assert undocumented == set(), f"rationale entries for unknown fields: {undocumented}"
