"""The CLI shape: three verbs, each with its required arguments (M4).

End-to-end runs of the verbs live in test_maturation.py, on a real small
PRRS run.
"""

import pytest

from mature.cli import main


def test_help_lists_all_three_stages(capsys):
    with pytest.raises(SystemExit) as excinfo:
        main(["--help"])
    assert excinfo.value.code == 0
    out = capsys.readouterr().out
    for stage in ("prepare", "shoot", "report"):
        assert stage in out


@pytest.mark.parametrize("stage", ["prepare", "shoot", "report"])
def test_a_stage_without_its_batch_directory_is_an_argparse_error(stage):
    with pytest.raises(SystemExit) as excinfo:
        main([stage])
    assert excinfo.value.code == 2


def test_unknown_command_is_an_argparse_error():
    with pytest.raises(SystemExit) as excinfo:
        main(["refine"])
    assert excinfo.value.code == 2
