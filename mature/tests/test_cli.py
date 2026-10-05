"""The CLI shell refuses fail-closed until its stage lands (M1..M4)."""

import pytest

from mature.cli import main


def test_each_stage_refuses_with_exit_code_2(capsys):
    for stage in ("prepare", "shoot", "report"):
        assert main([stage]) == 2
        err = capsys.readouterr().err
        assert "not implemented" in err
        assert stage in err


def test_help_lists_all_three_stages(capsys):
    with pytest.raises(SystemExit) as excinfo:
        main(["--help"])
    assert excinfo.value.code == 0
    out = capsys.readouterr().out
    for stage in ("prepare", "shoot", "report"):
        assert stage in out


def test_unknown_command_is_an_argparse_error():
    with pytest.raises(SystemExit) as excinfo:
        main(["refine"])
    assert excinfo.value.code == 2
