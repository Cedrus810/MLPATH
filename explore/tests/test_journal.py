from explore.journal import Account, Journal, repair_tail


def test_account_holds_caps_and_settles_actual():
    acc = Account("search", 100)
    assert acc.hold(60)
    assert not acc.hold(50)
    acc.settle(60, 7)
    assert (acc.spent, acc.held, acc.available) == (7, 0, 93)
    acc.charge(10)
    assert acc.available == 83


def test_journal_drops_a_torn_last_line_and_appends_cleanly(tmp_path):
    path = tmp_path / "journal.jsonl"
    Journal(path).append({"event": "dispatch", "task": "a"})
    with path.open("a", encoding="utf-8") as fh:
        fh.write('{"event": "sett')
    again = Journal(path)
    assert [(e["event"], e["task"]) for e in again.events] == [("dispatch", "a")]
    assert {"host", "t"} <= set(again.events[0])
    again.append({"event": "settle", "task": "a"})
    assert [e["event"] for e in Journal(path).events] == ["dispatch", "settle"]


def test_repair_tail_leaves_a_complete_file_alone(tmp_path):
    path = tmp_path / "x.jsonl"
    path.write_text('{"a": 1}\n', encoding="utf-8")
    assert repair_tail(path) is False
    assert path.read_text(encoding="utf-8") == '{"a": 1}\n'
    assert repair_tail(tmp_path / "missing.jsonl") is False
