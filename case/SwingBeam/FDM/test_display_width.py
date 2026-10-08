#!/usr/bin/env python3

import contextlib
import importlib.util
import io
from pathlib import Path

from display_width import display_width, pad_right, render_table


ROOT = Path(__file__).resolve().parent


def load_script(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / (name + ".py"))
    if spec is None or spec.loader is None:
        raise RuntimeError("cannot load %s" % name)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_display_width_handles_cjk_and_combining_marks():
    assert display_width("abc") == 3
    assert display_width("网格") == 4
    assert display_width("网格 cells") == 10
    assert display_width("：") == 2
    assert display_width("e\u0301") == 1


def test_pad_right_uses_terminal_cells():
    padded = pad_right("网格", 6)
    assert display_width(padded) == 6
    assert padded == "网格  "


def test_render_table_keeps_cjk_columns_aligned():
    lines = render_table(["阶数", "论文右列VE (Hz)"], [["1", "35.5200"]], {0, 1})
    assert [display_width(line) for line in lines] == [display_width(lines[0])] * 2
    assert lines[0].index("|") == lines[1].index("|")


def test_modal_comparison_prints_display_aligned_table():
    module = load_script("02_run_modal_virtual_elastic")
    stream = io.StringIO()
    with contextlib.redirect_stdout(stream):
        assert module.print_comparison([6.11, 35.52], [6.11, 35.52], [6.02, 34.98]) == 2
    lines = [line for line in stream.getvalue().splitlines() if "|" in line]
    assert len({display_width(line) for line in lines}) == 1


def test_modal_status_colons_share_the_same_display_column():
    module = load_script("02_run_modal_virtual_elastic")
    stream = io.StringIO()
    with contextlib.redirect_stdout(stream):
        module.print_status("流体块来源", "mesh")
        module.print_status("FSI 界面面片", "250")
        module.print_status("流体 min_vol", "8.7e-09")
    lines = stream.getvalue().splitlines()
    colon_columns = [display_width(line.split(" : ", 1)[0]) for line in lines]
    assert len(set(colon_columns)) == 1


def test_progress_bar_has_fixed_cell_width_and_final_newline():
    module = load_script("04_run_fdm")
    module._last_progress_width = 0
    stream = io.StringIO()
    with contextlib.redirect_stdout(stream):
        module._render_progress("求解", 0.25, 0.5, 2.0, 100.0, final=True)
    output = stream.getvalue()
    assert output.endswith("\n")
    assert display_width(output.rstrip("\n")) > module.BAR_WIDTH
