"""Terminal-cell formatting for the FDM workflow's own output."""

from __future__ import annotations

import unicodedata


def display_width(text: str) -> int:
    """Return the terminal-cell width of an ANSI-free string."""
    width = 0
    for char in text:
        if unicodedata.combining(char):
            continue
        if unicodedata.category(char).startswith("C"):
            continue
        width += 2 if unicodedata.east_asian_width(char) in ("W", "F") else 1
    return width


def pad_right(text: str, width: int) -> str:
    """Pad ``text`` on the right to exactly ``width`` terminal cells."""
    return text + " " * max(0, width - display_width(text))


def render_table(headers: list[str], rows: list[list[str]],
                 right_aligned: set[int] | None = None,
                 widths: list[int] | None = None) -> list[str]:
    """Render an ASCII-separated table using terminal-cell column widths."""
    all_rows = [headers, *rows]
    column_widths = widths or [max(display_width(row[index]) for row in all_rows)
                               for index in range(len(headers))]
    right = right_aligned or set()
    rendered = []
    for row in all_rows:
        cells = []
        for index, cell in enumerate(row):
            padding = max(0, column_widths[index] - display_width(cell))
            cells.append(" " * padding + cell if index in right
                         else cell + " " * padding)
        rendered.append(" | ".join(cells))
    return rendered
