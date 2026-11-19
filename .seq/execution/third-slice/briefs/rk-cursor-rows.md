# Lane RK: the keyboard cursor's rows

**Goal** (the author, 2026-11-19: "Yes, after HC"; found by NP): build the keyboard cursor's rows as
[grids-and-rows](../../../design/grids-and-rows.md), "The keyboard cursor" (EC-154) designs them: the cursor
`{ section, system, position, ySystem, yPosition, part }`, `part` `cell` or `header`; ← from the first column onto the
header and → back; ↑ ↓ step display rows (pinned band to scrolled rows, a summary row one stop); Page Up and Page Down a
window of rows (the cursor moves a page and reveals itself; NP's page commands scroll only when no cursor is set);
Cmd/Ctrl+↑ ↓ the first and last row; Shift+arrows extend the block; Enter focuses the row or cell, or opens a summary
row; Space (on key-up, not after a drag while held) toggles the cell in the selection, and on the header the row's
clade. Moves are arithmetic on display rows, never on loaded rows (rows never loaded are stepped over by index). The
cursor's drawing (a squircle per ux-foundations §5), the announcer's words (its row's label), and HC's hover card at
the cursor follow it. Grids: ↑ ↓ step the grid's y the same way.

**Tests** (fail on the parent): each move above; Page Down with pinned rows and summary rows; the first and last row;
Space on a cell and on a header; a row never loaded stepped over; the announcement and the card at the cursor; the
grid's y.

**Not in scope**: anything else. CS (p53's chromosome) runs beside you in `packages/seq/examples` and
`packages/seq/storybook`: stay out of its view files; a shared story file is a merge the orchestrator resolves.

**Model**: strong (opus). **Base**: `packages-seq` at the commit this brief is in: `git checkout -B seq/rk-cursor-rows
<that commit>`. Storybook port 6019 (your own; stop it after).

**Checks**: as every lane's (see [hc-hover-card.md](hc-hover-card.md)): unit, medium, tsc, typecheck:test, storybook
typecheck, test:perf, lint, format, size; folder durations with `uptime`. Look at Globins › Rows (Family, Lazy) in light
and dark.

**Headers**, **hot files** and **stop when** (a public API past the cursor's `ySystem`/`yPosition`/`part` and its
commands): as every lane's (see [qt-quick-lane.md](qt-quick-lane.md)).

**Report**: commits; checks; durations; sizes; what you saw; drift; followups.
