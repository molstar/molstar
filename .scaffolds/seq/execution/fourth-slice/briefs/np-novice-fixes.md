# Lane NP: the novice pass's bugs and stories (the fourth slice)

**Goal**: the findings [novice-findings.md](../novice-findings.md) places on fixes.md (bugs) and NP (stories), by row.
Reproduce each first; a finding that does not reproduce is reported as such, not "fixed". Fix with a test failing on the
parent where code changes; a story's words or composition need no test beyond its twin still passing.

**Bugs**: F5 (Back after a focus does not reach the structure, or the story says the camera is not followed), F7 (the
fake's rows painted empty in Chain maps), F13 (the hide slot and menu cut the label "Domains" in the default column),
F15 (after a header click, Cmd/Ctrl+Z does nothing and Back is disabled), F23 (the size legend's key: "Variants" twice
in the compact band, no circles in the side panel, no unit).

**Stories**: F8 (Chain maps hovers a zoomed window so a one-residue shift shows), F12 (Host data's B-factor legend; check
the end gap is the data), F17 (a ring's menu titled by its track), F18 (Host widget's words), F20 (the lollipop keeps
its label; the legend lists what is in view), F21 (the text says what Alt+double-click does), F27 (the host-drawn
navigator's labels do not overprint), F28 (its readout names the hovered range), F29 (the archive stories: a save and
an open button, a line of text above Refused).

**Out**: the design questions (F1, F3, F9, F11, F14, F22, F25) and the TODO rows: the orchestrator's.

**Model**: fast (sonnet). **Where**: a worktree on `seq/np-4` from the base named. Hot files are the orchestrator's.

**Report**: under 300 words: per finding the outcome and commit (or why not), tests failing on the parent, the stories
to look at, the `.v6/changelogs/` entry's text.
