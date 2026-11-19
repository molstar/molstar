# Lane DA: the defaults audit (fast, read-only)

**Goal** (the author, 2026-11-19: "we should have really good defaults for everything so it produces 'pretty' stuff by
default"; "Yes, audit next"): go through every story under Seq* (next)'s Examples and Features in Storybook, in light
and dark at dpr 2 (headless Chrome screenshots on a free port, e.g. 6028; stop the server after), and list what looks
poor *by default*: track heights and empty space (decided rule: layouts.md, "A track's height by default fits its
content, in steps" — propose its numbers: the lane, padding, cap), spacing between tracks and sections, label fit and
truncation (LL's rule is decided: inside if it fits, else beside), font sizes, colours and contrast, ruler density,
legend placement, chips and headers, anything that a reader would call unfinished.

**Output**: `.seq/execution/third-slice/drafts/defaults-audit.md` — **write it in your worktree and report its full
text** (you may not stage `.seq/`; the orchestrator copies it). Per finding: the story, what looks poor (one line),
the default that causes it (file:line or param), a proposed default with its reason, and whether a decided design
already covers it. Group by default (one proposal fixing many stories counts once). Rank by how much it improves the
first impression. Include the screenshot paths (in the session scratchpad
`/private/tmp/claude-501/-Users-davinci001--claude-spares-spare-BGbiLE/3fb68e09-fc5b-4c2f-a070-ae4f27a85f49/scratchpad/da-shots/`).

**Not in scope**: changing code. Read-only.

**Model**: fast (sonnet). **Base**: the commit the orchestrator gives.

**Report**: the draft's full text, the shots' folder, the stories visited and any that failed to render.
