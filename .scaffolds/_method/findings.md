# Findings: <work id>

<!-- The template for a work folder's findings, `execution/<work-id>/findings.md` (§used-before-done): append-only, one
     row per merged lane and one per gate, the novice pass and the author's use, so the fold-back can count which gate
     catches what. Written by the orchestrator at each merge (loop, Each step: record) and never rewritten; a correction
     is a new row. Findings are counted by kind: a bug, a missing design, a look. The Items cell answers every numbered item of the lane's brief,
     `<n> ✓ <proof: the test or output by name>` or `<n> ✗ not landed: <where it went: fixes.md, followups>`, separated
     by `;`; `_checks/features-planned.ts --brief <brief.md> --findings <findings.md>` fails an item answered neither
     way, and a lane with no row (the first cell is the brief's lane id). Delete this comment once filled. -->

| Lane or gate | Merged | Model | Items | Bug | Missing design | Look | Who looked | Notes |
|---|---|---|---|---|---|---|---|---|
| <LANE> | <date> (`<merge commit>`) | <strong or fast> | `1 ✓ <proof>; 2 ✗ not landed: <where it went>`, every item of the brief | <n (what)> | <n> | <n> | <who, at what: light and dark, the gestures> | <the review's refuted or unfixed findings> |
| GATE <n> | <date> | | | | | | | <the gate's findings> |
