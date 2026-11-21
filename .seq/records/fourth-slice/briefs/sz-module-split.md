# Lane SZ: the files over 1500 lines split by what changes together

**Task**: §module-size (oz, architecture rules) bounds a source file at 1500 lines; `.seq/_checks/file-size.ts` lists
the files already past it in `.seq/_checks/file-size-baseline.json`, which only shrinks. Split them:
`core/src/render/tracks.ts` (2882 at the rule), `core/src/render/rows.ts` (2520), `core/src/render/ranges.ts` (1543),
`viewer/src/viewer.ts` (1570); look at `viewer/src/models/view.ts` (1466) and say whether it is one concern. Read each
file first and name its concerns (for a layer: its layout, its draw, its hit-test, its label placement; for the
viewer: composition, the frame, rails and navigators, the UI mount) before moving code; split along them, never at a
line count. Public exports and their paths stay (a barrel re-exports what moved); behavior, draw lists and allocation
unchanged: the quick, medium and perf tiers pass with no test edited but for imports. Each new module gets the Mol\*
header (the moved code's authors kept) and a first comment saying its concern.

**Output**: the splits, one commit a file; every file under 1500 (aim under 800 where a concern allows); the `_spec`s'
module lists updated. The baseline is the orchestrator's (`file-size.ts --update` at the merge).

**Model**: strong. **Where**: a worktree, branch `seq/sz-module-split`, after FX1 and before the header lanes, while no other lane is in the files it splits (plan, Order and parallelism).

**Report**: under 250 words: each file's concerns and the modules it became (lines each), anything left whole and why,
the checks.
