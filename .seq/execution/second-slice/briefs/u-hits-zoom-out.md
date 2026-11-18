# Lane U: missing hit boxes, and zooming out

**Goal**, from the author's use of the stories (2026-11-18):
1. **Some marks have no hit box** (the author: "some boxes don't have a hit box"). In Examples/p53/Gene by region at
   the fit view, the hover outline appears on some domain boxes on the "Domains T1-domains" track and not on others
   (the upper-right box of that track, a domain carried onto the region through its transcript's codon map, shows
   none). Find which marks the hit index misses (a lane past the first, a piece of an item split by an intron, an item
   arriving in a later tile, a nested source's dynamic node), fix it with a test shown failing on the parent: every
   drawn mark of a ranges layer is hit at its centre, in every lane and piece, on that view and on Protein on its gene.
   A property test over random items and lanes is welcome.
2. **Zooming out is hard** (the author: "double click on something zooms it, but it's really hard to zoom out").
   Today: a double-click on an element focuses it; on empty space it fits; Cmd/Ctrl+wheel and pinch zoom. Add, each a
   binding or a key in the presets, configurable (the author: "the gestures are configurable anyway"):
   - **a view history of jumps** (session state, never the undo stack: ux-foundations section 7 "view state is not
     history"): a focus, a fit, a reveal, a bookmark's or a link's move pushes the viewport it left; **Back** returns
     (Alt+Left, the mouse's back button, Backspace when the viewer has focus) and **Forward** goes again (Alt+Right);
     a pan or zoom by hand does not push but a jump after it does; bounded (, 50);
   - **a double-click on what is already focused goes back** (the view equals that element's focus);
   - **Shift+double-click zooms out** about the pointer by the double-click's own factor; **-** and **=** (and **+**)
     zoom out and in about the cursor or the centre; **0** fits;
   - commands for each (`view.back`, `view.forward`, `view.zoomOut`, `view.zoomIn`, `view.fit` exists), so a host's
     toolbar reaches them; jumps animate as focus does (the motion of ux-foundations section 7), reduced motion cuts.
   Followers keep their rules (a follower leads: lane T); a back restores every leading section's range together.

Not in scope: collapsed introns, the reversed follower (the third slice), anything in `packages/core`.

**Read first**: `.seq/design/ux-foundations.md` (sections 6 to 8), `viewer-models.md`, `layout-and-draw-lists.md`
(the hit index); the viewer and core `_spec`s; `packages/seq/viewer/src/{models/view.ts,input.ts,bindings.ts,commands.ts,behaviors}`,
`packages/seq/core/src/render/{ranges.ts,marks.ts,hit*.ts}`; `packages/seq/AGENTS.md`.

**Model**: strong (opus).

**Base**: `packages-seq` at `fbdd8a105`. Your worktree may be cut from master: with a clean tree, first
`git checkout -B seq/u-hits fbdd8a105`, and check `git merge-base --is-ancestor fbdd8a105 HEAD`.

**Touches**: `packages/seq/viewer/src/**`, `packages/seq/core/src/{render,view}/**` for the hit index,
`packages/seq/storybook/**` (the Features stories' hints name Back and the zoom keys; the gesture tables), the
`_spec`s.

**Tests**: headless twins on a fake clock: the hit at every mark's centre (both views); back and forward over a focus,
a fit and a bookmark jump; a double-click on the focused element going back; Shift+double-click and the keys.

**Performance representative**: hit-testing stays allocation-free warm (the opt-in tier); a back step under 4 ms
(`itWithin`, after warm-ups).

**Configurable over fixed**: the bindings and keys, `view.historyLimit`, whether a double-click on the focused element
goes back, in the `_spec`s' Configuration tables.

**Checks**: each its own command, green only on exit 0 with a non-zero test count.
- `pnpm test:unit packages/seq/`
- `pnpm test:heavy packages/seq/`
- `npx tsc -b packages/seq/core packages/seq/viewer packages/seq/examples`
- `pnpm -C packages/seq/core typecheck:test`, the same for `viewer` and `examples`, and `pnpm -C packages/seq/storybook typecheck`
- `pnpm lint`, `pnpm format:check`

**Headers**: a new file has the Mol\* header `Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for
more info.` with `@author David Sehnal <david.sehnal@gmail.com>`; a modified file keeps its header, its year range
extended to 2026 and that author line added if missing. No mention of `packages/sequences` or `.seq/` in code, comments,
tests or `_spec`s; code cites its own `_spec` sections.

**Fixes land with a test shown failing on the fix's parent**: the report names the command.

**Hot files, never staged**: `pnpm-lock.yaml`, `scripts/workspace/inventory.json`, the root `tsconfig.json`,
`biome.json`, `vitest.config.ts`, `storybook/packages.ts`, `.v6/changelogs/` (write the entry in your report), every
`package.json` and `tsconfig*.json` (ask in the report), anything in `.seq/`.

**Stop when** the checks are green, or a check has failed twice after a fix, or the work needs a public API past the
agreed surface, or a target cannot be met without changing a design (report the options), or a permission is refused
(report it; do not work around it).

**Report**: the commits; each check with its exit code and test count; the hot-file changes you need; the changelog
entry's text; config items with their defaults; drift; followups. You cannot see a browser: say what should be looked
at.
