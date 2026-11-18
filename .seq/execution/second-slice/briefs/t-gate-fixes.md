# Lane T: stage 4's gate findings (multi-range selection, followers that lead, the hover at far zoom)

**Goal**, from the author's use of the p53 stories (2026-11-18):
1. **Several continuous ranges in one selection** (the author: "i am missing the ability to select multiple continuous
   ranges. shift + click always extends the selection"; approved as proposed, bindings configurable):
   - drag replaces the selection with the range (as now);
   - Cmd/Ctrl+drag **adds** the range, keeping the rest;
   - Alt/Option+drag **removes** the range from the selection;
   - Shift+click and Shift+drag extend **only the range the anchor is in**, keeping the others (today an extend
     replaces everything with anchor..element);
   - Cmd/Ctrl+click toggles the element (as now), the anchor moving there;
   - Shift+arrows extend the anchor's range; Cmd/Ctrl+Shift+arrows start a new range at the cursor;
   - commands for each (`selection.add`, `selection.subtract` over a range or an address), so the keyboard and a host
     reach them; Esc during an additive or subtractive drag restores the selection it began with.
   The selection's Location already holds several ranges a system: this is the gestures and the model's operations.
2. **A follower leads too** (the author: "i can't zoom ... or slide on the TP53 track .. so the binding is just one
   way"): a pan, zoom or pinch on a follower moves its leader through the follow's inverse, so the follower moves as
   asked and every section in step. `follow.followerMoves` gains `lead` as the default (`ignore` and `detach` kept).
   Where the inverse is not one run (the window over an intron, the ends crossing) the leader takes the hull of the
   preimage; a minus-strand follower's pan maps to the leader's opposite direction (what the follower shows moves the
   way the hand moves). Focus, reveal and the cursor on a follower keep their rules (viewer _spec "Follows").
3. **The hover blinks at far zoom** (the author: on the TP53 track at the fully zoomed-out view the hover box "blinks
   on and off" as the pointer moves; zoomed in it does not). Find the cause (hit-testing at many positions a pixel,
   the projected hover's minimum width, a hover cleared between frames, a coalescing race), fix it with a test shown
   failing on the parent: moving across the gene section at the fit view hovers continuously, in both sections.

Not in scope: a reversed follower (third slice), protein letters on the gene's axis (third slice), anything in
`packages/core`, the story steps.

**Read first**: `.seq/design/ux-foundations.md` (sections 6 to 8, the gesture table as now updated), `maps.md` (follow),
`viewer-models.md`; the viewer and core `_spec`s; `packages/seq/viewer/src/{models/selection.ts,models/view.ts,input.ts,bindings.ts,commands.ts}`; `packages/seq/AGENTS.md`.

**Model**: strong (opus).

**Base**: `packages-seq` at `BASE`. Your worktree may be cut from master: with a clean tree, first
`git checkout -B seq/t-gate BASE`, and check `git merge-base --is-ancestor BASE HEAD`.

**Touches**: `packages/seq/viewer/src/**`, `packages/seq/core/src/{render,query,view}/**` where hover or hit-testing
needs it, `packages/seq/storybook/**` (the Features/Selection story's hint names the new gestures; the protein-on-gene
stories pan and zoom from either section), the `_spec`s (the viewer's gesture and follow rules, Configuration).

**Two consumers**: a multi-range selection made in the viewer drawn in an SVG export of the same spec, and projected
onto the gene section in both.

**Tests**: headless twins on a fake clock for add, subtract, extend-one-of-several and Esc; a pan and a zoom on the
gene follower moving the protein leader (both strands, an intron in the window); the far-zoom hover.

**Performance representative**: a follower pan step and an additive drag step under 4 ms (`itWithin`, several after
warm-ups); hover warm allocation-free (the opt-in tier).

**Configurable over fixed**: the bindings (add, subtract, extend scope), `follow.followerMoves` (`lead` default), in the
`_spec`s' Configuration tables.

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

**Hot files, never staged**:
- `pnpm-lock.yaml`, `scripts/workspace/inventory.json`, the root `tsconfig.json`, `biome.json`, `vitest.config.ts`;
- `storybook/packages.ts`, `.v6/changelogs/` (write the changelog entry's text in your report);
- every `package.json` and `tsconfig*.json` (ask in the report);
- anything in `.seq/`.

**Stop when** the checks are green, or a check has failed twice after a fix, or the work needs a public API past the
agreed surface, or a target cannot be met without changing a design (report the measurement and the options), or a
permission is refused (report the action and why; do not work around it).

**Report**: the commits; each check with its exit code and test count; the hot-file changes you need; the changelog
entry's text; config items with their defaults; drift (what, where, options, your recommendation); followups. You
cannot see a browser: say what should be looked at.
