# Lane AY: accessibility

**Goal**: the plan's AY: **titles in exported SVG, the cursor and hover announced, forced colors**, as
[ux-foundations](../../../design/ux-foundations.md) (accessibility) and the viewer `_spec` (Input, announcements)
say, plus what today's lanes added: keyboard reach to the row columns' titles (RC: a title orders the rows by a click
only) and the highlight legend's Clear buttons (HL), the query box (QS) already labelled.
1. **Exported SVG**: a `<title>` and `<desc>` on the SVG (the spec's title, the view's description), a `<title>` on
   each track, each item (its label and extent) and each highlight; `role="img"` with `aria-labelledby`.
2. **Announcements**: the keyboard cursor's position (system, position, letter, the items there) and the hover's
   (when it rests) through one polite live region the mount owns, rate-limited; a selection's change announced once
   (its size and extent). Through the viewer's existing announcement path if one exists; else a small one in the
   mount, the `_spec` updated.
3. **Forced colors** (`forced-colors: active`): the canvas reads the system colors (`CanvasText`, `Canvas`,
   `Highlight`, `GrayText`) into its tokens when the media query matches, so the data's ink, the selection and the hover
   stay visible; tested by stubbing the media query.
4. **Keyboard to the row-column titles**: a command `rows.order` reachable from the keyboard (a key on the focused rows
   layer cycles the order fields, or the title is focusable in the DOM overlay if one exists): say which and why.

**Tests** (fail on the parent): the SVG's titles (axe-like checks on the exported SVG: every graphic labelled), the
live region's text after cursor moves and a rested hover, forced colors' tokens, keyboard ordering.

**Not in scope**: anything else. Other lanes run beside you; stay out of the datasets' builders.

**Model**: fast (sonnet). **Base**: `packages-seq` at the commit this brief is in: `git checkout -B seq/ay-accessibility
<that commit>`.

**Checks**: `pnpm test:unit packages/seq/`; `pnpm test:medium packages/seq/`; `npx tsc -b packages/seq/core
packages/seq/viewer packages/seq/examples`; `pnpm -C packages/seq/<core|viewer> typecheck:test`; `pnpm -C
packages/seq/storybook typecheck`; `pnpm lint`; `pnpm format:check`; `pnpm -C packages/seq/core test:perf`. Exit codes
explicit, non-zero test counts.

**Headers**, **hot files**, **stop when** and **report**: as every lane's (see `briefs/op-png-cli.md`'s), no new
package.
