# Lane HC: the hover card and the keys panel

**Goal** (the author, 2026-11-19: "Yes, lane HC now"; the novice pass's first and second findings,
[novice-findings.md](../novice-findings.md) rows 1, 2, 3, 5, 7, 8, 10, 18): build what
[ux-foundations](../../../design/ux-foundations.md) §9 now says ("The hover card", "The keys panel") and
[viewer-models](../../../design/viewer-models.md)'s behaviors list names (`hoverCard`).

1. **The hover card** (viewer): a `hoverCard` behavior and a mount overlay built from `@molstar/ui` (a surface, its
   tokens and squircle; `packages/ui`'s components, no new private styles).
   - **Text**: the announcer's description of the place, one source for both (viewer `announcer.ts`, the layers'
     `ItemDescription`): move or share what the announcer builds so the card and the live region read the same words.
     The position in every numbering of its system (ordinal, index, declared: NB's), the letter, the item (label,
     extent, category), a value track's value, a pair's two ends, a row's label, a grid cell's two positions and value
     (the grid's description gains its value if missing), a block's chain and copy.
   - **A dragged rectangle on a grid**: the card shows `gridSummary` (count, mean, min, max) while dragged.
   - **Placement**: beside the pointer, never over the hovered cell or item (above, else below), kept inside the view;
     after the hover rests (the announcer's rest); hidden while a drag pans or zooms, on leave, on Esc; at the keyboard
     cursor when it moves. One card. DOM over the canvas: not in the draw lists, not in exports.
   - **Config**: on in the default spec, off by a config item; the minimal viewer leaves it out (its size measured).
     Sections behind an internal interface (extensibility.md), not exported.
2. **The keys panel**: `?` and a toolbar button ("Keys and gestures (?)") open a `@molstar/ui` panel listing the
   keymap's commands and the gesture bindings, generated, grouped (move, zoom, select, edit, view), each with its key or
   gesture and its name; Esc closes it.
3. **Stories**: the existing example views show the card with no story code (it is the viewer's); remove the stories'
   own hover readout panels where the card now says the same (Antibody's, Complex's readout: keep a panel only where
   it shows more than the card), and say so in their descriptions. A Features story for the keys panel.

**Tests** (fail on the parent): the card's text for a letter in two numberings, an item, a value, a grid cell, a pair,
a row, a block copy (headless, on the HoverModel, no DOM); equal to the announcer's text; placement never over the
hovered box and inside the view (jsdom); hidden on drag and Esc; the rectangle's summary; the keys panel lists every
bound command once with its key; config off removes it. `test:perf`: the hover pass and pointer move allocate nothing
new (the card updates on rest, not per move).

**Not in scope**: the context menu; the NP lane's fixes (it runs beside you: Storybook stories and the core's tree,
ruler, ring and theme looks; stay out of those files).

**Model**: strong (opus). **Base**: `packages-seq` at the commit this brief is in: `git checkout -B seq/hc-hover-card
<that commit>`. **Storybook port**: 6018 (your own server; stop it after).

**Checks**: `pnpm test:unit packages/seq/` (folder durations: core, viewer and storybook sit near the 2 s guide; add no
slow quick test: DOM tests are quick only if small); `pnpm test:medium packages/seq/`; `npx tsc -b packages/seq/core
packages/seq/viewer packages/seq/examples`; `pnpm -C packages/seq/<core|viewer|examples> typecheck:test`; `pnpm -C
packages/seq/storybook typecheck`; `pnpm -C packages/seq/core test:perf`; `pnpm lint`; `pnpm format:check`; `pnpm -C
packages/seq/core size` before and after. Look at the card on every dataset's main view in light and dark.

**Headers**, **hot files** and **stop when** (a public API past `hoverCard`'s config item and the keys panel's command,
or a `@molstar/ui` change): as every lane's (see [qt-quick-lane.md](qt-quick-lane.md)).

**Report**: commits; checks; folder durations; sizes; what you saw; drift; followups.
