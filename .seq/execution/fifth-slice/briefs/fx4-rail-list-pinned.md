# Lane FX4: the rail's place as a list, and a pinned place

**Base**: `packages-seq` at the orchestrator's head (`git reset --hard packages-seq` first). Branch `seq/fx-4`. Model:
Sonnet. Follow `fx-fix-lane.md` (one commit per item, a test failing on its parent, a look in light and dark). Read
`.seq/design/navigators.md` ("The positions rail beside a declared navigator, and its place" and "A list, and a pinned
place", and the Configuration row), core's and the viewer's `_spec` (Navigators, Configuration), and FX3's code
(`railRows`, `railSections`, `RailPlace`, `RailRow`) first.

**Items**:

1. **`navigator.railPlace` takes one place or a list** (`RailPlace | RailPlace[]`): each place a rail on the one
   window (a drag on any moves the others in the same frame); the fallbacks first (`navigator` without a horizontal
   navigator → `ruler`; `ruler` without a ruler → `bottom`), then duplicates dropped. Twin tests: `['ruler', 'bottom']`
   draws two rails; a drag on the bottom one moves the top one's window in the same frame; `['ruler', 'ruler']` and a
   fallback onto a listed place give one rail.
2. **`pinned`**: a 12 px rail row pinned to the bottom edge of the view's visible frame, outside the scrolled rows (as
   the ruler is pinned at the top): it stays in view at any vertical scroll and takes no row in an export. Twin tests:
   scrolled to the top and to the bottom, the pinned rail's row is at the frame's bottom both times; an SVG export has
   no rail row. A story beside FamilyRailAtBottom (`FamilyRailPinned`, with `['ruler', 'pinned']`) on the Globins rows
   family.

Keep the default `ruler`; nothing moves for a reader who sets nothing. The viewer and core `_spec` follow the design's
text.

**Where**: a worktree on `seq/fx-4`. Headers updated (`@author David Sehnal <david.sehnal@gmail.com>`, 2026); hot files
never staged (`.seq/`, `.v6/changelogs/`, `packages/seq/TODO.md`, `package.json`s): list the edits. Do not start
storybook servers you leave running; stop any you start before you report.

**Checks** (all, with numbers): `pnpm test:unit packages/seq/`, `npx tsc -b`, `pnpm -s -C
packages/seq/{core,viewer,examples} typecheck:test`, `pnpm -C packages/seq/storybook typecheck`, `pnpm lint`,
`pnpm format:check`, both perf tiers, `pnpm -s -C packages/seq/core size`, `pnpm test:heavy packages/seq/`, the four
`.seq/_checks`.

**Report**: under 200 words: commits; checks with numbers; the story ids; the `.v6/changelogs/` text.
