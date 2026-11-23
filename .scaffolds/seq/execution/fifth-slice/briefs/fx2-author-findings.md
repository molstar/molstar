# Lane FX2: the author's findings on the genome stories

**Base**: `packages-seq` at the orchestrator's head (`git reset --hard packages-seq` first). Branch `seq/fx-2`. The last
agent task before the usage limit (the author, 2026-11-21): one lane, no reviewer after it, so the orchestrator's look
and checks are the review. Follow `fx-fix-lane.md` (one commit per item, a test failing on its parent, a look in light
and dark). Read core's and the viewer's `_spec` (Render: Ranges, Headers; viewer Models, Input) first.

**Items** (`fixes.md`, Open), in this order:

1. **The "…" actions button with no actions** (gff3--tp-53's Genes track): draw it only when at least one entry
   applies to that track; find why the fourth slice's fix of this is not holding; test every header kind.
2. **The SVG not repainted when the selection is cleared** (gff3--tp-53): every story's SVG redraws on every change
   the canvas does; a twin test: select, clear, the SVG holds no selection mark.
3. **The hover card far from the pointer** in a tall track (gff3--tp-53: hovering TP53-227 puts the card at the
   Transcripts track's top): anchor by the hovered item's row or the pointer.
4. **The Transcripts track's spacing**: blank lanes above the first and below the last transcript (the track taller
   than its lanes, lanes centred), the first row unlabelled, lanes uneven around labels: lanes from the top, evenly.
5. **A double-click on a long sequence jumps to bases** (BED cytobands, chr17): one zoom step, or fit the item under
   the pointer; never the closest zoom in one double-click (blocks are already fitted by GC3: keep that).

Leave the "bare number as a block name" item open (it needs a core repr surface).

**Where**: a worktree on `seq/fx-2`. Headers updated (`@author David Sehnal <david.sehnal@gmail.com>`, 2026); hot
files never staged (`.seq/`, `.v6/changelogs/`, `packages/seq/TODO.md`, `package.json`s): list the edits. Public API
changes allowed in `packages/seq`; anything else stops the item.

**Checks**: `pnpm test:unit packages/seq/`, `npx tsc -b`, `pnpm -s -C packages/seq/{core,viewer,examples}
typecheck:test`, `pnpm -C packages/seq/storybook typecheck`, `pnpm lint`, `pnpm format:check`, both perf tiers
(`pnpm -s -C packages/seq/{core,storybook} test:perf`), `pnpm -s -C packages/seq/core size`, the four `.seq/_checks`.

**Report**: under 250 words: commits per item; checks with numbers; the stories to look at; anything left and why;
the `.v6/changelogs/` text.
