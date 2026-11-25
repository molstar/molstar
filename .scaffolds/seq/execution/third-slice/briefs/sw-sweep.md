# Lane SW: today's small follow-ups

From the lanes merged on 2026-11-19 (`.scaffolds/seq/execution/third-slice/followups.md`), each small, each with a test that fails
on its parent:
1. **A pan keeps the content under the pointer** (LB): `SectionNav.pan` gets no pointer x, so a row of blocks pans
   about the view's middle and the content under the pointer drifts when a chip or a gap enters the view. Pass the
   pointer through `view.panPx` (`packages/seq/viewer/src/behaviors/navigate.ts`, `models/view.ts`) to the core's
   `rowPan`.
2. **A ring chip's label once** (LB): a chip on a ring shows its label outside and in its pill; drop the outside label
   when the pill draws it.
3. **A minimal viewer's size** (the author, 2026-11-19: "i think a minimal one could be smaller"; plan stage 7): a
   `viewer, minimal` row in `packages/seq/_measure/sizes.mjs` (a viewer composed with tracks, letters and ranges only:
   no navigators, stories, editing, records runtime, rows, grids), and find what keeps the unused parts from shaking
   out (side-effectful module tops, registries that import every layer type): list them with their sizes; fix the
   cheap ones (a `/*#__PURE__*/`, a lazy registry entry) where no public API changes, report the rest.
Not now: `rows.ts`'s private `elide` (DG is in rows).

**Model**: fast (sonnet). **Base**: `packages-seq` at the commit this brief is in: `git checkout -B seq/sw-sweep <that
commit>`.

**Checks**: `pnpm test:unit packages/seq/`; `pnpm test:medium packages/seq/`; `npx tsc -b packages/seq/core
packages/seq/viewer packages/seq/examples`; `pnpm -C packages/seq/<core|viewer> typecheck:test`; `pnpm -C
packages/seq/storybook typecheck`; `pnpm lint`; `pnpm format:check`; `pnpm -C packages/seq/core test:perf`; `pnpm -C
packages/seq/core size` before and after. Exit codes explicit, non-zero test counts.

**Headers**, **hot files**, **stop when** and **report**: as every lane's (see `briefs/op-png-cli.md`'s), no new
package; the report lists the sizes before and after and what keeps the minimal viewer large.
