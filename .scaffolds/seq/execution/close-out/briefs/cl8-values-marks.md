# Lane CL8: values marks (reference lines)

**Model**: Sonnet (a closed design, one layer param). **Base**: `packages-seq` at the orchestrator's head
(`git reset --hard packages-seq` first). Branch `seq/cl8`. Beside you: FX25 (TP53's lanes in `core/src/render/lanes.ts`
and `ranges.ts`, the variant card in `region-hover.ts`, `runtimes/vcf-genotypes.ts`, `data/protocol.ts`,
`tables-forms.ts`, stub labels in `ext-rondo/`, F9's pressed slot in the canvas headers' draw, the square's failure
text, the zoom hatch in the region reads). Never edit their files; the specs (`_spec/*.md`) are shared: add to your own
sections only. A file of theirs you need: stop and report.

**Why**: feature 130 ("Values for genomics": reference lines; points and the log scale are built), version 1's rows
(`drafts/version-1-rows.md`); the design agreed by the author, 2026-11-24 ("Agree all"). **No version bump**: CL6 made
`0.11.0`, unreleased; the `reference` param joins it, a line in core's Versions under 0.11.0.

**Read first**: `.scaffolds/seq/execution/close-out/drafts/values-marks.md` (every section, its Defaults,
Configuration, Consequences, Tensions and "In the specs"); `.scaffolds/seq/design/themes.md` ("Palettes and scales",
EC-215), `layouts.md` (the values look), `format.md` (versions); EC-56, 201, 215–217 in
`.scaffolds/seq/sources/seqstar/edge-cases.md`; `core/src/render/values-marks.ts` and the values layer; the core
spec's sections the design's "In the specs" names.

**Items**, a commit at each:

0. **The design in place**: `drafts/values-marks.md` copied to `.scaffolds/seq/design/values-marks.md` (its own file;
   header comment kept; `§n` self-citations written as "Decision n"), listed in `.scaffolds/seq/design/README.md`; its
   Storybook page under `Seq* (next)/Design` (`design/values-marks.mdx`). Proof: `principles-ids` and `links` green.
1. **The `reference` param** (Decision 1): up to four lines `{ at, label?, color?, style? }`, checked (the bound, `at`
   required, a color never a theme form), in `0.11.0`'s version line; the builder. Proof: language tests per check; a
   round trip.
2. **Through the layer's scale** (Decisions 2, 3): `ys` with the headroom, a log scale (`at <= 0`: no line), a closed
   domain's outside line not drawn and said once (`seq.render.reference.outside`), an open end including the lines
   (the union of EC-56). Proof: tests per case, EC-217 named.
3. **The label** (Decision 4): right end, above the line, 10 px, over the well's ground, only from 40 px of well; below
   near the top; the later of two overlapping dropped. Proof: layout tests on the label boxes.
4. **Marks, export, accessibility, cost** (Decision 5): drawn over the data and under the hover, never hit; solid a
   rect, dashed and dotted a dashed path, the SVG export the same; the accessible text "reference line at 7.3
   (genome-wide)"; a warm frame of each encoding with two lines at 0 B (the perf tier beside the encodings'). Proof:
   the draw list tests, an SVG test, the perf case.
5. **The story**: a Manhattan-like Genomes story with a line at 7.3 ("genome-wide") and one at 5 ("suggestive"), over
   data to about 6 on an open domain (the Provisional's reading), with its headless twin on a fake clock; listed in
   `examples.md`. Proof: the twin.
6. **The specs and records**: every row of "In the specs" in its section; the Configuration table; feature 130's row.
   Proof: `designs-in-specs` green with the design in `design/`.

**How**: nothing per frame; host neutral; reuse the draw lists' dash and the values layer's scale; module size
(§module-size); new code lazy where the default viewer does without it.

**Where**: a worktree on `seq/cl8`. Headers per the root `AGENTS.md`. Hot files never staged: `package.json`s,
`tsconfig*.json`, the lockfile, `biome.json`, `.v6/changelogs/`, `storybook/packages.ts`, `packages/seq/TODO.md`,
`packages/seq/_measure/sizes.mjs`, the work folder. Never a bare `git stash`. No network. Do not run the heavy lane.

**Checks**, green only on exit 0 with a non-zero test count: `pnpm test:unit packages/seq/ cli/seqx-render/`;
`pnpm test:medium packages/seq/`; `pnpm -s -C packages/seq/core test:perf` and `pnpm -s -C packages/seq/storybook
test:perf`; `pnpm -s -C packages/seq/core size` (every row: spec 263.6 of 265, headless 986.4 of 990, default 1,578.9
of 1,582, minimal 1,112.3 of 1,115: name a row that needs a raise and by how much); `npx tsc -b`; `pnpm lint` (5
warnings before); `pnpm format:check`; the scaffold checks in `.scaffolds/seq/AGENTS.md` (stage first).

**Stop when** green, or a check red twice after a fix, or a step is refused by the permission system (report the
action and why; don't work around it, not even by another tool).

**Report** (under 400 words): commits; each item with its proof; the checks with exit codes and counts; the size rows;
what to look at (the story, light and dark, at two well heights); drift, gaps and followups; the changelog text.
