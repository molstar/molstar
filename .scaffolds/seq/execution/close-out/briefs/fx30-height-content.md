# Lane FX30: `height: 'content'`, a track or section sized to its content with no cap

**Model**: Sonnet (one param value, decided). **Base**: `packages-seq` at the orchestrator's head (`git reset --hard
packages-seq` first). Branch `seq/fx30`. Beside you: FX28 (a bridge's canvas ribbons: `core/src/view/project.ts`,
`core/src/render/bridge*.ts`) and FX29 (the legend: `viewer/src/legends.ts`, `legend-eye-dom.ts`, `mount-legend.ts`,
`core/src/render/legends.ts`, headers' parts). Never edit their files; one you need: stop and report.

**Why**: the author, 2026-11-24 (Features › Tooltips › By Data: the canvas's Transcripts track scrolls at 240 px while
the SVG draws every lane): "maybe we should support 'auto' sizing for rendered tracks on canvas too?", then "Add a spec
param". Decided with them: one value `'content'` where `height` already is.

**The decision** (write it into the designs and specs as written):
- A track's `height` takes a px number (as now) or `'content'`: the track fits its content (the Track heights rule's
  lanes, steps and fixed heights) with no `track.maxPx` cap, so its lanes never scroll inside it; `rows` bands too (no
  `rows.maxPx`). A ring's track keeps `rings.trackMaxPx` (a ring has no page to grow).
- A section's `height` takes a px number or `'content'`: the section is its content with no `section.maxPx`, never
  scrolling inside (Section scroll).
- Absent `height` is unchanged (fits up to the caps). The canvas's own height already follows the layout
  (`canvas.height: 'content'`, capped by `canvas.maxHeight`): unchanged; a section past `canvas.maxHeight` scrolls the
  canvas as now.
- Language `0.11.0` (unreleased: no new bump), its version line; the check refuses any other string
  (`repr.schema` as for other enums); the builder takes it; an older reader without it reads the track as having no
  `height` (EC-244's lenient rule: a warning, its default sizing).
- Why: the cap protects a reader's page from a section's data, which is the host's default; a figure's author who wants
  every lane shown (an export-like page, a story twin) says so in the spec.

**Read first**: the core spec's Track heights, Section scroll and Versions; `.scaffolds/seq/design/layouts.md` (tracks,
heights, section scroll) and `format.md` (versions); `core/src/render/sizing.ts`, the tracks layout, the section
layout's cap; `lang/language.ts` (the `height` params; new code in its own file if `language.ts` is past 800 lines).

**Items**, a commit at each:

1. **The language**: `height: number | 'content'` on `track` and `section`, the check, the builder, the version line.
   Proof: language tests (accepted, a wrong string refused, a 0.10 spec unchanged).
2. **The layout**: `'content'` lifts the track's and the rows band's cap and the section's cap; a pan never resizes
   (the content rule unchanged); the SVG export the same as before. Proof: layout tests (a track of 40 lanes at
   `'content'` as tall as its lanes, no inner scroll; a section past 600 px at `'content'` unscrolled); a canvas
   against SVG height check on the By Data story's composition.
3. **The story**: Features › Tooltips › By Data's Transcripts track (and its section if needed) at `height: 'content'`,
   so canvas and SVG match; its twin. Proof: the twin.
4. **The records**: `layouts.md` (Track heights, Section scroll) and its Storybook page, the core spec's Track heights,
   Section scroll, Versions and Language, the Configuration notes. Proof: `designs-in-specs` green.

**How**: nothing per frame; module size; host neutral.

**Where**: a worktree on `seq/fx30`. Headers per the root `AGENTS.md`. Hot files never staged: `package.json`s,
`tsconfig*.json`, the lockfile, `biome.json`, `.v6/changelogs/`, `storybook/packages.ts`, `packages/seq/TODO.md`,
`packages/seq/_measure/sizes.mjs`, the work folder. Never a bare `git stash`. No network. Do not run the heavy lane.

**Checks**, green only on exit 0 with a non-zero test count: `pnpm test:unit packages/seq/ cli/seqx-render/`;
`pnpm test:medium packages/seq/`; `pnpm -s -C packages/seq/core test:perf`; `pnpm -s -C packages/seq/core size` (every
row: spec 264.3 of 265, headless 992.1 of 995, default 1,584.7 of 1,588, minimal 1,117.1 of 1,120); `npx tsc -b`;
`pnpm lint`; `pnpm format:check`; the scaffold checks in `.scaffolds/seq/AGENTS.md` (stage first).

**Stop when** green, or a check red twice after a fix, or a step is refused by the permission system (report it;
don't work around it).

**Report** (under 300 words): commits; each item with its proof; the checks with exit codes and counts; the size rows;
what to look at; the language's changelog line; the changelog text.
