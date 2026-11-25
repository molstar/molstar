# Lane CL6: track groups

**Model**: strong (a design's core: the language, the resolved view, the tracks layout). **Base**: `packages-seq` at
the orchestrator's head (`git reset --hard packages-seq` first). Branch `seq/cl6`. Beside you: CL5 (search: a new
`core/src/search/`, `viewer/src/models/search*`, the highlight channels in `viewer/src/models/` highlight and legend
code, `storybook/src/_stories/support/search-box.ts`, a plasmid story) and FX27 (the ring's wheel:
`viewer/src/behaviors/navigate.ts`, `viewer/src/input.ts`, `ext-rondo/`; `core/src/data/` intervals). Never edit their
files; the specs (`_spec/*.md`) and `viewer/src/config.ts` are shared: add to your own sections and entries only, never
rewrite another's lines. A file of theirs you need otherwise: stop and report. FX25's parked worktree holds
`core/src/render/` ranges and lanes and the variant card: the ranges layer's `lanes: one` (Decision 6) needs
`ranges.ts` or `lanes.ts`: make that item last, and stop and report before editing either.

**Why**: features 159 (groups and summaries), 70 and 141 (groups of tracks, EC-103), 185 (collapse and groups) and
128 (collapse), version 1's rows (`drafts/version-1-rows.md`); the design agreed by the author, 2026-11-24 ("Agree all":
the language change and the kind's name `group` as written). **This is the language lane**: it makes the bump to
`0.11.0`; CL8 adds the values layer's `reference` param within the same unreleased `0.11.0` after you.

**Read first**: `.scaffolds/seq/execution/close-out/drafts/track-groups.md` (every section, its Defaults,
Configuration, Consequences and "In the specs"); `.scaffolds/seq/design/format.md` (versions, validation),
`layouts.md` (tracks, heights), `track-header-parts.md`, `ux-foundations.md` (the header cursor, the matrix),
`viewer-models.md` (history kinds); EC-103 and EC-244 in `.scaffolds/seq/sources/seqstar/edge-cases.md`; the core and
viewer specs' sections the design's "In the specs" names.

**Items**, a commit at each:

0. **The design in place**: `drafts/track-groups.md` copied to `.scaffolds/seq/design/track-groups.md` (header comment
   kept; `§n` self-citations written as "Decision n"), listed in `.scaffolds/seq/design/README.md`; its Storybook page
   under `Seq* (next)/Design` (`design/track-groups.mdx`). Proof: `principles-ids` and `links` green.
1. **The `group` kind** (Decision 1): its params, `track.parents` gaining `group`, one level (`seq.check.parent`),
   `seq.check.group.summary`, the builder's `section.group`, the version `0.11.0` with its version history line; a
   `0.10` spec reads unchanged; the version-aware validate lists `group` as a kind an older reader flattens (EC-244).
   Proof: language tests per check; a `0.10` fixture round trip.
2. **The resolved view** (Decision 2): tracks flat in spec order, `ResolvedTrack.group`, `ResolvedSection.trackGroups`.
   Proof: resolve tests.
3. **Folding as session state** (Decisions 3, 8): the folded set in the section's viewport, `collapsed` the opening
   state, a toggle a hand's move, bulk and reveal jumps, reset on its patch; what a fold does to the cursor, the
   selection and the heights. Proof: viewer model tests on a fake clock, history steps counted.
4. **The header and the summary row** (Decisions 4, 5): the header line, indent and gaps, the hit target; the summary
   band, each layer kind's summary form (`LayerType.collapse`), `summary: ref` shown only while folded (EC-103),
   `none`. Proof: layout tests per kind; a test named for EC-103.
5. **Hidden, order, editing and commands** (Decisions 7, 9): `groups.toggle`, `collapseAll`, `expandAll`, `keep`;
   `track.hide`/`show`/`move` on a group's address; the header cursor's group stop, Space and Enter; the menu's
   `on: 'group'` entries; hidden groups in the host's list; groups open in text and rings
   (`seq.render.groups.layout`); `renderSvg`'s `groups`. Proof: a test per command and per key.
6. **The ranges layer's `lanes: one`** (Decision 6) and `layer.lanes`. Proof: tests on a ranges track in one lane with
   its labels off, folded and open. (See Beside you: stop before `ranges.ts`/`lanes.ts`.)
7. **The stories**: the p53 protein section in three groups (the design's Defaults), a story per fold state, each with
   its headless twin on a fake clock. Proof: the twins.
8. **The specs and records**: every row of "In the specs" in its section; the Configuration tables. Proof:
   `designs-in-specs` green with the design in `design/`.

**How**: nothing per frame (a fold changes the layout once, never a warm frame's work); host neutral; reuse the
clade folds' model (Decision 3); module size (§module-size: `language.ts`, `tracks.ts`, `view.ts` are over 800 lines
already: new code in files of its own); new code lazy where the default viewer can do without it at load.

**Where**: a worktree on `seq/cl6`. Headers per the root `AGENTS.md`. Hot files never staged: `package.json`s,
`tsconfig*.json`, the lockfile, `biome.json`, `.v6/changelogs/`, `storybook/packages.ts`, `packages/seq/TODO.md`,
`packages/seq/_measure/sizes.mjs`, the work folder (name what you need in the report). Never a bare `git stash`. No
network. Do not run the heavy lane (another lane measures beside you).

**Checks**, green only on exit 0 with a non-zero test count: `pnpm test:unit packages/seq/ cli/seqx-render/`;
`pnpm test:medium packages/seq/`; `pnpm -s -C packages/seq/core test:perf` and `pnpm -s -C packages/seq/storybook
test:perf`; `pnpm -s -C packages/seq/core size` (every row: spec 260.9 of 265, headless 978.7 of 980, default 1,557.3
of 1,560, minimal 1,096.6 of 1,100: name a row that needs a raise and by how much); `npx tsc -b`; `pnpm lint` (5
warnings before); `pnpm format:check`; the scaffold checks in `.scaffolds/seq/AGENTS.md` (stage first).

**Stop when** green, or a check red twice after a fix, or a step is refused by the permission system (report the
action and why; don't work around it, not even by another tool or the Edit tool).

**Report** (under 450 words): commits; each item with its proof; the checks with exit codes and counts; the size rows;
what to look at (p53's groups folded and open, the header cursor, light and dark); drift, gaps and followups; the
hot-file changes needed; the language's changelog line; the changelog text.
