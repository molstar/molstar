# Lane CL5: search

**Model**: strong (a design's core: matchers over untrusted patterns, a regex automaton, a viewer model and its public
surface). **Base**: `packages-seq` at the orchestrator's head (`git reset --hard packages-seq` first). Branch `seq/cl5`.
Beside you: CL6 (track groups: the language, `core/src/lang/`, `core/src/resolve/`, the tracks layout and headers in
`core/src/render/`, `viewer/src/models/view.ts`, the headers' behaviors, `examples/src/p53` and its views) and FX27 (the
ring's wheel: `viewer/src/behaviors/navigate.ts`, `viewer/src/input.ts`, `ext-rondo/`; `core/src/data/` intervals).
Never edit their files; the specs (`_spec/*.md`) and `viewer/src/config.ts` are shared: add to your own sections and
entries only, never rewrite another's lines. A file of theirs you need otherwise: stop and report. FX25's parked
worktree holds `core/src/render/` ranges and lanes and the variant card: avoid those too.

**Why**: features 182 (search and results) and the search parts of 180, 181 and 207, version 1's rows
(`drafts/version-1-rows.md`); the design agreed by the author, 2026-11-24 ("Agree all"; the API it adds included).

**Read first**: `.scaffolds/seq/execution/close-out/drafts/search.md` (the design: every section, its Defaults,
Configuration, Consequences and "In the specs"); `.scaffolds/seq/design/security.md` ("regular expressions in search"),
`extensibility.md` (registry entries), `addressing.md` (refs, reveal), `ux-foundations.md` (highlight channels, the
stripes rule, keys), `viewer-models.md`; `.scaffolds/seq/sources/seqstar/edge-cases.md` for EC-10, 12, 159–162 and 169;
the previous implementation's search in `packages/sequences` (read only: the kernels' shape); the core and viewer specs'
sections the design's "In the specs" names; `packages/seq/TODO.md`'s Bounds line.

**Items**, a commit at each:

0. **The design in place**: `drafts/search.md` copied to `.scaffolds/seq/design/search.md` (its header comment kept;
   its `§n` self-citations written as "Decision n", the principles' ids only after `§`), listed in
   `.scaffolds/seq/design/README.md`; its Storybook page under `Seq* (next)/Design`
   (`packages/seq/storybook/src/_stories/design/search.mdx`). Proof: `principles-ids` and `links` green (the
   orchestrator deletes the draft at the merge).
1. **Kernels and matchers in the core** (Decisions 1–4): exact, degenerate and regex over a system's stored letters,
   sync kernels to a visitor with nothing allocated a hit, `kernelMatcher`'s chunks and overlap, hits columnar with
   strands, palindromes once, circles (EC-10) and blocks. Proof: tests per matcher and per hit rule, a hit across a
   circle's origin named for EC-10; a perf-tier case at 0 bytes a hit.
2. **The regex automaton and the bounds** (Decisions 3, 6): the Glushkov automaton, never backtracking, its refusals
   and every bound with its code, counted past the cap. Proof: a test running a classic catastrophic pattern within
   `maxSteps`; a test per bound's code.
3. **The registry entry** (Decision 7): `search-matcher` in core's registry, built-ins unnamespaced, a missing one
   `seq.search.matcher.missing`. Proof: a test composing an application's matcher and one naming a missing matcher.
4. **`SearchModel` and the commands** (Decision 5): `viewer.search`, `search.run`/`next`/`previous`/`go`/`clear`/`save`
   with their kinds; `next` from the view (EC-159), `go` toggling (EC-161) and refusing outside (EC-160), `save` as one
   ref patch, one undo step (EC-162), refused past `maxSaved` or incomplete; stale on new letters; letters through the
   store, a pending region never searched as empty (EC-169). Proof: tests named per EC, on a fake clock.
5. **The channels, keys and words** (Decision 5, Defaults): `search` and `search-current` through `HighlightModel`,
   projected onto linked sections, the legend's channel removal running `search.clear`; `Mod+g`/`F3` and their
   Shift forms in the canvas scope; what is said, progress once after 1 s. Proof: tests for each direction of the
   legend link, the keys, and the announcements.
6. **The Storybook box and stories** (Decision 8): `storybook/src/_stories/support/search-box.ts` on the public
   surface only, the p53 "search" story and a plasmid story (a site across the origin, both strands), each with its
   headless twin on a fake clock. Proof: a test holding the box's imports to public exports; the twins.
7. **The specs and records**: every row of the design's "In the specs" in its spec section; the Configuration tables;
   `security.md`'s row marked built. Proof: `designs-in-specs` green with the design in `design/`.

**How**: engine, not the game: no chrome in a package (the box is the stories'); host neutral (the keys in the canvas
scope, never `Mod+f`); nothing per frame; reuse from `@molstar/core` first (`data/util`, `data/int`, `algorithms/`);
module size (§module-size), the matchers in their own chunk if the default viewer does not need them at load. No mock
exists: build to the design's Defaults; the orchestrator looks gesture by gesture before the merge.

**Where**: a worktree on `seq/cl5`. Headers per the root `AGENTS.md`. Hot files never staged: `package.json`s,
`tsconfig*.json`, the lockfile, `biome.json`, `.v6/changelogs/`, `storybook/packages.ts`, `packages/seq/TODO.md`,
`packages/seq/_measure/sizes.mjs`, the work folder (name what you need in the report). Never a bare `git stash`. No
network. Do not run the heavy lane (another lane measures beside you).

**Checks**, green only on exit 0 with a non-zero test count: `pnpm test:unit packages/seq/ cli/seqx-render/`;
`pnpm test:medium packages/seq/`; `pnpm -s -C packages/seq/core test:perf` and `pnpm -s -C packages/seq/storybook
test:perf`; `pnpm -s -C packages/seq/core size` (every row: headless 978.7 of 980, default 1,557.3 of 1,560, minimal
1,096.6 of 1,100: new code lazy where it can be; name a row that needs a raise and by how much); `npx tsc -b`;
`pnpm lint` (5 warnings before); `pnpm format:check`; the scaffold checks in `.scaffolds/seq/AGENTS.md` (stage first).

**Stop when** green, or a check red twice after a fix, or a step is refused by the permission system (report the
action and why; don't work around it, not even by another tool or the Edit tool).

**Report** (under 450 words): commits; each item with its proof; the checks with exit codes and counts; the size rows;
what to look at (the box on p53 and the plasmid, light and dark, gesture by gesture); drift, gaps and followups; the
hot-file changes needed; the changelog text.
