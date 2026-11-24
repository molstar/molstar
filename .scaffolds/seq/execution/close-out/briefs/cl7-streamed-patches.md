# Lane CL7: patches on streamed things

**Model**: strong (a design's core: patches, refs and the resolve pass over dynamic nodes). **Base**: `packages-seq` at
the orchestrator's head (`git reset --hard packages-seq` first). Branch `seq/cl7`. Beside you: CL9 (the viewer's
selection, input, navigator cursor and keys; `core/src/lang/{check,issues,language}.ts` and the version-aware
validate; the Antibody example) and CL13 (`core/src/theme/provenance.ts`, a folder-graph test, the lanes band's history,
`rowsValue` on navigators, `seq.render.pairs.filter`). Never edit their files; if you must, stop and report.

**Why**: version 1's planned row "Patches on streamed things" (95), built for version 1 (the author, 2026-11-23,
`drafts/version-1-rows.md`): p53's transcripts already stream, and a bookmark naming one is unsafe without it.

**Read first**: `.scaffolds/seq/design/changes.md` ("Patches may name streamed things"; its In the specs row: core ›
Resolving), `addressing.md` (refs, bookmarks and sets, dangling refused, rewritten on move), `dynamic-nodes.md`
("Naming before arrival"; the overlay holds dynamic nodes); `.scaffolds/seq/design/features.md` row 95 and row 184
(User-defined refs); `.scaffolds/seq/sources/seqstar/edge-cases.md` for the dynamic nodes' items (EC-176 to EC-180) and
refs; the core spec's Resolving and the viewer spec's AddressModel and History; the code in
`packages/seq/core/src/{resolve,lang,data}/` and `packages/seq/viewer/src/models/{address,spec}.ts`.

**Items**, a commit at each:

0. **The checks first**: what a bookmark or ref naming a dynamic node does today (refused as dangling? resolved only
   while held?); report it. A contradiction with the design stops the lane: a drift note.
1. **A patch may name a streamed thing**: a bookmark or ref naming a dynamic node (`/data/genes/(node :key T1)`) is
   accepted as a patch whether or not the node is held; it is not "dangling" while its source can still serve the key.
   Proof: tests on p53's transcripts by region: bookmark a transcript, evict its region, the bookmark stays.
2. **The ref waits while the node is not held**: what names it is pending (a layer pending, a jump said "waiting for
   T1"), resolves when the node arrives, and says once if the source answers that the key does not exist. Proof: tests
   with a delayed source (fake clock), and with a key the source refuses.
3. **Undo works whether or not the node is held**: undoing a bookmark on an evicted node, redoing it. Proof: tests.
4. **The specs and designs**: core Resolving, the viewer's AddressModel; the design's "Built" note lists the ticked
   items; the features row 95 moved to built (`features.md`, `features.mdx`). Proof: `designs-in-specs` green; the rows.

**How**: host neutral (the undo stack is the host's history through the viewer); reuse the store's pending machinery;
module size (§module-size); nothing per frame (allocation tiers green); bounded work (a waiting ref holds no data); EC
ids in test titles; stories with twins where something is drawn (p53: transcripts by region).

**Where**: a worktree on `seq/cl7`. Out of scope: CL9's and CL13's files. Headers per the root `AGENTS.md`. Hot files
never staged: `package.json`s, `tsconfig*.json`, the lockfile, `biome.json`, `.v6/changelogs/`, `storybook/packages.ts`,
`packages/seq/TODO.md`, `packages/seq/_measure/sizes.mjs`, the work folder. The language's grammar does not change; a
public API change beyond the design stops the lane. Never a bare `git stash`. No network.

**Checks**, green only on exit 0 with a non-zero test count: `pnpm test:unit packages/seq/ cli/seqx-render/`; `pnpm
test:medium packages/seq/`; `pnpm -s -C packages/seq/core test:perf`; `pnpm -s -C packages/seq/core size` (every row);
`npx tsc -b`; `pnpm lint`; `pnpm format:check`; the scaffold checks in `.scaffolds/seq/AGENTS.md` (stage first).

**Stop when** green, or a check red twice after a fix, or a step is refused by the permission system (report the
action and why; don't work around it, not even by another tool).

**Report** (under 400 words): commits; each item with its proof; the checks with exit codes and counts; the size rows;
what to look at; drift, gaps and followups; the changelog text.
