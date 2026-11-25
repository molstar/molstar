# Lane FS-A1: the records runtime and the core additions

**Goal**: build what [files.md](../../../design/files.md) (agreed 2026-11-19, "Ok as recommended") lists for FS-A1 in its
Consequences, Lanes:

1. **MessagePack moved down** (base-library.md's row): `@molstar/io`'s `common/msgpack` (encoder, decoder) to
   `@molstar/core/data/msgpack`, the tree's codec (`io/repr/msgpack.ts`) to `@molstar/core/repr/codecs/msgpack` beside
   the JSON and S-expression codecs; `@molstar/io` re-exports both at their old paths (BinaryCIF, g3d and the base
   sequences unchanged, their tests green); the codec's column extensions exported as a decoder hook for a plain value.
2. **A delimited-record tokenizer** at `@molstar/core/data/text/delimited` (files.md, The file's own syntax): RFC 4180,
   plain string cells, one-character delimiter, quote and comment, a BOM skipped, LF/CRLF/CR; io's CSV parser the
   guide, io's reader unchanged.
3. **`seqx-core`**: `AssetReader.readBytes` (and the policy adapter's, and `node-assets.ts` in examples);
   `OpenContext.limits` (protocol.md, Runtimes); `data/records.ts`: the runtime for TSV, CSV, JSON and MessagePack in
   both layouts with the transform params of the `source` node (`fields`, `records`, `letters`, `system`, `base`,
   `ends`, `types`, `attribution`, `licence`) in language 0.4.0, its `seq.source.records.*` issues and bounds; the
   `seq` runtime reading MessagePack; the default spec registering the records runtime and the viewer's source sync
   passing `limits`; the `secondary` class table reading `Helix`, `Beta strand`, `Turn`.
4. **Tests** as files.md's Consequences list for seqx-core (the conformance suite; one table in every encoding as a
   property; the transform's property; typed SoA columns kept by identity; each issue; fuzz without a throw; bounds
   before parsing); the 10,000-row open quick, 100,000 under 300 ms heavy (`itWithin`).
5. **`_spec`s**: seqx-core's Data ("Records"), Protocol, Language, Configuration; core's and io's where they describe
   MessagePack.

**Not in scope**: GFF3, `dataFilter`'s trailing rule, p53's entry, the Files stories (FS-A2); held trees and pairs
(FS-B). **Agreed surface**: the two core modules and io's re-exports, the decoder hook, `readBytes`, `limits`, the
source params above. Anything past it stops the lane.

**Model**: strong (opus). **Base**: `packages-seq` at the commit this brief is in: `git checkout -B seq/fs-a1-records
<that commit>`. **The language stays 0.4.0.** Nothing else runs in `packages/seq` or `packages/core/src/data` beside
you.

**Minimal code**: reuse what exists in `@molstar/core` (AGENTS.md, Reuse before writing helpers); no private copies.

**Checks**: `pnpm test:unit packages/seq/ packages/core/ packages/io/`; `pnpm test:medium packages/seq/`; `npx tsc -b`;
`pnpm -C packages/seq/<core|viewer|examples> typecheck:test`; `pnpm -C packages/seq/storybook typecheck`; `pnpm lint`;
`pnpm format:check`; `pnpm check:workspace` (after the orchestrator applies the exports); `pnpm -C packages/seq/core
size` before and after. Green only on exit 0 with a non-zero test count.

**Headers**: a new file has the Mol\* header `Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file
for more info.` with `@author David Sehnal <david.sehnal@gmail.com>`; a moved file keeps its header and authors, year to
2026, that author added; a modified file likewise. No mention of `.scaffolds/seq/` in code; code cites its own `_spec`.

**Hot files, never staged**: `pnpm-lock.yaml`, `scripts/workspace/inventory.json`, the root `tsconfig.json`,
`biome.json`, `vitest.config.ts`, `storybook/packages.ts`, `.v6/changelogs/` (write the entries' text in the report),
every `package.json` and `tsconfig*.json` (write the exact export entries for `@molstar/core` and `@molstar/io` in the
report; for local runs, edit them unstaged), anything in `.scaffolds/seq/`.

**Stop when** green, or a check fails twice after a fix, or the work needs a public API past the agreed surface, or a
permission is refused.

**Report**: commits; each check with exit code and test count; sizes before and after; the hot-file entries; drift;
followups.
