# Lane PG: a real PCG in @molstar/core

**Goal** (the author, 2026-11-19: "Fix PCG in core"; [base-library.md](../../../design/base-library.md), "a fixed
seeded generator"): `PCG` in `packages/core/src/data/util/hash-functions.ts:588` keeps 32 bits of state but uses
PCG's 64-bit output step (xorshift by 18 and 27, rotate by the top 5 bits): its output is biased (CS measured letters
G 33%, C 20% over 1 Mb) and consecutive draws are correlated. Make it a real PCG32 (XSH RR over 64-bit state, with
`BigInt`-free 64-bit arithmetic on two 32-bit halves, or a documented equivalent), keeping its public shape
(`new PCG(seed)`, `int()`, `float()`, any other method it has). **The agreed surface** is that class's output; nothing
else in `@molstar/core`.

**Tests** (fail on the parent): the first outputs for a seed against the PCG32 reference values (pcg-random.org's
`pcg32` for the same seed and stream); floats in eighths within 1% over 1M draws; consecutive draws uncorrelated
(|r| < 0.01); speed no worse than 2x the old (a quick-lane representative, the bound heavy under `itWithin`).

**Then regenerate what depends on it, once, as diffs to review**: Mol*'s `util/color/distinct.ts` (say what changes for
a user: the distinct palettes' colors; its tests' expected values), `packages/core`'s seeded tests
(`data/immutable/_test/model.ts`, `data/int/_test/half-open-ranges.test.ts`, `repr/codec/_test/fuzz.test.ts`), Seq*'s
seeded fixtures and tests (`packages/seq/core/src/**/_test` users of `PCG`, the examples' generators: `tp53-gene.seq`
with `SEQ_RECORDS=update`, the globins scale, the plasmid, the chromosome's gene count range in
`examples/src/_test/chromosome.test.ts`). A test whose expected value changes: update it only where the new value is
as valid (a seeded sample), and list each in the report; a test that now fails for another reason: stop and report.

**Not in scope**: anything else. Runs alone in `packages/core` and the seeded files (the orchestrator starts it when
CF has merged).

**Model**: strong (opus). **Base**: `packages-seq` at the commit the orchestrator gives.

**Checks**: `pnpm test:unit packages/core packages/seq/`, `pnpm test:medium packages/seq/`, `npx tsc -b packages/core
packages/seq/core packages/seq/viewer packages/seq/examples`, the typecheck:test scripts, `pnpm lint`, `pnpm
format:check`, then `pnpm build` and `pnpm check:workspace` (the known `glass-wordmark` failure is master's).

**Headers**, **hot files** and **stop when**: as every lane's (see [qt-quick-lane.md](qt-quick-lane.md)).

**Report**: commits; checks; the reference values matched; the statistics before and after; every regenerated file and
test with what changed; the text for `.v6/plans/breaking-v6-changes.md` and a `.v6/changelogs/` entry.
