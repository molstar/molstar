# Lane S1D: the spec pass for the ring and the catalogs

**Model**: Sonnet. Shared instructions: [s1-spec-pass.md](s1-spec-pass.md) (read it first). Branch `seq/s1d`.

**Specs you edit**: ext-rondo (and design-only rows for the catalogs).

**Designs, in order**: `ring-extension`, `examples`, `edge-cases`, `sources`, `testing`, `packages`.

**Beside you**: S1C (the core spec). For the catalogs (examples, edge-cases, sources, testing, packages) most rows are design-only (a catalog, a plan): say so in one row per group, and move a rule the code follows (e.g. testing's lanes) only if no spec has it; a rule for the whole repository stays where the root AGENTS.md has it.

**Also, the examples' spec** (from S1B's review, 2026-11-23): `@molstar/seqx-examples` (`packages/seq/examples`) has
no `_spec`, so rules its code follows stayed design-only: small committed fixtures, generated ones from a seed,
fixtures read by their relative location, a test per view, the builder modules per example. Write
`packages/seq/examples/_spec/README.md` (the package's spec: what it holds, its rules and their reasons, terse, in the
other specs' voice), point `examples` and `examples-and-stories` rows at it as `examples › <section>`, and list in your
report that `designs-in-specs.ts` must learn the `examples` spec (the orchestrator adds it before you start).
