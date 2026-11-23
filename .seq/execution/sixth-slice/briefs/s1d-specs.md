# Lane S1D: the spec pass for the ring and the catalogs

**Model**: Sonnet. Shared instructions: [s1-spec-pass.md](s1-spec-pass.md) (read it first). Branch `seq/s1d`.

**Specs you edit**: ext-rondo (and design-only rows for the catalogs).

**Designs, in order**: `ring-extension`, `examples`, `edge-cases`, `sources`, `testing`, `packages`.

**Beside you**: S1C (the core spec). For the catalogs (examples, edge-cases, sources, testing, packages) most rows are design-only (a catalog, a plan): say so in one row per group, and move a rule the code follows (e.g. testing's lanes) only if no spec has it; a rule for the whole repository stays where the root AGENTS.md has it.
