# Lane SP: seq data ports (strong, design only)

**Task**: draft `drafts/seq-data-ports.md` (shape: [design-drafts.md](../../third-slice/briefs/design-drafts.md)):
Seq\*'s nodes read from a host's in-memory representation (sequences, entities, components, units) and give
annotations back as the host's, as registry entries by representation kind ([ports.md](../../../design/ports.md)'s
rule); the layering against `@molstar/io` and `@molstar/model` (a package of its own or not); whether MC's converter
is shared; a fake host's adapter.

**The author's direction (2026-11-19)**: "i suppose we can provide default ports which 'read data independently' so we
have the same mechanism?": Seq\*'s own reading (its format runtimes over the assets port, `files.md`, `sources.md`)
becomes the **default adapter of the same port**, so built-in reading and a host's data take one path: a host
replaces the adapter, nothing else changes. Design to that: what the port's surface is so both fit (reading by region,
streaming, waiting, issues, identity of nodes across a host's updates), what the default adapter wraps today's runtimes
into, and what a host adapter (the fake host; later the new mol viewer) supplies. Say what changes for today's
`sources.md` protocol and store (ideally nothing but a name) and what does not fit the one mechanism, as gaps with
recommendations.

**Reads**: ports.md, sources.md, files.md, molecules.md (gap 11, Tensions), maps.md, viewer-composition.md,
extensibility.md; the revised structure syncing draft in `../drafts/structure-syncing.md` if DR has written it (DR
runs beside you: do not edit its files); core `src/data/` and `src/resolve/` for the protocol and store as built.

**Output**: `.seq/execution/fourth-slice/drafts/seq-data-ports.md` only. No code, no commits.

**Model**: strong (opus). **Where**: the main checkout, read only except that file.

**Report**: under 300 words: the decision in brief, the port's surface, the gaps with recommendations.
