## Seq*: genome language 0.6.0

- Language 0.6.0: assemblies (a genome's sequences by name, with an order, an include filter of safe patterns, derived
  sequences and per-sequence overrides), sections of assemblies, a map between assemblies, and the sequence, locus and
  block forms; the builder's genome helpers (`assembly`, `genome`, `colinearBlocks`).
- Genome resolve is loaded when used: a spec with an assembly is pending until the chunk lands, and a failed load is
  reported, not left pending. Declared params nothing reads yet are reported (`seq.resolve.param.unread`).
