# Maps: kind, sourcing, encoding; structures

## Intent

Maps carry data from one coordinate system to another: a protein to its gene's codons, a row to its alignment's
columns, a design to its wild type, a sequence to a structure's residues, one genome to another. They must be small to
write, cheap to hold, and the same mechanism whatever their size or origin.

## Decision

The author, 2026-11-17: "could map be computational or compressed so we dont always spell every position out?"; then
"the "computed" … should still be a codon mapping .. the "computed" is how it's sourced compared to a table… a
"computed" one mapped to a runtime during "viewer composition" is also important. also we should have a "compressed
version" of the mappings too … not always we need full expanded table"; and "for the data maps, we need ability to map
to structural data .. basically SIFTS, but in abstract form."

- **Three questions, answered separately** ([building-blocks.md](building-blocks.md), Map):
  - **kind**, what a map means: `identity`, `offset`, `numbering`, `segments`, `pairwise`, `alignment`, `codon`,
    `structure`;
  - **sourcing**, where it comes from: inline (the kind's compact form), a source (a table or file through the source
    protocol, by region), or computed (a runtime the application binds when it composes the viewer). Computed is never
    a kind: a codon map computed from a gene model is a `codon` map.
  - **encoding**, how a resolved map is held: implicit (arithmetic), runs (segments), or expanded (per position, for
    the region asked), always the cheapest that answers the map interface.
- **One interface**: every resolved map answers the image of a range forward and back as runs, the boundary rules and
  the status of each element (paired, conflict, unpaired). Paths compose encodings without expanding them; projection
  reads only the interface.
- **Structure maps, SIFTS in abstract form**:
  - A structure system is a chain instance (entry, model, assembly operator, chain); its positions are its residues in
    the structure's order, each with a key (sequence index where the structure has one; author number and insertion
    code) and an observed flag.
  - A `structure` map relates a sequence (UniProt, or a construct) to that chain: segments, each with its provenance
    (SIFTS, an alignment, the author's numbering), and per-residue status: paired, conflict (an engineered mutation, a
    microheterogeneity), unpaired (an expression tag, a linker), unobserved.
  - A chain may map to several sequences (a chimera, a fusion), and a sequence to many chains (copies, entries); each
    chain instance stays apart by the copy rule.
  - Sourcing as for any map: inline (small, for tests and examples), a source (a SIFTS-style residue table), or
    computed (a runtime that looks SIFTS up or aligns the chain); the binding to atoms is the structure runtime's (Mol\*'s
    adapter), so the core reads no structure file.

## Consequences

- A variant on a gene reaches a residue of a structure through a codon map and a structure map, with no special case.
- The p53 dataset carries a structure map with an expression tag, unobserved loops and author numbering with an offset
  ([examples.md](examples.md)); its conformance cases include a chimeric chain and copies.
- Edge cases EC-3 to EC-37 apply to every encoding: a test runs each on implicit, runs and expanded forms of the same map
  and expects the same answers.
