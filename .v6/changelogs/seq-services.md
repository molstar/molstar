## Seq* (next): services in Storybook

- Seq* storybook: Examples/Services composes UniProt, PDBe SIFTS and AlphaFold DB from recorded fixtures with a live
  switch through the asset policy's hosts; a lab's predictions onto a structure, addressing beside a fake structure
  viewer (hover and selection both ways through the spec's per-chain sequences), and a view of your own file (seq,
  FASTA, Stockholm, records, mmCIF; GenBank refused, saying so; a file checked by size and content before it is read).
- Seq* examples: `services/{endpoints,sifts,composed,predictions,addressing}`; SIFTS read into structure maps keyed by
  `struct_asym_id`; `afdbConfidence` shared.
