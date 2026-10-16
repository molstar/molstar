# The Mol\* adapter's test fixtures

Recorded by hand (tests never fetch). The Prism sequence tests, Model Archive's Prism test and the Seq\* Storybook's
Mol\* demo (`packages/sequences/storybook`, its stories and their tests) read them too.

- `structures/`: the structure fixtures, recorded 2026-10-12 by
  `packages/sequences/storybook/scripts/record-fixtures.ts` (by hand only; tests never fetch). PDBe's BinaryCIF, which
  carries the updated mmCIF's SIFTS categories (wwPDB CC0 1.0, SIFTS via PDBe): `1tup.bcif` (p53 core domain with DNA,
  chains C–E, 342 KB), `2ocj.bcif` (p53 core domain, four chains, 326 KB), `2xwr.bcif` (the DNA-binding domain with an
  extended N terminus, P04637 89–293, 256 KB), `1c26.bcif` (the tetramerization domain, 325–356, chain A, the tetramer
  by its assembly, 130 KB), from `https://www.ebi.ac.uk/pdbe/entry-files/download/<id>.bcif`; and
  `AF-P04637-F1-model_v6.bcif` (AlphaFold DB's P04637 model, pLDDT in the B-factors and `ma_qa_metric_local`, CC BY 4.0,
  172 KB) from `https://alphafold.ebi.ac.uk/files/`. 1.2 MB together.
- `structures/ins1.cif`: synthetic, written by hand (no source, no license): the mmCIF text of insertion-codes.test.ts,
  one chain H of nine residues numbered 50, 51, 52, 52A, 52B, 53, 54, 55, 55A (the last unobserved), SIFTS onto a
  made-up accession Q0INS1's 70–78, a CA atom per observed residue. 2.5 KB. For the `Seq*/File formats` stories too.
- The PAE test reads Seq\*'s `af-barnase-barstar-pae.json` and the load-path tests P04637's sequence (`p53.json`) from
  `packages/sequences/core/_test/fixtures`, as data.
