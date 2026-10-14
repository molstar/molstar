# The Mol* adapter's test fixtures

Recorded by hand (tests never fetch); the `Seq*/Mol*` Storybook demo reads them too.

- `structures/`: the structure fixtures, recorded 2026-10-12 by
  `scripts/record-fixtures.ts` (by hand only; tests never fetch). PDBe's BinaryCIF, which carries
  the updated mmCIF's SIFTS categories (wwPDB CC0 1.0, SIFTS via PDBe): `1tup.bcif` (p53 core domain with DNA, chains
  C–E, 342 KB), `2ocj.bcif` (p53 core domain, four chains, 326 KB), `2xwr.bcif` (the DNA-binding domain with an
  extended N terminus, P04637 89–293, 256 KB), `1c26.bcif` (the tetramerization domain, 325–356, chain A, the tetramer
  by its assembly, 130 KB), from `https://www.ebi.ac.uk/pdbe/entry-files/download/<id>.bcif`; and
  `AF-P04637-F1-model_v6.bcif` (AlphaFold DB's P04637 model, pLDDT in the B-factors and `ma_qa_metric_local`, CC BY
  4.0, 172 KB) from `https://alphafold.ebi.ac.uk/files/`. 1.2 MB together.
- `structures/ins1.cif`: synthetic, written by hand (no source, no license): the mmCIF text of insertion-codes.test.ts,
  one chain H of nine residues numbered 50, 51, 52, 52A, 52B, 53, 54, 55, 55A (the last unobserved), SIFTS onto a
  made-up accession Q0INS1's 70–78, a CA atom per observed residue. 2.5 KB. For the `Seq*/File formats` stories.
- `p53/`: P04637's UniProt entry (`P04637.json`, https://rest.uniprot.org/uniprotkb/P04637.json; UniProtKB CC BY 4.0)
  and its Proteins API variants (`variation-P04637.json`, https://www.ebi.ac.uk/proteins/api/variation/P04637; CC BY
  4.0), recorded 2026-10-12 by the same script and trimmed to what the format's readers take: the entry's features but
  its natural variants, mutagenesis and isoforms (type, location, description); each variant inside the sequence
  (type, begin, end, wild type, alternative sequence, consequence). 17 KB and 353 KB. For the `Seq*/Mol*` demo.
  `P04637-ptms.tsv` is derived from `P04637.json` by the same script (no download): each modified residue and
  cross-link's 1-based position, its modification `kind` (phosphorylation, acetylation, methylation, ubiquitination,
  sumoylation, glycosylation, other: read from UniProt's description) and its description, 36 sites.
