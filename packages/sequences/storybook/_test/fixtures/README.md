# The Mol\* demo's fixtures

Recorded by hand (tests never fetch); the `Seq*/Cases/Mol*` and `Seq*/File formats` stories, their smoke entries and the
demo's tests (`_test/molstar`) read them. The structures are the Mol\* adapter's, in Prism
(`packages/prism/viewer/src/sequence/molstar/_test/fixtures`).

- `p53/`: P04637's UniProt entry (`P04637.json`, https://rest.uniprot.org/uniprotkb/P04637.json; UniProtKB CC BY 4.0)
  and its Proteins API variants (`variation-P04637.json`, https://www.ebi.ac.uk/proteins/api/variation/P04637; CC BY
  4.0), recorded 2026-10-12 by `scripts/record-fixtures.ts` and trimmed to what the format's readers take: the entry's
  features but its natural variants, mutagenesis and isoforms (type, location, description); each variant inside the
  sequence (type, begin, end, wild type, alternative sequence, consequence). 17 KB and 353 KB. For the `Seq*/Mol*` demo.
  `P04637-ptms.tsv` is derived from `P04637.json` by the same script (no download): each modified residue and
  cross-link's 1-based position, its modification `kind` (phosphorylation, acetylation, methylation, ubiquitination,
  sumoylation, glycosylation, other: read from UniProt's description) and its description, 36 sites.
