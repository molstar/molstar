/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Records @molstar/seq-molstar-poc's fixtures (packages/sequences/molstar-poc/_test/fixtures/structures and p53; VENDORED.md, seqstar). By hand only,
// never in tests: it downloads public files (wwPDB CC0 via PDBe; AlphaFold DB, UniProtKB and the Proteins API CC BY 4.0).
// P04637's UniProt entry and its Proteins API variants are trimmed to what the readers take (the format's `features`):
// the entry's features but its natural variants, mutagenesis and isoforms; each variant's position, letters and
// consequence (no description, so the consequence is the row's label).
// Usage: node packages/sequences/molstar-poc/scripts/record-fixtures.ts [file name…] (P04637-ptms.tsv only when named)
import { mkdirSync, readFileSync, writeFileSync } from 'node:fs';
import { join } from 'node:path';

const OUT = join(import.meta.dirname, '../_test/fixtures/structures');
const ONLY = process.argv.slice(2);
const FILES: Record<string, string> = {
  // PDBe's BinaryCIF carries the updated mmCIF's SIFTS categories (pdbx_sifts_xref_db, pdbx_sifts_unp_segments).
  '1tup.bcif': 'https://www.ebi.ac.uk/pdbe/entry-files/download/1tup.bcif',
  '2ocj.bcif': 'https://www.ebi.ac.uk/pdbe/entry-files/download/2ocj.bcif',
  // 2XWR is the DNA-binding domain with an extended N terminus (P04637 89-293), not the tetramerization domain.
  '2xwr.bcif': 'https://www.ebi.ac.uk/pdbe/entry-files/download/2xwr.bcif',
  // The tetramerization domain (P04637 325-356): chain A deposited, the tetramer by its assembly.
  '1c26.bcif': 'https://www.ebi.ac.uk/pdbe/entry-files/download/1c26.bcif',
  'AF-P04637-F1-model_v6.bcif': 'https://alphafold.ebi.ac.uk/files/AF-P04637-F1-model_v6.bcif',
};

const P53 = join(import.meta.dirname, '../_test/fixtures/p53');
type F = Record<string, unknown> & { type?: string };
const TRIMS: Record<string, [string, (v: { features: F[] }) => unknown]> = {
  'P04637.json': [
    'https://rest.uniprot.org/uniprotkb/P04637.json',
    (v) => ({
      primaryAccession: 'P04637',
      features: v.features
        .filter((f) => !['Natural variant', 'Mutagenesis', 'Alternative sequence'].includes(f.type ?? ''))
        .map(({ type, location, description }) => ({ type, location, description })),
    }),
  ],
  'variation-P04637.json': [
    'https://www.ebi.ac.uk/proteins/api/variation/P04637',
    (v) => ({
      accession: 'P04637',
      // Variants inside the sequence (a frameshift's reach past its end is not a position of P04637).
      features: v.features
        .filter((f) => Number(f.end) <= 393 && Number(f.begin) >= 1)
        .map(({ type, begin, end, wildType, alternativeSequence, consequenceType }) => ({
          type,
          begin,
          end,
          wildType,
          alternativeSequence,
          consequenceType,
        })),
    }),
  ],
};
mkdirSync(P53, { recursive: true });
for (const [name, [url, trim]] of Object.entries(TRIMS)) {
  if (ONLY.length && !ONLY.includes(name)) continue;
  const res = await fetch(url, { headers: { Accept: 'application/json' } });
  if (!res.ok) throw new Error(`${url}: ${res.status}`);
  const text = JSON.stringify(trim((await res.json()) as { features: F[] }));
  writeFileSync(join(P53, name), text);
  console.log(`${name}\t${text.length} bytes\t${url}`);
}

/**
 * P04637's PTM sites with their kind, derived from the recorded entry (no download): a TSV of `pos` (1-based), `kind`
 * and `label`, for a track colored by modification type. The kind is read from UniProt's description here, in the
 * fixture's conversion, never in code that draws or colors.
 */
const PTMS = 'P04637-ptms.tsv';
function ptmKind(type: string, description: string): string {
  const d = description.toLowerCase();
  if (type === 'Glycosylation') return 'glycosylation';
  if (type === 'Cross-link')
    return d.includes('sumo') ? 'sumoylation' : d.includes('ubiquitin') ? 'ubiquitination' : 'other';
  if (d.startsWith('phospho')) return 'phosphorylation';
  if (d.includes('acetyl')) return 'acetylation';
  if (d.includes('methyl')) return 'methylation';
  return 'other';
}
if (ONLY.includes(PTMS)) {
  const entry = JSON.parse(readFileSync(join(P53, 'P04637.json'), 'utf8')) as {
    features: { type: string; description?: string; location: { start: { value: number } } }[];
  };
  const rows = entry.features
    .filter((f) => ['Modified residue', 'Cross-link', 'Glycosylation'].includes(f.type))
    .map((f) => `${f.location.start.value}\t${ptmKind(f.type, f.description ?? '')}\t${f.description ?? f.type}`);
  writeFileSync(join(P53, PTMS), `pos\tkind\tlabel\n${rows.join('\n')}\n`);
  console.log(`${PTMS}\t${rows.length} sites, from P04637.json`);
}

mkdirSync(OUT, { recursive: true });
for (const [name, url] of Object.entries(FILES)) {
  if (ONLY.length && !ONLY.includes(name)) continue;
  const res = await fetch(url);
  if (!res.ok) throw new Error(`${url}: ${res.status}`);
  const bytes = new Uint8Array(await res.arrayBuffer());
  writeFileSync(join(OUT, name), bytes);
  console.log(`${name}\t${bytes.length} bytes\t${url}`);
}
