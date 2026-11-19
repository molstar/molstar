/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Records the examples' published fixtures, by hand, never in tests (examples.md, "Recorded by hand"): fetches each
// entry from its service once, trims it to what the views read (members and records removed, none renamed or moved, so
// a spec's paths read the live response too), writes the fixture and its `sources.json` row (URL, date, licence,
// attribution, what was trimmed, checksums). Run it from the package folder:
//
//   npx tsc -b packages/seq/examples   (from the repository root: the trna entry imports the built `lib/`)
//   node --experimental-strip-types scripts/record.ts p53
//
// The entries: `p53`, UniProt P04637's REST entry (CC-BY-4.0) and the Proteins API's variants of it, those of
// UniProt's own and ClinVar's sources; `trna`, PDB 1EHZ's mmCIF from PDBe (CC0-1.0: chain A's sequence, modified
// residues and base pairs as `seq` files) and Rfam RF00005's seed (CC0-1.0, Stockholm as published).

import { createHash } from 'node:crypto';
import { mkdirSync, readFileSync, writeFileSync } from 'node:fs';
import { dirname, join } from 'node:path';
import { fileURLToPath } from 'node:url';
import { writeDerived } from '../lib/trna/generated.js';
import { moleculeSeq, pairsSeq, trimEntry as trimCif } from '../src/trna/record-trim.ts';
import {
  EntryDroppedTypes,
  EntryMembers,
  FeatureDropped,
  KeptSources,
  trimEntry,
  trimVariation,
  VariantMembers,
} from '../src/files/record-trim.ts';

const Folder = join(dirname(fileURLToPath(import.meta.url)), '..');

/** The licences a fixture may carry (examples.md, "Licences allowed"). */
export const AllowedLicences = [
  'CC0-1.0',
  'CC-BY-4.0',
  'LicenseRef-NCBI-public-data',
  'LicenseRef-seqx-generated',
] as const;

const sha256 = (b: Uint8Array | string) => createHash('sha256').update(b).digest('hex');

interface Row {
  readonly file: string;
  readonly entry: string;
  readonly url: string;
  readonly recorded: string;
  readonly licence: readonly string[];
  readonly attribution: string;
  readonly trimmed: string;
  readonly derived: string;
  readonly bytesDownloaded: number;
  readonly bytesFixture: number;
  readonly sha256Download: string;
  readonly sha256Fixture: string;
}

/**
 * A response as the service sent it. With `SEQ_RECORD_CACHE=<folder>` a response already downloaded (`entry.json`,
 * `variation.json`) is trimmed instead of fetched again: the services are asked once.
 */
async function fetchText(url: string, accept: string, cached: string): Promise<string> {
  const cache = process.env.SEQ_RECORD_CACHE;
  if (cache) return readFileSync(join(cache, cached), 'utf8');
  const r = await fetch(url, { headers: { accept } });
  if (!r.ok) throw new Error(`${url}: ${r.status}`);
  return r.text();
}

async function p53(): Promise<Row[]> {
  const dir = join(Folder, 'fixtures', 'p53');
  mkdirSync(dir, { recursive: true });
  const day = new Date().toISOString().slice(0, 10);
  const entryUrl = 'https://rest.uniprot.org/uniprotkb/P04637.json';
  const variationUrl = 'https://www.ebi.ac.uk/proteins/api/variation/P04637';
  const entryText = await fetchText(entryUrl, 'application/json', 'entry.json');
  const variationText = await fetchText(variationUrl, 'application/json', 'variation.json');
  const entry = JSON.stringify(trimEntry(JSON.parse(entryText)));
  const variation = JSON.stringify(trimVariation(JSON.parse(variationText)));
  // checked in compact (the formatter skips the fixtures' JSON); the sizes and the checksums are of the files as written
  const files = [join(dir, 'P04637.json'), join(dir, 'P04637-variation.json')];
  writeFileSync(files[0], `${entry}\n`);
  writeFileSync(files[1], `${variation}\n`);
  const written = files.map((f) => readFileSync(f, 'utf8'));
  const bytes = (s: string) => new TextEncoder().encode(s).byteLength;
  return [
    {
      file: 'p53/P04637.json',
      entry: 'UniProtKB P04637 (P53_HUMAN), REST entry',
      url: entryUrl,
      recorded: day,
      licence: ['CC-BY-4.0'],
      attribution: 'UniProt Consortium, UniProtKB P04637, CC BY 4.0',
      trimmed: `members kept: ${EntryMembers.join(', ')}; features without ${FeatureDropped.join(', ')} and without the ${EntryDroppedTypes.join(', ')} features (the variants file's)`,
      derived: 'none',
      bytesDownloaded: bytes(entryText),
      bytesFixture: bytes(written[0]),
      sha256Download: sha256(entryText),
      sha256Fixture: sha256(written[0]),
    },
    {
      file: 'p53/P04637-variation.json',
      entry: 'Proteins API variation, P04637',
      url: variationUrl,
      recorded: day,
      licence: ['CC-BY-4.0'],
      attribution:
        'UniProt Consortium, Proteins API, CC BY 4.0 (the variants are UniProt’s: sourceType uniprot or mixed); the xrefs to ClinVar are names only',
      trimmed: `single-residue changes of UniProt’s own (sourceType uniprot or mixed, not those only large-scale studies report) inside the canonical sequence (begin equal to end, from 1 to its length, with their wild type and mutated type) with an xref of ${KeptSources.join(' or ')}, one xref a source by its name, of those two only; each significance by its type; members kept: ${VariantMembers.join(', ')}; records only large-scale studies report and every other source’s xrefs dropped`,
      derived: 'none',
      bytesDownloaded: bytes(variationText),
      bytesFixture: bytes(written[1]),
      sha256Download: sha256(variationText),
      sha256Fixture: sha256(written[1]),
    },
  ];
}

async function trna(): Promise<Row[]> {
  const dir = join(Folder, 'fixtures', 'trna');
  mkdirSync(dir, { recursive: true });
  const day = new Date().toISOString().slice(0, 10);
  const cifUrl = 'https://www.ebi.ac.uk/pdbe/entry-files/download/1ehz.cif';
  const seedUrl = 'https://rfam.org/family/RF00005/alignment/stockholm';
  const cif = await fetchText(cifUrl, 'text/plain', '1EHZ.cif');
  const seed = await fetchText(seedUrl, 'text/plain', 'RF00005.sto');
  // the release the seed belongs to, from Rfam's own README of the current release
  const readmeUrl = 'https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/README';
  const release = /^Release\s+(\S+)/m.exec(await fetchText(readmeUrl, 'text/plain', 'rfam-README.txt'))?.[1];
  if (!release) throw new Error(`${readmeUrl}: no release line`);
  const trimmed = trimCif(cif);
  const bytes = (s: string) => new TextEncoder().encode(s).byteLength;
  const files: [string, string][] = [
    ['1EHZ.seq', moleculeSeq(trimmed)],
    ['1EHZ-pairs.seq', pairsSeq(trimmed)],
    ['RF00005.sto', seed],
  ];
  for (const [f, text] of files) writeFileSync(join(dir, f), text);
  const pdb = {
    entry: 'PDB 1EHZ (Shi, Moore 2000), yeast tRNA-Phe at 1.93 Å, chain A',
    url: cifUrl,
    licence: ['CC0-1.0'],
    attribution: 'wwPDB, entry 1EHZ (CC0 1.0), as served by PDBe',
    bytesDownloaded: bytes(cif),
    sha256Download: sha256(cif),
  };
  // what is computed from the two files, committed as code the builders read synchronously (derive.ts)
  const derivedFile = join(Folder, 'src', 'trna', 'derived-data.ts');
  writeDerived(derivedFile, seed, trimmed.letters);
  const derivedText = readFileSync(derivedFile, 'utf8');
  const generated = {
    recorded: day,
    licence: ['LicenseRef-seqx-generated'],
    url: '',
    bytesDownloaded: 0,
    sha256Download: sha256(`${files[0][1]}${seed}`),
  };
  return [
    {
      ...generated,
      file: 'src/trna/derived-data.ts',
      entry: '1EHZ’s place in the family and the covariation of its positions',
      attribution:
        'computed by src/trna/derive.ts from trna/1EHZ.seq and trna/RF00005.sto (sha256Download: the two joined)',
      trimmed: 'none',
      derived:
        'the seed row of highest identity to 1EHZ (a pairwise alignment, core’s globalAlign), 1EHZ’s positions’ columns in it, mutual information over those columns less the Miller–Madow correction (at least 30 rows, 0.05 bits), the row’s letters and its SS_cons pairs',
      bytesFixture: bytes(derivedText),
      sha256Fixture: sha256(derivedText),
    },
    {
      ...generated,
      file: 'src/trna/data.ts#scaleGrid',
      entry: 'the dot plot’s 4,000 × 4,000 scale variant',
      attribution: 'generated by src/trna/data.ts on core’s PCG (seeds 4000 and 4001), never stored',
      trimmed: 'none',
      derived: 'letters and values from a seed',
      bytesFixture: 0,
      sha256Download: '',
      sha256Fixture: '',
    },
    {
      ...pdb,
      file: 'trna/1EHZ.seq',
      recorded: day,
      trimmed:
        'chain A’s canonical sequence (parent letters, `_entity_poly`) and its modified residues (`_pdbx_struct_mod_residue`: a range each, its code as `label`, its parent letter as `category`); nothing else of the entry',
      derived: 'none',
      bytesFixture: bytes(files[0][1]),
      sha256Fixture: sha256(files[0][1]),
    },
    {
      ...pdb,
      file: 'trna/1EHZ-pairs.seq',
      recorded: day,
      trimmed:
        'the entry’s base-pair table (`ndb_struct_na_base_pair`, 27 pairs of chain A) as a range on each pair’s first base with its partner; kind from the Saenger class (19, 20 Watson–Crick; 28 wobble; else non-canonical), tertiary the six pairs outside the four stems that hold the L-shape (8–14, 15–48, 18–55, 19–56, 26–44, 54–58)',
      derived: 'kind and tertiary flag',
      bytesFixture: bytes(files[1][1]),
      sha256Fixture: sha256(files[1][1]),
    },
    {
      file: 'trna/RF00005.sto',
      entry: `Rfam RF00005 (tRNA) seed, release ${release}`,
      url: seedUrl,
      recorded: day,
      licence: ['CC0-1.0'],
      attribution: 'Rfam (EMBL-EBI), family RF00005, CC0 1.0',
      trimmed: 'none: the seed as published',
      derived: 'none',
      bytesDownloaded: bytes(seed),
      bytesFixture: bytes(seed),
      sha256Download: sha256(seed),
      sha256Fixture: sha256(seed),
    },
  ];
}

async function main(name: string | undefined) {
  if (name !== 'p53' && name !== 'trna') throw new Error('usage: record.ts p53|trna');
  const rows = await (name === 'p53' ? p53() : trna());
  const file = join(Folder, 'fixtures', name, 'sources.json');
  let existing: Row[] = [];
  try {
    existing = JSON.parse(readFileSync(file, 'utf8')) as Row[];
  } catch {
    // none yet
  }
  const merged = [...existing.filter((e) => !rows.some((r) => r.file === e.file)), ...rows];
  for (const r of merged) console.log(`${r.file}: ${r.bytesDownloaded} bytes downloaded, ${r.bytesFixture} trimmed`);
  writeFileSync(file, `${JSON.stringify(merged, null, 2)}\n`);
}

if (process.argv[1] && fileURLToPath(import.meta.url) === process.argv[1]) await main(process.argv[2]);
