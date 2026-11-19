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
//   node --experimental-strip-types scripts/record.ts p53
//
// The one entry now: `p53`, UniProt P04637's REST entry (CC-BY-4.0) and the Proteins API's variants of it, those of
// UniProt's own and ClinVar's sources.

import { execFileSync } from 'node:child_process';
import { createHash } from 'node:crypto';
import { mkdirSync, readFileSync, writeFileSync } from 'node:fs';
import { dirname, join } from 'node:path';
import { fileURLToPath } from 'node:url';

type Json = null | boolean | number | string | Json[] | { [key: string]: Json };
type Obj = { [key: string]: Json };

const Folder = join(dirname(fileURLToPath(import.meta.url)), '..');

/** The licences a fixture may carry (examples.md, "Licences allowed"). */
export const AllowedLicences = [
  'CC0-1.0',
  'CC-BY-4.0',
  'LicenseRef-NCBI-public-data',
  'LicenseRef-seqx-generated',
] as const;

/** The entry's members the views read; every other (comments, references, cross references…) is removed. */
const EntryMembers = ['entryType', 'primaryAccession', 'uniProtkbId', 'organism', 'sequence', 'features'];
/** A feature's members that are evidence, not what a view reads. */
const FeatureDropped = ['evidences', 'featureCrossReferences', 'ligand'];
/** Natural variants are the variants fixture's (a thousand rows); every other feature type stays. */
const EntryDroppedTypes = ['Natural variant'];

export function trimEntry(entry: Obj): Obj {
  const out: Obj = {};
  for (const k of EntryMembers) if (k in entry) out[k] = entry[k];
  out.features = (entry.features as Obj[])
    .filter((f) => !EntryDroppedTypes.includes(f.type as string))
    .map((f) => Object.fromEntries(Object.entries(f).filter(([k]) => !FeatureDropped.includes(k))));
  return out;
}

/** The sources of the variants kept: UniProt's own (CC-BY-4.0) and ClinVar's (no restriction). */
const KeptSources = ['UniProt', 'ClinVar'];
const VariantMembers = [
  'type',
  'alternativeSequence',
  'begin',
  'end',
  'xrefs',
  'consequenceType',
  'wildType',
  'mutatedType',
  'sourceType',
  'clinicalSignificances',
];
const ResponseMembers = [
  'accession',
  'entryName',
  'proteinName',
  'geneName',
  'organismName',
  'taxid',
  'sequence',
  'features',
];

export function trimVariation(response: Obj): Obj {
  const out: Obj = {};
  for (const k of ResponseMembers) if (k in response) out[k] = response[k];
  const features: Obj[] = [];
  // the canonical sequence's positions: the response also holds the changes of other transcripts' proteins
  const length = (response.sequence as string).length;
  for (const f of response.features as Obj[]) {
    const xrefs = ((f.xrefs ?? []) as Obj[]).filter((x) => KeptSources.includes(x.name as string));
    // the views read single-residue changes of a source of known licence
    const at = Number(f.begin);
    if (!xrefs.length || f.begin !== f.end || !f.wildType || !f.mutatedType || at < 1 || at > length) continue;
    const kept: Obj = {};
    for (const k of VariantMembers) if (k in f) kept[k] = f[k];
    // one xref a source, its name and id (the links are not data the views read)
    const first = new Map<string, Obj>();
    for (const x of xrefs)
      if (!first.has(x.name as string)) first.set(x.name as string, { name: x.name as string, id: x.id as string });
    kept.xrefs = [...first.values()];
    features.push(kept);
  }
  out.features = features;
  return out;
}

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
  // checked in as the repository's formatter writes JSON, so a format check leaves them as they are; the sizes and the
  // checksums are of the files as they stand
  const files = [join(dir, 'P04637.json'), join(dir, 'P04637-variation.json')];
  writeFileSync(files[0], entry);
  writeFileSync(files[1], variation);
  execFileSync('pnpm', ['exec', 'biome', 'format', '--write', ...files], {
    cwd: join(Folder, '..', '..', '..'),
    stdio: 'ignore',
  });
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
      licence: ['CC-BY-4.0 (UniProt’s records)', 'LicenseRef-NCBI-public-data (ClinVar’s records)'],
      attribution: 'UniProt Consortium, Proteins API, CC BY 4.0; NCBI ClinVar, no restriction on use',
      trimmed: `single-residue changes inside the canonical sequence (begin equal to end, from 1 to its length, with their wild type and mutated type) with an xref of ${KeptSources.join(' or ')}, one xref a source (name and id) of those two only; members kept: ${VariantMembers.join(', ')}; COSMIC, TCGA, gnomAD, dbSNP, Ensembl and ClinGen records and xrefs dropped`,
      derived: 'none',
      bytesDownloaded: bytes(variationText),
      bytesFixture: bytes(written[1]),
      sha256Download: sha256(variationText),
      sha256Fixture: sha256(written[1]),
    },
  ];
}

async function main(name: string | undefined) {
  if (name !== 'p53') throw new Error('usage: record.ts p53');
  const rows = await p53();
  const file = join(Folder, 'fixtures', 'p53', 'sources.json');
  let existing: Row[] = [];
  try {
    existing = JSON.parse(readFileSync(file, 'utf8')) as Row[];
  } catch {
    // none yet
  }
  const merged = [...existing.filter((e) => !rows.some((r) => r.file === e.file)), ...rows];
  for (const r of merged) console.log(`${r.file}: ${r.bytesDownloaded} bytes downloaded, ${r.bytesFixture} trimmed`);
  writeFileSync(file, `${JSON.stringify(merged, null, 2)}\n`);
  execFileSync('pnpm', ['exec', 'biome', 'format', '--write', file], {
    cwd: join(Folder, '..', '..', '..'),
    stdio: 'ignore',
  });
}

if (process.argv[1] && fileURLToPath(import.meta.url) === process.argv[1]) await main(process.argv[2]);
