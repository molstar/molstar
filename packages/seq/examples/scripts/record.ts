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
// (`globins` imports `@molstar/core` and `@molstar/seqx-core` built: `npx tsc -b packages/seq/core` first.) The
// datasets: `p53`, UniProt P04637's REST entry (CC-BY-4.0) and the Proteins API's variants of it, those of UniProt's
// own and ClinVar's sources; `globins`, the Pfam PF00042 seed (CC0-1.0) cut to 60 rows, UniProt's organisms of its
// rows (CC-BY-4.0) and PDB 4HHB (CC0-1.0).

import { createHash } from 'node:crypto';
import { mkdirSync, readFileSync, writeFileSync } from 'node:fs';
import { gunzipSync } from 'node:zlib';
import { dirname, join } from 'node:path';
import { fileURLToPath } from 'node:url';
import {
  AllowedLicences,
  EntryDroppedTypes,
  EntryMembers,
  FeatureDropped,
  KeptSources,
  trimEntry,
  trimVariation,
  VariantMembers,
} from '../src/files/record-trim.ts';
import { HemeContacts } from '../src/globins/hemes.ts';
import {
  ContactDistance,
  cifChains,
  columnsFromUniprot,
  cutSeed,
  hemeContacts,
  hemoglobinSeq,
  SeedRows,
  seedAccessions,
  seedVersion,
  UniprotAttribution,
} from '../src/globins/recording.ts';

const Folder = join(dirname(fileURLToPath(import.meta.url)), '..');

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
  return new TextDecoder().decode(await fetchBytes(url, accept, cached));
}

async function fetchBytes(url: string, accept: string, cached: string): Promise<Uint8Array> {
  const cache = process.env.SEQ_RECORD_CACHE;
  if (cache) return new Uint8Array(readFileSync(join(cache, cached)));
  const r = await fetch(url, { headers: { accept } });
  if (!r.ok) throw new Error(`${url}: ${r.status}`);
  return new Uint8Array(await r.arrayBuffer());
}

const bytesOf = (s: string) => new TextEncoder().encode(s).byteLength;

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

/**
 * The Globins: the Pfam seed (InterPro's, gzipped Stockholm) cut to the rows nearest 4HHB's α and β chains; the
 * UniProt entries of the rows kept (their organisms and classes, a columns file joined by row key); 4HHB's mmCIF as
 * `seq` (its chains' letters and helices, the heme a position) and its heme contacts checked against `hemes.ts`.
 */
async function globins(): Promise<Row[]> {
  const dir = join(Folder, 'fixtures', 'globins');
  mkdirSync(dir, { recursive: true });
  const day = new Date().toISOString().slice(0, 10);
  const seedUrl = 'https://www.ebi.ac.uk/interpro/api/entry/pfam/PF00042/?annotation=alignment:seed';
  const cifUrl = 'https://files.rcsb.org/download/4HHB.cif';
  const seedGz = await fetchBytes(seedUrl, '*/*', 'seed.raw');
  // gzipped as served; fetch inflates it when the response says its encoding: the download's size and hash are always
  // of the inflated text, so they do not depend on which
  const inflated = seedGz[0] === 0x1f && seedGz[1] === 0x8b ? gunzipSync(seedGz) : seedGz;
  const seed = new TextDecoder().decode(inflated);
  const version = seedVersion(seed);
  if (!version) throw new Error('the seed has no #=GF AC line');
  const cif = await fetchText(cifUrl, 'text/plain', '4hhb.cif');
  const chains = cifChains(cif);
  const human = ['A', 'B'].map((a) => chains.find((c) => c.asym === a)!.letters);
  const cut = cutSeed(seed, human);
  const accessions = [...seedAccessions(cut.text).values()];
  const uniprotUrl = `https://rest.uniprot.org/uniprotkb/accessions?accessions=${accessions.join(',')}&format=tsv&fields=accession,id,organism_name,length,lineage&size=500`;
  const uniprot = await fetchText(uniprotUrl, 'text/plain', 'uniprot.tsv');
  // the contacts the views draw are code (derived): the entry must still give them
  const contacts = hemeContacts(cif);
  for (const [name, asym] of [
    ['hba', 'A'],
    ['hbb', 'B'],
  ] as const)
    if (JSON.stringify(contacts.get(asym)) !== JSON.stringify(HemeContacts[name]))
      throw new Error(`4HHB chain ${asym}'s heme contacts are ${contacts.get(asym)}, not hemes.ts's`);
  const { encodeSeq } = await import('@molstar/seqx-core/data');
  const files = {
    seed: join(dir, 'PF00042-seed.sto'),
    columns: join(dir, 'family-columns.tsv'),
    hb: join(dir, '4hhb.seq'),
  };
  writeFileSync(files.seed, cut.text);
  writeFileSync(files.columns, columnsFromUniprot(uniprot, cut.text));
  writeFileSync(files.hb, `${encodeSeq(hemoglobinSeq(cif))}\n`);
  const written = Object.fromEntries(Object.entries(files).map(([k, f]) => [k, readFileSync(f, 'utf8')]));
  return [
    {
      file: 'globins/PF00042-seed.sto',
      entry: `Pfam ${version} (Globin) seed, as InterPro serves it`,
      url: seedUrl,
      recorded: day,
      licence: ['CC0-1.0'],
      attribution: 'Pfam (EMBL-EBI), PF00042 seed alignment, CC0',
      trimmed: `cut from ${cut.kept.length + cut.removed.length} rows to ${SeedRows}: the rows nearest 4HHB's α or β chain by identity kept (HBA_HUMAN and HBB_HUMAN are not rows of this seed); removed with their #=GS and #=GR lines: ${cut.removed.join(', ')}; #=GF SQ said again (${cut.kept.length}); every other line as published; the download's size and SHA-256 are of the inflated text`,
      derived:
        'the tree over its rows is computed in code (neighbor joining on identity, src/globins/tree.ts) and written inline in the views’ specs',
      bytesDownloaded: inflated.byteLength,
      bytesFixture: bytesOf(written.seed),
      sha256Download: sha256(inflated),
      sha256Fixture: sha256(written.seed),
    },
    {
      file: 'globins/family-columns.tsv',
      entry: `UniProtKB, the ${accessions.length} entries of the seed's rows kept`,
      url: uniprotUrl,
      recorded: day,
      licence: ['CC-BY-4.0'],
      attribution: UniprotAttribution,
      trimmed: 'columns kept: Organism (as OS) and Length (as uniprot_length); the lineage reduced to its class',
      derived:
        'id: the seed row whose #=GS AC names the entry; class: the lineage’s name of rank class (empty where none)',
      bytesDownloaded: bytesOf(uniprot),
      bytesFixture: bytesOf(written.columns),
      sha256Download: sha256(uniprot),
      sha256Fixture: sha256(written.columns),
    },
    {
      file: 'globins/4hhb.seq',
      entry: 'PDB 4HHB, human deoxyhaemoglobin (Fermi, Perutz, Shaanan, Fourme 1984), mmCIF',
      url: cifUrl,
      recorded: day,
      licence: ['CC0-1.0'],
      attribution: 'wwPDB, PDB 4HHB, CC0',
      trimmed:
        'chains A (α) and B (β) only, C and D their copies: their letters (pdbx_poly_seq_scheme) and helices (struct_conf, label_seq_id); the heme a position; coordinates, waters and phosphates dropped',
      derived: `the hemes' contacting residues (heavy atoms within ${ContactDistance} Å, the heme's own chain) are code (src/globins/hemes.ts), checked against the entry here`,
      bytesDownloaded: bytesOf(cif),
      bytesFixture: bytesOf(written.hb),
      sha256Download: sha256(cif),
      sha256Fixture: sha256(written.hb),
    },
  ];
}

const Datasets: Readonly<Record<string, () => Promise<Row[]>>> = { p53, globins };

async function main(name: string | undefined) {
  if (!name || !Datasets[name]) throw new Error(`usage: record.ts ${Object.keys(Datasets).join(' | ')}`);
  const rows = await Datasets[name]();
  const file = join(Folder, 'fixtures', name, 'sources.json');
  let existing: Row[] = [];
  try {
    existing = JSON.parse(readFileSync(file, 'utf8')) as Row[];
  } catch {
    // none yet
  }
  for (const r of rows)
    for (const l of r.licence)
      if (!(AllowedLicences as readonly string[]).includes(l))
        throw new Error(`${r.file}: the licence ${l} is refused`);
  const merged = [...existing.filter((e) => !rows.some((r) => r.file === e.file)), ...rows];
  for (const r of merged) console.log(`${r.file}: ${r.bytesDownloaded} bytes downloaded, ${r.bytesFixture} trimmed`);
  writeFileSync(file, `${JSON.stringify(merged, null, 2)}\n`);
}

if (process.argv[1] && fileURLToPath(import.meta.url) === process.argv[1]) await main(process.argv[2]);
