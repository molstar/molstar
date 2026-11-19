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
// The entries: `p53`, UniProt P04637's REST entry (CC-BY-4.0) and the Proteins API's variants of it, those of
// UniProt's own and ClinVar's sources; `plasmid`, INSDC L09137.2 (pUC19c) as ENA serves it in FASTA, converted to the
// `seq` format, its letters also as a module of the examples (a cyclic sequence read from a source has no length until
// it is opened, so a view that wraps the origin takes them inline), and the rotated variant's features and sites,
// derived from it, as GFF3.

import { createHash } from 'node:crypto';
import { mkdirSync, readFileSync, writeFileSync } from 'node:fs';
import { dirname, join } from 'node:path';
import { fileURLToPath } from 'node:url';
import {
  EntryDroppedTypes,
  EntryMembers,
  FeatureDropped,
  KeptSources,
  trimEntry,
  trimVariation,
  VariantMembers,
} from '../src/files/record-trim.ts';
import {
  gff3Text,
  PlasmidAccession,
  PlasmidFiles,
  PlasmidLength,
  PlasmidName,
  recordModule,
  rotatedAnnotations,
  seqText,
} from '../src/plasmid/derive.ts';

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

/** The fixtures of the Plasmid: the record converted, and the files derived from it by `src/plasmid/derive.ts`. */
async function plasmid(): Promise<Row[]> {
  const dir = join(Folder, 'fixtures', 'plasmid');
  mkdirSync(dir, { recursive: true });
  const day = new Date().toISOString().slice(0, 10);
  const url = `https://www.ebi.ac.uk/ena/browser/api/fasta/${PlasmidAccession}`;
  const text = await fetchText(url, 'text/plain', `${PlasmidAccession}.fa`);
  const lines = text.split('\n');
  const letters = lines
    .filter((l) => !l.startsWith('>'))
    .join('')
    .trim()
    .toUpperCase();
  if (letters.length !== PlasmidLength || /[^ACGT]/.test(letters))
    throw new Error(`${PlasmidAccession}: not 2,686 bp ACGT`);
  const rotated = rotatedAnnotations(letters);
  const bytes = (s: string) => new TextEncoder().encode(s).byteLength;
  const generated = 'LicenseRef-seqx-generated';
  const write = (path: string, content: string, row: Partial<Row>): Row => {
    writeFileSync(join(Folder, path), content);
    return {
      file: path.replace('fixtures/', ''),
      entry: '',
      url: 'none',
      recorded: day,
      licence: [generated],
      attribution: 'generated by the examples package',
      trimmed: 'none',
      derived: '',
      bytesDownloaded: 0,
      bytesFixture: bytes(content),
      sha256Download: '',
      sha256Fixture: sha256(content),
      ...row,
    };
  };
  const published = seqText(PlasmidName, letters);
  writeFileSync(join(Folder, 'src', 'plasmid', 'record.ts'), recordModule(letters));
  const derived = (what: string) => `${what}, computed in code from the record’s letters (src/plasmid/derive.ts)`;
  return [
    write(PlasmidFiles.published, published, {
      entry: `INSDC ${PlasmidAccession}, Cloning vector pUC19c, complete sequence (2,686 bp, circular)`,
      url,
      licence: ['LicenseRef-NCBI-public-data'],
      attribution: `INSDC (ENA, NCBI) ${PlasmidAccession}; Yanisch-Perron C, Vieira J, Messing J. Gene 33:103-119 (1985)`,
      trimmed:
        'the FASTA record’s letters only, converted to the seq format because FASTA has no cyclic flag: one cyclic system named pUC19c (the record has no feature annotation but source)',
      derived: 'none',
      bytesDownloaded: bytes(text),
      sha256Download: sha256(text),
    }),
    write(PlasmidFiles.rotatedFeatures, gff3Text(rotated.features, rotated.sites), {
      entry:
        'the features and sites of the rotated variant (origin at 2000, inside bla), a generated primer site across it',
      derived: derived(
        'bla and lacZα as open reading frames, the cloning site from EcoRI to HindIII, the origin by its signature, every site of the listed enzymes, all turned to the new origin; the primer site is generated',
      ),
    }),
  ];
}

async function main(name: string | undefined) {
  if (name !== 'p53' && name !== 'plasmid') throw new Error('usage: record.ts p53 | plasmid');
  const rows = await (name === 'p53' ? p53() : plasmid());
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
