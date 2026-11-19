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
// UniProt's own and ClinVar's sources; `complex`, PDB 1FIN (CC0-1.0) converted to `seq` with its derived contacts, the
// ligands of six CC0 entries of CDK2 (AlphaFill's transplants are not taken: its usage policy is no stated CC BY 4.0
// or looser licence, the data's decisions, 1), UniProt P24941 and P20248 trimmed, and AlphaFold DB's per-residue
// confidence of both chains as published (CC-BY-4.0).

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
  ContactCutoff,
  contactsSeq,
  type Ligand,
  moleculeSeq,
  trim1fin,
  trimLigand,
} from '../src/complex/record-trim.ts';

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

/** The CDK2 entries the ligands come from, with the component taken from each (CC0-1.0, RCSB). */
export const LigandEntries: readonly (readonly [entry: string, comp: string])[] = [
  ['1AQ1', 'STU'],
  ['2VTA', 'LZ1'],
  ['6Q3B', 'PYZ'],
  ['2C6O', '4SP'],
  ['3PXQ', '2AN'],
  ['3TI1', 'B49'],
];

async function complex(): Promise<Row[]> {
  const dir = join(Folder, 'fixtures', 'complex');
  mkdirSync(dir, { recursive: true });
  const day = new Date().toISOString().slice(0, 10);
  const bytes = (s: string) => new TextEncoder().encode(s).byteLength;
  const rows: Row[] = [];
  const write = (name: string, text: string) => {
    writeFileSync(join(dir, name), text);
    return readFileSync(join(dir, name), 'utf8');
  };
  const pdb = (id: string) => `https://files.rcsb.org/download/${id}.cif`;
  const fin = await fetchText(pdb('1FIN'), 'text/plain', '1FIN.cif');
  const t = trim1fin(fin);
  const ligands: Ligand[] = [];
  const sources: { id: string; text: string }[] = [];
  for (const [id, comp] of LigandEntries) {
    const text = await fetchText(pdb(id), 'text/plain', `${id}.cif`);
    sources.push({ id, text });
    ligands.push(trimLigand(text, id, comp, t.cdk2));
  }
  const molecule = write('1FIN.seq', moleculeSeq(t));
  const contacts = write('contacts.seq', contactsSeq(t, ligands));
  const cc0 = (
    file: string,
    entry: string,
    url: string,
    download: string,
    fixture: string,
    trimmed: string,
    derived: string,
  ) =>
    rows.push({
      file: `complex/${file}`,
      entry,
      url,
      recorded: day,
      licence: ['CC0-1.0'],
      attribution: `RCSB PDB ${entry}, CC0 1.0`,
      trimmed,
      derived,
      bytesDownloaded: bytes(download),
      bytesFixture: bytes(fixture),
      sha256Download: sha256(download),
      sha256Fixture: sha256(fixture),
    });
  cc0(
    '1FIN.seq',
    'PDB 1FIN (cyclin A–CDK2–ATP, Jeffrey et al. 1995)',
    pdb('1FIN'),
    fin,
    molecule,
    'converted to seq: each entity’s letters (CDK2 chains A, C; the cyclin A construct, chains B, D) and its chains’ helices and strands (struct_conf, struct_sheet_range) as ranges with chain and category; the coordinates dropped',
    'none',
  );
  cc0(
    'contacts.seq',
    `PDB 1FIN and ${LigandEntries.map(([e, c]) => `${e} (${c})`).join(', ')}`,
    [pdb('1FIN'), ...LigandEntries.map(([e]) => pdb(e))].join(' '),
    fin + sources.map((s) => s.text).join(''),
    contacts,
    'converted to seq: ranges on CDK2’s residues, a row a contact',
    `contacts within ${ContactCutoff} Å between heavy atoms: CDK2 with cyclin A in each copy of 1FIN (the cyclin chain it touches most), each copy’s ATP with its CDK2 chain, and each ligand with the CDK2 chain it touches most (its residues by author number, UniProt P24941’s)`,
  );
  const uniprot: Record<string, string> = {};
  for (const acc of ['P24941', 'P20248']) {
    const url = `https://rest.uniprot.org/uniprotkb/${acc}.json`;
    const text = await fetchText(url, 'application/json', `${acc}.json`);
    uniprot[acc] = JSON.parse(text).sequence.value;
    const fixture = write(`${acc}.json`, `${JSON.stringify(trimEntry(JSON.parse(text)))}\n`);
    rows.push({
      file: `complex/${acc}.json`,
      entry: `UniProtKB ${acc}, REST entry`,
      url,
      recorded: day,
      licence: ['CC-BY-4.0'],
      attribution: `UniProt Consortium, UniProtKB ${acc}, CC BY 4.0`,
      trimmed: `members kept: ${EntryMembers.join(', ')}; features without ${FeatureDropped.join(', ')} and without the ${EntryDroppedTypes.join(', ')} features`,
      derived: 'none',
      bytesDownloaded: bytes(text),
      bytesFixture: bytes(fixture),
      sha256Download: sha256(text),
      sha256Fixture: sha256(fixture),
    });
    const name = `AF-${acc}-F1-confidence_v6.json`;
    const afUrl = `https://alphafold.ebi.ac.uk/files/${name}`;
    const af = await fetchText(afUrl, 'application/json', name);
    const kept = write(name, af.endsWith('\n') ? af : `${af}\n`);
    rows.push({
      file: `complex/${name}`,
      entry: `AlphaFold DB AF-${acc}-F1, model version 6, per-residue confidence (pLDDT)`,
      url: afUrl,
      recorded: day,
      licence: ['CC-BY-4.0'],
      attribution: `AlphaFold DB (Jumper et al. 2021; Varadi et al. 2024), AF-${acc}-F1 v6, CC BY 4.0`,
      trimmed: 'none (as published)',
      derived: 'none',
      bytesDownloaded: bytes(af),
      bytesFixture: bytes(kept),
      sha256Download: sha256(af),
      sha256Fixture: sha256(kept),
    });
  }
  // what the views hold in their specs (a pairs node is inline rows only), compared with the files by the fixtures test
  const derived = {
    letters: { cdk2: t.cdk2, construct: t.cyclin, P24941: uniprot.P24941, P20248: uniprot.P20248 },
    cyclinFirst: t.cyclinFirst,
    chains: t.chains,
    contacts: t.contacts.map((c) => [c.copy, c.cdk2, c.cyclin]),
    ligands: [...t.atp, ...ligands].map((l) => ({
      comp: l.comp,
      entry: l.entry,
      chain: l.chain,
      resolution: l.resolution,
      atoms: l.atoms,
      contacts: l.contacts.map((c) => [c.residue, c.atom]),
    })),
  };
  writeFileSync(
    join(Folder, 'src', 'complex', 'derived-data.ts'),
    `${Header}\nimport type { Derived } from './derived.js';\n\n// Written by \`scripts/record.ts complex\` from PDB 1FIN and the ligands' entries, compared with \`fixtures/complex/contacts.seq\`\n// by the fixtures test. Do not edit.\n\nexport const DerivedData: Derived = ${JSON.stringify(derived)};\n`,
  );
  return rows;
}

const Header = `/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */
`;

async function main(name: string | undefined) {
  if (name !== 'p53' && name !== 'complex') throw new Error('usage: record.ts p53|complex');
  const rows = name === 'p53' ? await p53() : await complex();
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
