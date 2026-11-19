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
// (`globins` imports `@molstar/core` and `@molstar/seqx-core` built: `npx tsc -b packages/seq/core` first.) The
// datasets: `p53`, UniProt P04637's REST entry (CC-BY-4.0) and the Proteins API's variants of it, those of UniProt's
// own and ClinVar's sources; `globins`, the Pfam PF00042 seed (CC0-1.0) cut to 60 rows, UniProt's organisms of its
// rows (CC-BY-4.0) and PDB 4HHB (CC0-1.0); `plasmid`, INSDC L09137.2 (pUC19c) as ENA serves it in FASTA, converted to
// the `seq` format, its letters also as a module of the examples (a cyclic sequence read from a source has no length
// until it is opened, so a view that wraps the origin takes them inline), and the rotated variant's features and sites,
// derived from it, as GFF3; `trna`, PDB 1EHZ's mmCIF from PDBe (CC0-1.0: chain A's sequence, modified residues and
// base pairs as `seq` files) and Rfam RF00005's seed (CC0-1.0, Stockholm as published); `complex`, PDB 1FIN (CC0-1.0)
// converted to `seq` with its derived contacts, the ligands of six CC0 entries of CDK2 (AlphaFill's transplants are not
// taken: its usage policy states no CC BY 4.0 or looser licence, the data's decisions, 1), UniProt P24941 and P20248
// trimmed, and AlphaFold DB's per-residue confidence of both chains as published (CC-BY-4.0).

import { createHash } from 'node:crypto';
import { mkdirSync, readFileSync, writeFileSync } from 'node:fs';
import { gunzipSync } from 'node:zlib';
import { dirname, join } from 'node:path';
import { fileURLToPath } from 'node:url';
import { writeDerived } from '../lib/trna/generated.js';
import { moleculeSeq, pairsSeq, trimEntry as trimCif } from '../src/trna/record-trim.ts';
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
import {
  ContactCutoff,
  contactsSeq,
  type Ligand,
  moleculeSeq as complexSeq,
  trim1fin,
  trimLigand,
} from '../src/complex/record-trim.ts';
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
  const molecule = write('1FIN.seq', complexSeq(t));
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
    // as published, byte for byte (no newline added)
    const kept = write(name, af);
    rows.push({
      file: `complex/${name}`,
      entry: `AlphaFold DB AF-${acc}-F1, model version 6, per-residue confidence (pLDDT)`,
      url: afUrl,
      recorded: day,
      licence: ['CC-BY-4.0'],
      attribution: `AlphaFold DB (Jumper et al. 2021; Varadi et al. 2024), AF-${acc}-F1 v6, CC BY 4.0`,
      trimmed: 'none: as published, byte for byte',
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

const Datasets: Readonly<Record<string, () => Promise<Row[]>>> = { p53, globins, plasmid, trna, complex };

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
