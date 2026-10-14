/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The `Seq*/File formats` stories' documents and loaders (seqspec:sources#formats-and-services), one builder per story,
// no React and no Mol* (the smoke test imports this; the stories register the Mol* adapter themselves). Each builder
// gives the document, the loaders that answer its sources, the refs of the nodes that load the data (printed as the
// document has them) and the data as given: a file's text, or the requests the loaders answered (`Fetched`, filled as
// they are asked). Recorded data are the committed fixtures or text written here from them; live loaders fetch.

import { createSeqViewSpec } from '@molstar/seq-core/builder';
import type { DocNode, Loaders, SeqViewSpecDoc, StreamProvider } from '@molstar/seq-core/format';
import { archiveLoaders, readDocument, readSvsx, toJSON, toMsgPack } from '@molstar/seq-core/format';
import { bigBedProvider, bigWigProvider, twoBitProvider, vcfProvider } from '@molstar/seq-core/stream/gmod';
import { tilesStream } from '@molstar/seq-core/stream/tiles';
import { tileServer } from '@molstar/seq-core/stream/tiles/server';
import SEED_ALN from '../../../../core/_test/fixtures/pf00042-seed-8.aln?raw';
import { BRS, DESIGN_7, GLOBINS, INSERT, P53, VH, VH_KABAT } from '../../../../core/_test/format/load-fixtures.js';
import { absoluteDoc, SVSX_FILES } from '../../../../core/_test/format/svsx-fixture.js';
import { DEMO_SYSTEMS, demoGenome, memoryReader } from '../../../../core/_test/stream/gmod/demo-files.js';

// ---- Common ----

/** A request a loader answered: its URI, its size, the start of what came back. */
export interface Fetched {
  readonly uri: string;
  readonly bytes: number;
  readonly head: string;
}
/** The data as given: a file's text (or its start), a request, a list of reads. */
export interface Given {
  readonly label: string;
  readonly text: string;
}
/** A story's document and what it reads. */
export interface FormatCase {
  readonly doc: SeqViewSpecDoc;
  readonly loaders: Loaders & { readonly fetched: Fetched[] };
  /** The refs of the nodes that load the data, in the order to print them. */
  readonly shown: readonly string[];
  readonly given: readonly Given[];
}

const enc = new TextEncoder();
const dec = new TextDecoder();
/** The start of `bytes` as a reader would see it: text up to a few lines, else hex. */
export function headOf(bytes: Uint8Array, lines = 8): string {
  const sample = bytes.subarray(0, 4096);
  const printable = sample.every((b) => b === 9 || b === 10 || b === 13 || (b >= 32 && b < 127) || b >= 128);
  if (!printable) {
    const hex = Array.from(bytes.subarray(0, 48), (b) => b.toString(16).padStart(2, '0'));
    const rows: string[] = [];
    for (let k = 0; k < hex.length; k += 16) rows.push(hex.slice(k, k + 16).join(' '));
    return `${rows.join('\n')}\n… (binary, ${bytes.length.toLocaleString('en-US')} bytes)`;
  }
  return excerpt(dec.decode(sample), lines, bytes.length);
}
/** The first `lines` lines of `text` (each cut at 120 characters), and what was left out. */
export function excerpt(text: string, lines = 8, size = text.length): string {
  const all = text.split('\n');
  const kept = all.slice(0, lines).map((l) => (l.length > 120 ? `${l.slice(0, 117)}…` : l));
  const more = all.length - kept.length;
  return more > 0 || size > text.length
    ? `${kept.join('\n')}\n… (${more > 0 ? `${more} more lines, ` : ''}${size.toLocaleString('en-US')} bytes)`
    : kept.join('\n');
}

/** Loaders over recorded files by URI or host file; every request answered is logged. Anything else fails. */
export function recordedLoaders(files: Readonly<Record<string, string | Uint8Array>>): FormatCase['loaders'] {
  const fetched: Fetched[] = [];
  return {
    fetched,
    async load(request) {
      const name = request.uri ?? request.hostFile ?? '';
      const f = files[name];
      if (f === undefined) throw new Error(`no recorded file for ${name}`);
      const bytes = typeof f === 'string' ? enc.encode(f) : f;
      fetched.push({ uri: name, bytes: bytes.length, head: headOf(bytes) });
      return bytes;
    },
  };
}
/** Loaders that fetch each URI from the network (the live mode of the service stories); logged as the recorded ones. */
export function liveLoaders(): FormatCase['loaders'] {
  const fetched: Fetched[] = [];
  return {
    fetched,
    async load(request, signal) {
      const uri = request.uri ?? '';
      const r = await fetch(uri, { signal });
      if (!r.ok) throw new Error(`${uri}: HTTP ${r.status} ${r.statusText}`);
      const bytes = new Uint8Array(await r.arrayBuffer());
      fetched.push({ uri, bytes: bytes.length, head: headOf(bytes) });
      return bytes;
    },
  };
}

/** A value as the code block prints it: short lists inline, long ones and long strings abridged. */
function printed(v: unknown, indent: string): string {
  if (typeof v === 'string') return JSON.stringify(v.length > 90 ? `${v.slice(0, 87)}…` : v);
  if (v === null || typeof v !== 'object') return JSON.stringify(v);
  const inner = `${indent}  `;
  if (Array.isArray(v)) {
    if (v.every((x) => x === null || typeof x !== 'object')) {
      const shown = v.slice(0, 12).map((x) => printed(x, inner));
      return `[${shown.join(', ')}${v.length > 12 ? `, … ${v.length - 12} more` : ''}]`;
    }
    return `[\n${v.map((x) => inner + printed(x, inner)).join(',\n')}\n${indent}]`;
  }
  const entries = Object.entries(v as Record<string, unknown>).filter(([, x]) => x !== undefined);
  if (entries.length === 0) return '{}';
  return `{\n${entries.map(([k, x]) => `${inner}${JSON.stringify(k)}: ${printed(x, inner)}`).join(',\n')}\n${indent}}`;
}
/** The nodes with these refs as the document holds them (children included), abridged for reading. */
export function nodesText(doc: SeqViewSpecDoc, refs: readonly string[]): string {
  const found = new Map<string, DocNode>();
  const walk = (n: DocNode) => {
    if (n.ref && refs.includes(n.ref)) found.set(n.ref, n);
    for (const c of n.children ?? []) walk(c);
  };
  walk(doc.root);
  return refs
    .map((r) => found.get(r))
    .filter((n): n is DocNode => n !== undefined)
    .map((n) => printed(n, ''))
    .join(',\n');
}

const fasta = (header: string, seq: string) => `>${header}\n${seq.replace(/(.{60})/g, '$1\n').replace(/\n$/, '')}\n`;

// ---- Sequences ----

export const SEQUENCES_URI = 'https://example.org/files/barnase-barstar.fasta';

/** Two FASTA records by `record` (and an empty one, noted), a FASTA in place (`value`), and residues as given. */
export function sequencesCase(): FormatCase {
  const file =
    fasta('barnase P00648 Ribonuclease, Bacillus amyloliquefaciens', BRS.barnase.sequence) +
    fasta('barstar P11540 Barstar, Bacillus amyloliquefaciens', BRS.barstar.sequence) +
    '>unfinished a header whose sequence was never written\n';
  const inline = fasta('pET-insert a DNA construct', INSERT);
  const b = createSeqViewSpec({ title: 'Sequences' });
  const src = b.data.source({ ref: 'fasta_file', uri: SEQUENCES_URI, format: 'fasta' });
  const barnase = b.data.sequence({ ref: 'barnase', source_ref: src, record: 'barnase' });
  const barstar = b.data.sequence({ ref: 'barstar', source_ref: src, record: 'barstar' });
  const inSrc = b.data.source({ ref: 'fasta_value', value: inline, format: 'fasta' });
  const ins = b.data.sequence({ ref: 'insert', source_ref: inSrc });
  const raw = b.data.sequence({ ref: 'design', value: DESIGN_7 });
  const v = b.view({ ref: 'sequences', label: 'Sequences' });
  for (const s of [barnase, barstar, ins, raw]) v.block({ sequence: s });
  v.track({ label: '' }).ruler();
  v.track({ label: 'Residues' }).letters();
  return {
    doc: b.getDocument(),
    loaders: recordedLoaders({ [SEQUENCES_URI]: file }),
    shown: ['fasta_file', 'barnase', 'barstar', 'fasta_value', 'insert', 'design'],
    given: [
      { label: `${SEQUENCES_URI} (recorded)`, text: excerpt(file, 12) },
      { label: 'fasta_value: the file in place', text: inline.trimEnd() },
      { label: 'design: residues as given', text: DESIGN_7 },
    ],
  };
}

// ---- Alignments and trees ----

export const ALIGNMENT_FORMATS = ['fasta', 'a3m', 'stockholm', 'clustal'] as const;
export type AlignmentFormat = (typeof ALIGNMENT_FORMATS)[number];
const ALIGNMENT_ROWS = 8;
const ALIGNMENT_URI: Record<AlignmentFormat, string> = {
  fasta: 'https://example.org/files/PF00042-seed-8.fasta',
  a3m: 'https://example.org/files/PF00042-seed-8.a3m',
  stockholm: 'https://example.org/files/PF00042-seed-8.sto',
  clustal: 'https://example.org/files/PF00042-seed-8.aln',
};
const TREE_URI = 'https://example.org/files/PF00042-seed-8.nwk';
/** A tree over the eight rows, written for the story (not inferred): leaves by the rows' names before the `/`. */
const GLOBIN_NWK =
  '((HBA_HUMAN:0.12,(A0A096MFI1_POEFO:0.21,A0A1I7SGZ4_BURXY:0.33):0.08):0.1,((A0A010R001_9PEZI:0.25,A0A178ZLR6_9EURO:0.22):0.14,(A0A0N4SXE9_BRUPA:0.3,(A0A1Y1I094_KLENI:0.18,A0A2I0WAL8_9ASPA:0.27):0.09):0.11):0.1);';

/** The seed's first rows with their Pfam names (`ID/start-end`), the accession before the `_`. */
function globinRows() {
  return GLOBINS.slice(0, ALIGNMENT_ROWS).map((r) => {
    const residues = r.gapped.replace(/[-.]/g, '').length;
    const id = r.label;
    return { id, name: `${id}/${r.start}-${r.start + residues - 1}`, gapped: r.gapped };
  });
}
/**
 * The rows as aligned FASTA, A3M (the first row's gap columns as inserts: lower case, their gaps dropped), Stockholm or
 * Clustal (the core's fixture, written from the same rows in CLUSTAL W 2.1's layout).
 */
export function alignmentText(format: AlignmentFormat): string {
  if (format === 'clustal') return SEED_ALN;
  const rows = globinRows();
  if (format === 'fasta') return rows.map((r) => fasta(r.name, r.gapped)).join('');
  if (format === 'a3m') {
    const insert = [...rows[0].gapped].map((c) => c === '-');
    const a3m = (g: string) => [...g].map((c, i) => (!insert[i] ? c : c === '-' ? '' : c.toLowerCase())).join('');
    return rows.map((r) => fasta(r.name, a3m(r.gapped))).join('');
  }
  const width = Math.max(...rows.map((r) => r.name.length)) + 2;
  const accession = (id: string) => (id === 'HBA_HUMAN' ? 'P69905' : id.slice(0, id.indexOf('_')));
  return [
    '# STOCKHOLM 1.0',
    '#=GF ID   Globin',
    '#=GF AC   PF00042.29',
    "#=GF DE   Globin (the seed's first eight rows)",
    ...rows.map((r) => `#=GS ${r.name.padEnd(width)}AC ${accession(r.id)}`),
    '',
    ...rows.map((r) => `${r.name.padEnd(width)}${r.gapped}`),
    // Cut before its closing `//`: read, with a note.
    '',
  ].join('\n');
}

/** Pfam's PF00042 seed (its first eight rows) in one of three formats, a Newick tree over them. */
export function alignmentCase(format: AlignmentFormat): FormatCase {
  const text = alignmentText(format);
  const b = createSeqViewSpec({ title: 'PF00042, eight rows' });
  const src = b.data.source({ ref: 'alignment_file', uri: ALIGNMENT_URI[format], format });
  const aln = b.data.alignment({ ref: 'globins', source_ref: src });
  const treeSrc = b.data.source({ ref: 'tree_file', uri: TREE_URI, format: 'newick' });
  aln.tree({ ref: 'globins_tree', source_ref: treeSrc });
  const accessions =
    format === 'stockholm'
      ? b.data.annotation({ ref: 'accessions', geometry: 'rows', target: aln, from_source: aln, tag: 'GS:AC' })
      : null;
  const v = b.view({ ref: 'family', label: 'Family' }).block({ alignment: aln });
  v.track({ label: 'Members' }).rows(
    accessions ? { data_columns: [{ data: accessions, field: 'AC', label: 'Accession', type: 'text' }] } : {},
  );
  return {
    doc: b.getDocument(),
    loaders: recordedLoaders({ [ALIGNMENT_URI[format]]: text, [TREE_URI]: GLOBIN_NWK }),
    shown: ['alignment_file', 'tree_file', 'globins', ...(accessions ? ['accessions'] : [])],
    given: [
      {
        label: `${ALIGNMENT_URI[format]} (recorded)`,
        text: excerpt(text, format === 'fasta' || format === 'a3m' ? 10 : 16),
      },
      { label: `${TREE_URI} (written for the story)`, text: GLOBIN_NWK },
    ],
  };
}

// ---- Tables ----

export const TABLE_FORMATS = ['tsv', 'csv', 'json_rows', 'json_columns'] as const;
export type TableFormat = (typeof TABLE_FORMATS)[number];
/** Kabat's CDRs of a heavy chain, by Kabat numbers (H2 and H3 hold 52A and 100A–C). */
const CDRS: [number, number, string][] = [
  [31, 35, 'CDR-H1'],
  [50, 65, 'CDR-H2'],
  [95, 102, 'CDR-H3'],
];
/** Kyte–Doolittle per residue over a centered window of 7: per Kabat label, its number and insertion code apart. */
function hydropathy(): { residue: number[]; ins: string[]; kd: number[] } {
  const KD: Record<string, number> = { A: 1.8, R: -4.5, N: -3.5, D: -3.5, C: 2.5, Q: -3.5, E: -3.5, G: -0.4, H: -3.2 };
  Object.assign(KD, { I: 4.5, L: 3.8, K: -3.9, M: 1.9, F: 2.8, P: -1.6, S: -0.8, T: -0.7, W: -0.9, Y: -1.3, V: 4.2 });
  const out = { residue: [] as number[], ins: [] as string[], kd: [] as number[] };
  for (let i = 3; i < VH.length - 3; i++) {
    let s = 0;
    for (let k = i - 3; k <= i + 3; k++) s += KD[VH[k]] ?? 0;
    const label = String(VH_KABAT[i]);
    const m = /^(\d+)([A-Z]?)$/.exec(label) as RegExpExecArray;
    out.residue.push(Number(m[1]));
    out.ins.push(m[2]);
    out.kd.push(Math.round((s / 7) * 100) / 100);
  }
  return out;
}
/** Columns as a file of the format: delimited text, or JSON in either table form. */
function tableFile(format: TableFormat, columns: Record<string, readonly (string | number)[]>): string {
  const names = Object.keys(columns);
  const n = columns[names[0]].length;
  if (format === 'json_columns') return JSON.stringify(columns);
  if (format === 'json_rows')
    return `[\n${Array.from({ length: n }, (_, i) => JSON.stringify(Object.fromEntries(names.map((k) => [k, columns[k][i]])))).join(',\n')}\n]`;
  const sep = format === 'tsv' ? '\t' : ',';
  const cell = (v: string | number) => (typeof v === 'string' && v.includes(sep) ? `"${v}"` : String(v));
  return [names.join(sep), ...Array.from({ length: n }, (_, i) => names.map((k) => cell(columns[k][i])).join(sep))]
    .join('\n')
    .concat('\n');
}
const TABLE_EXT: Record<TableFormat, string> = { tsv: 'tsv', csv: 'csv', json_rows: 'json', json_columns: 'json' };

/** Trastuzumab's VH in Kabat numbering, its CDRs and a value per residue read from table files, columns renamed by `fields`. */
export function tablesCase(format: TableFormat): FormatCase {
  const kind = format.startsWith('json') ? 'json' : format;
  const cdrUri = `https://example.org/files/vh-cdrs.${TABLE_EXT[format]}`;
  const kdUri = `https://example.org/files/vh-hydropathy.${TABLE_EXT[format]}`;
  const cdrFile = tableFile(format, {
    first: CDRS.map((c) => c[0]),
    last: CDRS.map((c) => c[1]),
    name: CDRS.map((c) => c[2]),
  });
  const kdFile = tableFile(format, hydropathy());
  const b = createSeqViewSpec({ title: 'Trastuzumab VH, tables' });
  const vh = b.data.sequence({ ref: 'vh', value: VH, label: 'Trastuzumab VH' });
  vh.numbering({ name: 'kabat', method: 'labels', labels: [...VH_KABAT], default: true });
  const cdrs = b.data.annotation({
    ref: 'cdrs',
    geometry: 'ranges',
    target: vh,
    uri: cdrUri,
    format: kind,
    numbering: 'kabat',
    fields: { beg: 'first', end: 'last', label: 'name' },
  } as never);
  const kd = b.data.annotation({
    ref: 'hydropathy',
    geometry: 'values',
    target: vh,
    uri: kdUri,
    format: kind,
    numbering: 'kabat',
    fields: { pos: 'residue', ins_code: 'ins', value: 'kd' },
  } as never);
  const v = b.view({ ref: 'vh_view', label: 'VH' }).block({ sequence: vh });
  v.track({ label: '' }).ruler();
  v.track({ label: 'Residues' }).letters();
  v.track({ label: 'CDRs' }).ranges({ data: cdrs, labels: true });
  v.track({ label: 'Hydropathy', height: 2 }).values({
    data: kd,
    color: {
      theme: 'field_palette',
      params: { field: 'value', palette: { kind: 'continuous', colors: 'RdBu', value_domain: [-3, 3], reverse: true } },
    },
  } as never);
  return {
    doc: b.getDocument(),
    loaders: recordedLoaders({ [cdrUri]: cdrFile, [kdUri]: kdFile }),
    shown: ['cdrs', 'hydropathy'],
    given: [
      { label: `${cdrUri} (written for the story: Kabat's CDR definitions)`, text: excerpt(cdrFile, 6) },
      { label: `${kdUri} (computed for the story: Kyte–Doolittle, a window of 7)`, text: excerpt(kdFile, 10) },
    ],
  };
}

// ---- Structures (the Mol* adapter's formats) ----

export const STRUCTURE_FORMATS = ['bcif', 'mmcif'] as const;
export type StructureForm = (typeof STRUCTURE_FORMATS)[number];
export const AF_P53_URI = 'https://alphafold.ebi.ac.uk/files/AF-P04637-F1-model_v6.bcif';
export const INS1_URI = 'https://example.org/files/ins1.cif';
/** AlphaFold DB's bands of pLDDT (very low, low, confident, very high). */
const AFDB_BANDS = {
  kind: 'discrete',
  mode: 'absolute',
  colors: [
    ['#FF7D45', 0],
    ['#FFDB13', 50],
    ['#65CBF3', 70],
    ['#0053D6', 90],
  ],
};

/**
 * A structure file read by URI through `molstar:structure`, and a CIF category of it as a table (`from_source`,
 * `category_name`): an AlphaFold DB model's BinaryCIF with its pLDDT (`ma_qa_metric_local`), or a small mmCIF with
 * insertion codes and its B-factors from `atom_site` by author numbers (`auth_seq_id` with `pdbx_PDB_ins_code`).
 */
export function structuresCase(form: StructureForm, files: { af: Uint8Array; ins1: string }): FormatCase {
  const b = createSeqViewSpec({ title: form === 'bcif' ? 'AF-P04637-F1' : 'INS1' });
  if (form === 'bcif') {
    const src = b.data.source({ ref: 'model_file', uri: AF_P53_URI, format: 'bcif' });
    const m = b.data.complex({ ref: 'model', provider: 'molstar:structure', source_ref: src, origin: 'predicted' });
    const plddt = b.data.annotation({
      ref: 'plddt',
      geometry: 'values',
      from_source: m,
      category_name: 'ma_qa_metric_local',
      numbering: 'label',
      fields: { target: 'label_asym_id', pos: 'label_seq_id', value: 'metric_value' },
    } as never);
    const v = b.view({ ref: 'model_view', label: 'Model' }).block({ unit: 'model/A' });
    v.track({ label: '' }).ruler();
    v.track({ label: 'Residues' }).letters();
    v.track({ label: 'pLDDT (ma_qa_metric_local)', height: 2 }).values({
      data: plddt,
      encoding: 'line',
      domain: [0, 100],
      color: { theme: 'field_palette', params: { field: 'value', palette: AFDB_BANDS } },
    } as never);
    return {
      doc: b.getDocument(),
      loaders: recordedLoaders({ [AF_P53_URI]: files.af }),
      shown: ['model_file', 'model', 'plddt'],
      given: [],
    };
  }
  const src = b.data.source({ ref: 'cif_file', uri: INS1_URI, format: 'mmcif' });
  const m = b.data.complex({
    ref: 'ins1',
    provider: 'molstar:structure',
    source_ref: src,
    options: { unit_name: 'auth_asym_id' },
  } as never);
  const bf = b.data.annotation({
    ref: 'b_factors',
    geometry: 'values',
    from_source: m,
    category_name: 'atom_site',
    numbering: 'auth',
    fields: { target: 'auth_asym_id', pos: 'auth_seq_id', ins_code: 'pdbx_PDB_ins_code', value: 'B_iso_or_equiv' },
  } as never);
  const v = b.view({ ref: 'chain', label: 'Chain H' }).block({ unit: 'ins1/H' });
  v.track({ label: '' }).ruler();
  v.track({ label: 'Residues' }).letters();
  v.track({ label: 'B-factor (atom_site)' }).values({ data: bf } as never);
  const atoms = files.ins1.slice(files.ins1.indexOf('loop_\n_atom_site.'));
  return {
    doc: b.getDocument(),
    loaders: recordedLoaders({ [INS1_URI]: files.ins1 }),
    shown: ['cif_file', 'ins1', 'b_factors'],
    given: [{ label: `${INS1_URI} (synthetic, written by hand): its atom_site`, text: atoms.trim() }],
  };
}

// ---- Services ----

const UNIPROT = 'https://rest.uniprot.org/uniprotkb';
export const P53_URIS = {
  fasta: `${UNIPROT}/P04637.fasta`,
  features: `${UNIPROT}/P04637.json`,
  variants: 'https://www.ebi.ac.uk/proteins/api/variation/P04637',
  pdbe: 'https://www.ebi.ac.uk/pdbe/entry-files/download/2ocj.bcif',
  afdb: AF_P53_URI,
  pae: 'https://alphafold.ebi.ac.uk/files/AF-P04637-F1-predicted_aligned_error_v6.json',
} as const;
/** P04637's FASTA as UniProt serves it (the sequence from the committed p53.json). */
export const P53_FASTA = fasta(
  'sp|P04637|P53_HUMAN Cellular tumor antigen p53 OS=Homo sapiens OX=9606 GN=TP53 PE=1 SV=4',
  P53.sequence as string,
);

/** The nodes that load the UniProt story's data: each inline `source` is the node `<reader>.src`. */
export const UNIPROT_SHOWN = ['p53.src', 'p53', 'features.src', 'features', 'variants.src', 'variants'];
/** `uniprot`'s three resources for P04637: the sequence (with its `uniprot:P04637` numbering), features, variants. */
export function uniprotDoc(): SeqViewSpecDoc {
  const b = createSeqViewSpec({ title: 'P04637 from UniProt' });
  const p53 = b.data.sequence({ ref: 'p53', source: { service: 'uniprot', resource: 'sequence', id: 'P04637' } });
  const feats = b.data.annotation({
    ref: 'features',
    geometry: 'ranges',
    target: p53,
    source: { service: 'uniprot', resource: 'features', id: 'P04637' },
  });
  const variants = b.data.annotation({
    ref: 'variants',
    geometry: 'ranges',
    target: p53,
    source: { service: 'uniprot', resource: 'variants', id: 'P04637' },
  });
  const byCategory = { theme: 'field_palette', params: { field: 'category' } };
  const v = b.view({ ref: 'p53_view', label: 'P04637' }).block({ sequence: p53 });
  v.track({ label: '' }).ruler();
  v.track({ label: 'Residues' }).letters();
  v.track({ label: 'Domains and regions' }).ranges({
    data: feats,
    filter: { field: 'category', values: ['DNA binding', 'Region', 'Motif', 'Compositional bias'] },
    labels: true,
    color: byCategory,
  } as never);
  v.track({ label: 'Sites' }).points({
    data: feats,
    filter: { field: 'category', values: ['Modified residue', 'Cross-link', 'Binding site', 'Site'] },
    color: byCategory,
  } as never);
  v.track({ label: 'Variants' }).variants({ data: variants } as never);
  return b.getDocument();
}
/** The UniProt story's recorded responses: the FASTA, and the committed, trimmed entry and variation JSON. */
export function uniprotFiles(features: string, variants: string): Record<string, string> {
  return { [P53_URIS.fasta]: P53_FASTA, [P53_URIS.features]: features, [P53_URIS.variants]: variants };
}

/**
 * A synthetic predicted aligned error for AF-P04637-F1 (393 residues), in AlphaFold DB's JSON shape: low within the
 * DNA-binding (94–292) and tetramerization (325–356) domains, high elsewhere and between them. Not a prediction: the
 * real file is about 400 KB, so the recorded mode makes this one; the live mode reads AlphaFold DB's.
 */
export function syntheticP53Pae(): string {
  const n = (P53.sequence as string).length;
  const domain = (i: number) => (i >= 93 && i < 292 ? 1 : i >= 324 && i < 356 ? 2 : 0);
  const rows: number[][] = [];
  for (let i = 0; i < n; i++) {
    const row: number[] = [];
    for (let j = 0; j < n; j++) {
      const [a, c] = [domain(i), domain(j)];
      const base = i === j ? 0.3 : a > 0 && a === c ? 2 + Math.min(4, Math.abs(i - j) / 40) : a > 0 && c > 0 ? 24 : 28;
      const wobble = ((Math.imul(i * 31 + j, 2654435761) >>> 0) % 100) / 50;
      row.push(Math.round(Math.min(31, base + wobble)));
    }
    rows.push(row);
  }
  return JSON.stringify([{ predicted_aligned_error: rows, max_predicted_aligned_error: 31.75 }]);
}

/** The nodes that load the PDBe and AlphaFold DB story's data. */
export const PDBE_AFDB_SHOWN = ['pdbe_2ocj', '2ocj', '2ocj.sifts', 'afdb_model', 'af', 'afdb_pae', 'pae', 'coverage'];
/**
 * `pdbe` and `afdb`: PDBe's 2OCJ (its BinaryCIF carries SIFTS: the chains' `uniprot:P04637` numbering) mapped onto
 * P04637 by numbering, AlphaFold DB's model of P04637 with its pLDDT (`molstar:plddt`) and its predicted aligned error
 * (`molstar:pae`), and 2OCJ's coverage on the model through P04637 (`molstar:sifts_coverage`).
 */
export function pdbeAfdbDoc(): SeqViewSpecDoc {
  const b = createSeqViewSpec({ title: 'P04637: PDBe and AlphaFold DB' });
  const p53 = b.data.sequence({ ref: 'p53', source: { service: 'uniprot', resource: 'sequence', id: 'P04637' } });
  const pdbSrc = b.data.source({ ref: 'pdbe_2ocj', service: 'pdbe', resource: 'structure', id: '2ocj' });
  b.data.complex({ ref: '2ocj', provider: 'molstar:structure', source_ref: pdbSrc });
  b.data.mapping({
    ref: '2ocj.sifts',
    from: p53,
    to: '2ocj/entity:1',
    method: 'by_numbering',
    numbering: 'uniprot:P04637',
  } as never);
  const afSrc = b.data.source({ ref: 'afdb_model', service: 'afdb', resource: 'structure', id: 'AF-P04637-F1' });
  const af = b.data.complex({ ref: 'af', provider: 'molstar:structure', source_ref: afSrc, origin: 'predicted' });
  b.data.mapping({ ref: 'af.uniprot', from: p53, to: 'af/entity:1', method: 'identity' } as never);
  const paeSrc = b.data.source({ ref: 'afdb_pae', service: 'afdb', resource: 'pae', id: 'AF-P04637-F1' });
  const plddt = b.data.computed({ ref: 'plddt', provider: 'molstar:plddt', options: { of: af } } as never);
  const pae = b.data.computed({ ref: 'pae', provider: 'molstar:pae', options: { of: af, source: paeSrc } } as never);
  const coverage = b.data.computed({
    ref: 'coverage',
    provider: 'molstar:sifts_coverage',
    options: {
      target: 'p53',
      accession: 'P04637',
      structures: [{ complex: '2ocj', source: { service: 'pdbe', id: '2ocj' } }],
    },
  } as never);
  const v = b.view({ ref: 'model_view', label: 'AF-P04637-F1' }).block({ unit: 'af/A' });
  v.track({ label: '' }).ruler();
  v.track({ label: 'Residues' }).letters();
  v.track({ label: 'pLDDT' }).values({ data: plddt } as never);
  v.track({ label: '2OCJ chains (SIFTS)' }).ranges({ data: coverage, split: 'chain' } as never);
  v.group({ label: 'Predicted aligned error', collapsed: true })
    .track({ label: 'PAE' })
    .dot_plot({
      data: pae,
      color: {
        theme: 'field_palette',
        params: { field: 'value', palette: { kind: 'continuous', colors: 'Greens', reverse: true } },
      },
      domain: [0, 30],
    } as never);
  return b.getDocument();
}
/** The PDBe and AlphaFold DB story's recorded responses: the structures' bytes, P04637's FASTA, the synthetic PAE. */
export function pdbeAfdbFiles(bcif: { pdbe: Uint8Array; afdb: Uint8Array }): Record<string, string | Uint8Array> {
  return {
    [P53_URIS.fasta]: P53_FASTA,
    [P53_URIS.pdbe]: bcif.pdbe,
    [P53_URIS.afdb]: bcif.afdb,
    [P53_URIS.pae]: syntheticP53Pae(),
  };
}

// ---- Streamed genome files ----

/** The host's code for the genome story: a provider per file, by the sources' stream ids. */
export const GENOME_HOST_CODE = `import { bigBedProvider, bigWigProvider, fetchRanges, twoBitProvider, vcfProvider } from "@molstar/seq-core/stream/gmod";

// Each file read by byte ranges: fetchRanges(url) in a page; here memoryReader(bytes), the synthetic files in memory.
const providers = {
  genome: twoBitProvider(read("demo.2bit")),
  gc: bigWigProvider(read("demo.gc.bw")),
  genes: bigBedProvider(read("demo.genes.bb"), { category: "strand", vocabulary: ["+", "-"] }),
  variants: vcfProvider(read("demo.vcf.gz"), indexBytes /* demo.vcf.gz.tbi, read whole */, { info: ["CLNSIG"] }),
};
const loaders = { load: () => Promise.reject(new Error("no files")), stream: (src) => providers[src.stream] ?? null };`;

/** Each genome file's reader, its reads logged. */
export interface GenomeFiles {
  readonly reads: Readonly<Record<'2bit' | 'bigwig' | 'bigbed' | 'vcf', [number, number][]>>;
}

/**
 * A synthetic genome (demo-files.ts) streamed from its 2bit, BigWig, BigBed and tabix VCF over in-memory byte ranges:
 * the letters' source is a complex (a unit per sequence), the others annotations on it; chrS1 in view.
 */
export function genomeCase(): FormatCase & GenomeFiles {
  const g = demoGenome();
  const r = {
    '2bit': memoryReader(g.twobit),
    bigwig: memoryReader(g.bigwig),
    bigbed: memoryReader(g.bigbed),
    vcf: memoryReader(g.vcf),
  };
  const providers: Record<string, StreamProvider> = {
    genome: twoBitProvider(r['2bit']),
    gc: bigWigProvider(r.bigwig),
    genes: bigBedProvider(r.bigbed, { category: 'strand', vocabulary: ['+', '-'] }),
    variants: vcfProvider(r.vcf, g.tbi, { info: ['CLNSIG'] }),
  };
  const b = createSeqViewSpec({ title: 'A synthetic genome, streamed' });
  const genome = b.data.source({ ref: 'genome_2bit', stream: 'genome', data_kind: 'letters' });
  b.data.complex({ ref: 'demo', source_ref: genome, label: 'Synthetic genome' });
  b.data.source({ ref: 'gc_bigwig', stream: 'gc', data_kind: 'values', coordinate_system: 'demo' });
  const gc = b.data.annotation({ ref: 'gc', geometry: 'values', source_ref: 'gc_bigwig', label: 'GC %' });
  b.data.source({ ref: 'genes_bigbed', stream: 'genes', data_kind: 'items', coordinate_system: 'demo' });
  const genes = b.data.annotation({ ref: 'genes', geometry: 'ranges', source_ref: 'genes_bigbed' });
  b.data.source({ ref: 'variants_vcf', stream: 'variants', data_kind: 'items', coordinate_system: 'demo' });
  const variants = b.data.annotation({ ref: 'variants', geometry: 'ranges', source_ref: 'variants_vcf' });
  const v = b.view({ ref: 'genome', label: 'chrS1' }).block({ unit: 'demo/chrS1' });
  v.focus({ target: 'demo/chrS1', start: 1_500, stop: 8_000 });
  v.track({ label: '' }).ruler();
  v.track({ label: 'Letters (2bit)' }).letters({ color: b.colors.scheme('nucleotide') });
  v.track({ label: 'Genes (BigBed)', height: 3 }).ranges({
    data: genes,
    labels: true,
    color: { theme: 'field_palette', params: { field: 'category' } },
  } as never);
  v.track({ label: 'Variants (tabix VCF)' }).variants({
    data: variants,
    // By the provider's vocabulary (SNV, insertion, deletion, …): a streamed palette takes the manifest's order.
    color: { theme: 'field_palette', params: { field: 'category' } },
    glyph: {
      theme: 'field_map',
      params: {
        field: 'category',
        map: { SNV: 'circle', deletion: 'triangle', insertion: 'diamond' },
        default: 'square',
      },
    },
  } as never);
  v.track({ label: 'GC % (BigWig)', height: 2 }).values({
    data: gc,
    color: { theme: 'field_palette', params: { field: 'value', palette: { kind: 'continuous', colors: 'Viridis' } } },
  } as never);
  const fetched: Fetched[] = [];
  const sizes = DEMO_SYSTEMS.map(([n, l]) => `${n} ${l.toLocaleString('en-US')}`).join(', ');
  return {
    doc: b.getDocument(),
    loaders: {
      fetched,
      load: () => Promise.reject(new Error('no files: every source here is streamed')),
      stream: (src) => (src.stream ? (providers[src.stream] ?? null) : null),
    },
    shown: ['genome_2bit', 'demo', 'gc_bigwig', 'genes_bigbed', 'variants_vcf'],
    given: [
      {
        label: `demo.2bit, ${g.twobit.length.toLocaleString('en-US')} bytes (${sizes}), as FASTA`,
        text: excerpt(g.text.fasta, 4, g.text.fasta.length),
      },
      {
        label: `demo.gc.bw, ${g.bigwig.length.toLocaleString('en-US')} bytes, as wiggle`,
        text: excerpt(g.text.wig, 5),
      },
      {
        label: `demo.genes.bb, ${g.bigbed.length.toLocaleString('en-US')} bytes, as BED12`,
        text: excerpt(g.text.bed, 6),
      },
      {
        label: `demo.vcf.gz, ${g.vcf.length.toLocaleString('en-US')} bytes (and its .tbi, ${g.tbi.length} bytes), inflated`,
        text: excerpt(g.text.vcf, 9),
      },
    ],
    reads: { '2bit': r['2bit'].reads, bigwig: r.bigwig.reads, bigbed: r.bigbed.reads, vcf: r.vcf.reads },
  };
}

// ---- The tiles server ----

export const TILES_URI = 'https://example.org/tiles/family';
const TILES_ROWS = 4_000;
const TILES_COLUMNS = 360;

const TILES_HOST_CODE = `import { tilesStream } from "@molstar/seq-core/stream/tiles";
import { tileServer } from "@molstar/seq-core/stream/tiles/server";

// The reference server in the page: procedural rows, a column statistic served (stats/conservation/…).
const server = tileServer({ rowCount: 4000, columnCount: 360, stats: { conservation: (column) => … } });
// The client's transport per source URI; tilesStream() alone fetches under the URI over HTTP (fetchTransport).
const loaders = { stream: tilesStream((uri) => (uri === "${TILES_URI}" ? server.transport : null)) };`;

/** An alignment of 4,000 rows by `seqstar-tiles/1` from the in-process reference server, its conservation served. */
export function tilesCase(): FormatCase & { readonly calls: string[] } {
  const server = tileServer({
    rowCount: TILES_ROWS,
    columnCount: TILES_COLUMNS,
    stats: { conservation: (c) => 0.5 + 0.45 * Math.sin(c / 9) * Math.cos(c / 31) },
  });
  const tiles = tilesStream((uri) => (uri === TILES_URI ? server.transport : null));
  const b = createSeqViewSpec({ title: 'A family by tiles' });
  const src = b.data.source({ ref: 'tiles', uri: TILES_URI, protocol: 'seqstar-tiles/1' });
  const aln = b.data.alignment({
    ref: 'family',
    source_ref: src,
    row_count: TILES_ROWS,
    column_count: TILES_COLUMNS,
  });
  const cons = b.data.computed({ ref: 'cons', provider: 'column_stats', alignment: aln, measure: 'conservation' });
  const v = b.view({ ref: 'tiles_view', label: 'Family' }).block({ alignment: aln });
  v.track({ label: 'Conservation (stats tiles)' }).values({ data: cons, domain: [0, 1] } as never);
  v.track({ label: 'Rows (tiles, labels, rows)' }).rows();
  return {
    doc: b.getDocument(),
    loaders: { fetched: [], load: () => Promise.reject(new Error('no files')), stream: (s) => tiles(s) },
    shown: ['tiles', 'family', 'cons'],
    given: [
      {
        label: 'The host: the tiles client over the reference server (HTTP GET under the URI by default)',
        text: TILES_HOST_CODE,
      },
    ],
    calls: server.calls,
  };
}

// ---- Seq*'s own documents ----

export const DOCUMENT_FORMATS = ['svsj', 'svsm'] as const;
export type DocumentFormat = (typeof DOCUMENT_FORMATS)[number] | 'svsx';

/**
 * The `.svsx` fixture (`kinase.svsx`, its document naming its files by relative paths) read, and opened over its own
 * files: `archiveLoaders` carry the archive's base, so a viewer validates and loads it as `compile` does; nothing else
 * answers. The paths it read are logged as the recorded loaders log theirs.
 */
export async function archiveCase(bytes: Uint8Array) {
  const svsx = await readSvsx(bytes);
  if (!svsx.doc) throw new Error(svsx.issues.map((i) => i.message).join('; '));
  const fetched: Fetched[] = [];
  const archive = archiveLoaders(
    svsx.files,
    { load: (r) => Promise.reject(new Error(`only the archive's files (${r.uri})`)) },
    svsx.base,
  );
  const loaders: FormatCase['loaders'] = {
    ...archive,
    fetched,
    async load(request, signal) {
      const b = await archive.load(request, signal);
      fetched.push({ uri: request.uri ?? '', bytes: b.length, head: headOf(b) });
      return b;
    },
  };
  const files = [...svsx.files].map(([p, b]) => `${p.padEnd(34)} ${b.length.toLocaleString('en-US')} bytes`).join('\n');
  const c: FormatCase = {
    doc: svsx.doc,
    loaders,
    shown: ['fa', 'abl', 'sites'],
    given: [
      { label: `kinase.svsx, ${bytes.length.toLocaleString('en-US')} bytes: its files`, text: files },
      { label: `${svsx.main}, read`, text: toJSON(svsx.doc, { pretty: true }) },
    ],
  };
  return { c, svsx };
}

/** The `.svsx` fixture's document (absolute URIs) written as `.svsj` or `.svsm`, its bytes, and read back. */
export function documentCase(format: (typeof DOCUMENT_FORMATS)[number]): FormatCase & { readonly bytes: Uint8Array } {
  const original = absoluteDoc();
  const bytes = format === 'svsj' ? enc.encode(toJSON(original, { pretty: true })) : toMsgPack(original);
  const read = readDocument(bytes);
  if (!read.doc) throw new Error(read.issues.map((i) => i.message).join('; '));
  return {
    doc: read.doc,
    loaders: recordedLoaders(SVSX_FILES),
    shown: ['fa', 'abl', 'sites'],
    given: [{ label: `index.${format}, ${bytes.length.toLocaleString('en-US')} bytes`, text: headOf(bytes, 14) }],
    bytes,
  };
}
