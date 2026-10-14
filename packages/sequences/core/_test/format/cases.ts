/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The acceptance cases, written with the builder as the spec writes them (case 11 is a mapping
// table, not a document). Constants the spec leaves to fixtures are load-fixtures.ts's. Each case function sits between
// `// #region caseN` and `// #endregion`: the Seq* stories show that source as the case's builder code.

import { basePairs, createSeqViewSpec, translate } from '@molstar/seq-core/builder';
import { Q } from '@molstar/seq-core/query';
import type { SeqViewSpecDoc, SvsRegistries } from '../../src/format/index.js';
import {
  ABL1A,
  BARNASE_MATURE,
  BRS,
  CODING_CDS,
  DESIGN_7,
  GK_COLUMN,
  GLOBINS,
  HOMOLOGS_STO,
  IEP,
  ISOFORM_SEGMENTS,
  KINASE_ROWS,
  P1,
  PF00042,
  RNAFOLD_3,
  STI_NAMES,
  TRANSCRIPT,
  VH_IMGT,
  VH_KABAT,
  VL_IMGT,
  VL_KABAT,
} from './load-fixtures.js';

const P53_REGIONS = [
  { start: 0, stop: 44, label: 'Transactivation' },
  { start: 101, stop: 292, label: 'DNA binding' },
  { start: 324, stop: 356, label: 'Tetramerization' },
];

// #region case1
export function case1(): SeqViewSpecDoc {
  const b = createSeqViewSpec({ title: 'p53' });
  const p53 = b.data.sequence({ ref: 'p53', source: { service: 'uniprot', resource: 'sequence', id: 'P04637' } });
  const regions = b.data.annotation({ ref: 'regions', geometry: 'ranges', target: p53, value: P53_REGIONS });
  const zinc = b.data.annotation({
    ref: 'zinc',
    geometry: 'ranges',
    target: p53,
    label: 'Zinc ligands',
    numbering: 'uniprot:P04637',
    value: { beg: [176, 179, 238, 242] },
  });
  const hot = b.data.annotation({
    ref: 'hotspots',
    geometry: 'ranges',
    target: p53,
    numbering: 'uniprot:P04637',
    value: { beg: [175, 245, 248, 249, 273, 282] },
    provenance: { source: 'literature', method: 'curated', note: 'most frequent somatic missense sites' },
  });
  const text = b.view({ ref: 'text', label: 'Text', layout: 'text' });
  const seq = text.track({ label: 'Sequence' });
  seq.letters();
  seq.emphasis({ data: hot });
  text.track({ label: 'Regions' }).ranges({
    data: regions,
    display: 'underline',
    color: { theme: 'field_palette', params: { field: 'label', palette: { kind: 'categorical', colors: 'Set2' } } },
  });
  text.track({ label: 'Zinc' }).points({ data: zinc });
  // UniProt's features and the Proteins API's variants, as case 13 declares them (loaded when the view opens).
  const feats = b.data.annotation({
    ref: 'p53-features',
    geometry: 'ranges',
    target: p53,
    source: { service: 'uniprot', resource: 'features', id: 'P04637' },
    provenance: { source: 'UniProtKB/Swiss-Prot', method: 'curated' },
  });
  const variants = b.data.annotation({
    ref: 'p53-variants',
    geometry: 'ranges',
    target: p53,
    source: { service: 'uniprot', resource: 'variants', id: 'P04637' },
    provenance: { source: 'Proteins API', method: 'imported', note: 'polymorphisms and somatic hotspots' },
  });
  // The Tracks view, ordered like a protein browser: regions, UniProt motifs, ligands and hotspots, variants, disorder.
  const tracks = b.view({ ref: 'tracks', label: 'Tracks' });
  tracks.track({ label: '' }).ruler();
  tracks.track({ label: 'Sequence' }).letters();
  tracks.track({ label: 'Regions' }).ranges({
    data: regions,
    labels: true,
    color: { theme: 'field_palette', params: { field: 'label', palette: { kind: 'categorical', colors: 'Set2' } } },
  });
  tracks.track({ label: 'Motifs' }).ranges({
    data: feats,
    filter: { field: 'category', values: ['Motif'] },
    labels: true,
    color: { theme: 'uniform', params: { color: '#3aa58f' } },
  });
  tracks.track({ label: 'Zinc' }).points({ data: zinc });
  tracks.track({ label: 'Hotspots' }).points({ data: hot });
  tracks.track({ label: 'Variants' }).variants({ data: variants });
  tracks.track({ label: 'Disorder' }).values({
    encoding: 'line',
    annotation: {
      geometry: 'values',
      target: p53,
      uri: 'https://example.org/p53-iupred.csv',
      format: 'csv',
      numbering: 'ordinal',
      fields: { pos: 'position', value: 'iupred' },
    },
  });
  return b.getDocument();
}
// #endregion

// #region case2
export function case2(): SeqViewSpecDoc {
  const b = createSeqViewSpec({ title: 'Trastuzumab Fv' });
  const fasta = b.data.source({ ref: 'fasta', host_file: 'trastuzumab.fasta', format: 'fasta' });
  const vh = b.data.sequence({ ref: 'vh', source_ref: fasta, record: 'VH' });
  const vl = b.data.sequence({ ref: 'vl', source_ref: fasta, record: 'VL' });
  for (const [s, kabat, imgt] of [
    [vh, VH_KABAT, VH_IMGT],
    [vl, VL_KABAT, VL_IMGT],
  ] as const) {
    s.numbering({ name: 'kabat', method: 'labels', labels: [...kabat], default: true });
    s.numbering({ name: 'imgt', method: 'labels', labels: [...imgt] });
  }
  const cdrs = b.data.annotation({
    ref: 'cdrs',
    geometry: 'ranges',
    numbering: 'imgt',
    provenance: { source: 'ANARCI', method: 'numbering', references: ['doi:10.1093/bioinformatics/btv552'] },
    value: {
      target: ['vh', 'vh', 'vh', 'vl', 'vl', 'vl'],
      beg: [27, 56, 105, 27, 56, 105],
      end: [38, 65, 117, 38, 65, 117],
      label: ['CDR-H1', 'CDR-H2', 'CDR-H3', 'CDR-L1', 'CDR-L2', 'CDR-L3'],
    },
  });
  const h1 = b.data.annotation({
    ref: 'cdr-h1-kabat',
    geometry: 'ranges',
    target: vh,
    numbering: 'kabat',
    value: { beg: [31], end: ['35B'], label: ['CDR-H1 (Kabat)'] },
  });
  const fv = b.view({ ref: 'fv', label: 'Fv' }).block({ sequence: vh }).block({ sequence: vl });
  fv.track({ label: '' }).ruler();
  fv.track({ label: 'Residues' }).letters();
  fv.track({ label: 'CDRs' }).ranges({
    data: [cdrs, h1],
    color: { theme: 'field_palette', params: { field: 'label', palette: { kind: 'categorical', colors: 'Set2' } } },
  });
  return b.getDocument();
}
// #endregion

export const CLADE_A = [0, 2];
export const CLADE_B = [1, 3];
export const ARG_COLUMN = 25;

// #region case3
export function case3(): SeqViewSpecDoc {
  const b = createSeqViewSpec({ title: 'SH2' });
  const sh2 = b.data.alignment({ ref: 'sh2', source: { host_file: 'PF00017.seed.sto', format: 'stockholm' } });
  sh2.tree({ source: { host_file: 'PF00017.nwk', format: 'newick' } });
  const stat = (ref: string, rows?: number[]) =>
    b.data.computed({
      ref,
      provider: 'column_stats',
      alignment: sh2,
      measure: 'conservation',
      row_subset: rows ?? null,
    });
  const [all, a, c] = [stat('cons'), stat('cons-a', CLADE_A), stat('cons-b', CLADE_B)];
  const profile = b.data.computed({ ref: 'profile', provider: 'profile', alignment: sh2 });
  const sites = b.data.annotation({
    ref: 'sites',
    geometry: 'ranges',
    target: sh2,
    numbering: 'row',
    source: { host_file: 'sh2-sites.tsv', format: 'tsv' },
    fields: { target: 'row', beg: 'position' },
    links: [{ label: 'UniProt', url: 'https://identifiers.org/uniprot:{accession}' }],
  });
  const arg = b.data.annotation({ ref: 'flvr-arg', geometry: 'ranges', target: sh2, value: { start: [ARG_COLUMN] } });
  const meta = b.data.annotation({
    ref: 'species',
    geometry: 'rows',
    target: sh2,
    source: { host_file: 'sh2-rows.tsv', format: 'tsv' },
  });
  const v = b.view({ ref: 'family', label: 'Family' }).block({ alignment: sh2 });
  v.track({ label: 'Logo' }).logo({ data: profile });
  v.track({ label: 'Conservation' }).values({ data: all });
  v.track({ label: 'Clade A' }).values({ data: a });
  v.track({ label: 'Clade B' }).values({ data: c });
  v.track({ label: 'FLVR arginine' }).points({
    data: arg,
    glyph: { theme: 'uniform', params: { glyph: { char: '★', label: 'FLVR Arg' } } },
  });
  const members = v.track({ label: 'Members' });
  members.rows({
    color: { theme: 'conservation_shading', params: { data: all, threshold: 0.4 } },
    data_columns: [{ data: meta, field: 'species', label: 'Species', type: 'text' }],
  });
  members.points({ data: sites, filter: { field: 'category', values: ['phosphosite'] } });
  return b.getDocument();
}
// #endregion

// #region case4
export function case4(): SeqViewSpecDoc {
  const b = createSeqViewSpec({ title: 'tRNA' });
  const trna = b.data.alignment({ ref: 'trna', source: { host_file: 'RF00005.seed.sto', format: 'stockholm' } });
  const ss = b.data.annotation({
    ref: 'ss-cons',
    geometry: 'pairs',
    target: trna,
    from_source: trna,
    tag: 'GC:SS_cons',
  });
  const fold = basePairs(b.data, {
    ref: 'ss-fold',
    target: trna.row(3),
    structure: RNAFOLD_3,
    notation: 'wuss',
    provenance: { source: 'RNAfold', method: 'prediction', parameters: { temperature: 37 } },
  });
  const cons = b.data.computed({ ref: 'cons', provider: 'column_stats', alignment: trna, measure: 'conservation' });
  const v = b.view({ ref: 'rna', label: 'Alignment' }).block({ alignment: trna });
  v.track({ label: 'Consensus structure' }).pairs({ data: ss, display: 'arcs', side: 'above' });
  v.track({ label: 'Members' }).rows({ color: { theme: 'letter', params: { colors: 'PurinePyrimidine' } } });
  v.track({ label: 'Row 3, RNAfold' }).pairs({ data: fold, display: 'arcs', side: 'below' });
  v.track({ label: 'Conservation' }).values({ data: cons });
  return b.getDocument();
}
// #endregion

export const CUT_ECORI = 4;
export const CUT_XHOI = 64;

// #region case5
export function case5(): SeqViewSpecDoc {
  const b = createSeqViewSpec({ title: 'pET insert' });
  const ins = b.data.sequence({
    ref: 'insert',
    alphabet: 'dna',
    source: { host_file: 'pET-insert.fasta', format: 'fasta' },
  });
  const feats = b.data.annotation({
    ref: 'features',
    geometry: 'ranges',
    target: ins,
    source: { host_file: 'pET-insert.tsv', format: 'tsv' },
  });
  const cuts = b.data.annotation({
    ref: 'cuts',
    geometry: 'ranges',
    target: ins,
    label: 'Restriction sites',
    value: { boundary: [CUT_ECORI, CUT_XHOI], label: ['EcoRI', 'XhoI'] },
  });
  const map = b.view({ ref: 'map', label: 'Map' });
  map.track({ label: '' }).ruler();
  map.track({ label: 'Features' }).ranges({ data: feats });
  map
    .track({ label: 'Sites' })
    .ranges({ data: cuts, glyph: { theme: 'uniform', params: { glyph: { char: '✂', label: 'cut site' } } } });
  const text = b.view({ ref: 'seq', label: 'Sequence', layout: 'text' });
  text.track({ label: 'Sequence' }).letters();
  text.track({ label: 'Features' }).ranges({ data: feats, display: 'underline' });
  return b.getDocument();
}
// #endregion

const STI_CHAIN = 'C';
/** Hydrogen bonds and other contacts of the crystal STI and the four poses (the spike's 1IEP contacts). */
function contactRows() {
  const rows: Record<string, unknown[]> = { target: [], beg: [], b_target: [], b_beg: [], category: [], direction: [] };
  (IEP.poses as { contacts: [number, string, number, number][] }[]).forEach((pose, k) => {
    for (const [res, category, atom] of pose.contacts) {
      rows.target.push('1iep/A');
      rows.beg.push(res);
      rows.b_target.push(k === 0 ? `1iep/${STI_CHAIN}` : `poses/${k - 1}/ligand`);
      rows.b_beg.push(STI_NAMES[atom]);
      rows.category.push(category);
      rows.direction.push(category === 'hbond' ? (res % 2 ? 'ab' : 'ba') : null);
    }
  });
  return rows;
}

/** Case 6's data nodes (case 10 reuses them). */
function abl(b: ReturnType<typeof createSeqViewSpec>) {
  const x = b.data.complex({
    ref: '1iep',
    provider: 'molstar:structure',
    source: { service: 'pdbe', resource: 'structure', id: '1iep' },
    options: { chains: ['A', STI_CHAIN] },
  });
  const kinase = x.unit('A');
  const sti = x.unit(STI_CHAIN);
  const poses = b.data.complex({
    ref: 'poses',
    provider: 'molstar:structure',
    origin: 'docked',
    source: { host_file: '1iep_ligand_vina_poses.sdf', format: 'sdf' },
    options: { model: 'all' },
  });
  const POSE_UNITS = [0, 1, 2, 3].map((k) => poses.unit(`${k}/ligand`));
  const a1 = b.data.sequence({ ref: 'abl1a', source: { service: 'uniprot', resource: 'sequence', id: 'P00519-1' } });
  const b1 = b.data.sequence({ ref: 'abl1b', source: { service: 'uniprot', resource: 'sequence', id: 'P00519-2' } });
  b.data.mapping({ ref: 'isoforms', from: a1, to: b1, method: 'segments', segments: ISOFORM_SEGMENTS });
  b.data.mapping({ ref: 'abl1a-1iep', from: a1, to: kinase, method: 'align' });
  const contacts = b.data.annotation({
    ref: 'contacts',
    geometry: 'pairs',
    semantic_type: 'interactions:interaction',
    numbering: 'auth',
    b_numbering: 'atom',
    value: contactRows(),
  });
  const pocket = b.data.computed({
    ref: 'pocket',
    provider: 'paired',
    pairs: [contacts],
    partner: [sti, ...POSE_UNITS],
  });
  // `fingerprint()` is the interactions extension's (I1): written here as the node it adds.
  const fp = b.data.computed({
    ref: 'fp',
    provider: 'interactions:fingerprint',
    options: { pairs: [contacts], ligands: [sti, ...POSE_UNITS], system: a1 },
  });
  return { kinase, sti, a1, contacts, pocket, fp };
}
// The interactions extension's shapes (`KindStyle`, which the core cannot import; a viewer test holds the two equal).
const interactionGlyphs = () => ({
  theme: 'field_map',
  params: { map: { hbond: 'diamond', hydrophobic: 'circle', pi_stacking: 'square', cation_pi: 'star' } },
});

// #region case6
export function case6(): SeqViewSpecDoc {
  const b = createSeqViewSpec({ title: 'Imatinib in ABL' });
  const { kinase, sti, a1, contacts, pocket, fp } = abl(b);
  const chains = b
    .view({ ref: 'chains', label: 'Chains', numbering: 'auth' })
    .block({ unit: kinase })
    .block({ unit: sti });
  chains.track({ label: '' }).ruler();
  chains.track({ label: 'Residues' }).letters();
  chains
    .track({ ref: 'chains.contacts', label: 'Contacts (STI)' })
    .pairs({ data: contacts, glyph: interactionGlyphs(), arcs: 'hover' });
  // The pocket's columns named by number and letter above the fingerprint (1a's numbering is 1IEP's auth: T315).
  const pk = b.view({ ref: 'poses-view', label: 'Poses', positions: pocket }).block({ sequence: a1 });
  pk.track({ label: '' }).ruler();
  pk.track({ label: 'Pocket' }).letters();
  pk.track({ label: 'Fingerprint' }).layer({ type: 'interactions:fingerprint', data: fp, fallback: 'matrix' });
  return b.getDocument();
}
// #endregion

/** Contacts within 4 Å between barnase A and barstar D, in auth numbering (P00648 48 is 1; P11540 2 is 1). */
export function contactsAD() {
  const c = BRS.contacts as [number, number, number][];
  return {
    target: c.map(() => '1brs/A'),
    beg: c.map((r) => r[0] - 47),
    b_target: c.map(() => '1brs/D'),
    b_beg: c.map((r) => r[1] - 1),
    distance: c.map((r) => r[2]),
  };
}
export const AMINO_ACIDS = [...'ACDEFGHIKLMNPQRSTVWY'];

// #region case7
export function case7(): SeqViewSpecDoc {
  const b = createSeqViewSpec({ title: 'Barnase-barstar' });
  const s = b.data.complex({
    ref: '1brs',
    provider: 'molstar:structure',
    source: { service: 'pdbe', resource: 'structure', id: '1brs' },
  });
  const barnase = b.data.sequence({
    ref: 'barnase',
    source: { service: 'uniprot', resource: 'sequence', id: 'P00648' },
  });
  const barstar = b.data.sequence({
    ref: 'barstar',
    source: { service: 'uniprot', resource: 'sequence', id: 'P11540' },
  });
  b.data.mapping({ from: barnase, to: s.sequence('entity:1'), method: 'segments', segments: [[47, 0, 110]] });
  b.data.mapping({ from: barstar, to: s.sequence('entity:2'), method: 'align' });
  const feats = ['P00648', 'P11540'].map((acc, k) =>
    b.data.annotation({
      ref: `feat-${acc}`,
      geometry: 'ranges',
      target: [barnase, barstar][k],
      source: { service: 'uniprot', resource: 'features', id: acc },
    }),
  );
  const ss = b.data.computed({ ref: 'ss', provider: 'molstar:secondary_structure', options: { of: s } });
  const contacts = b.data.annotation({ ref: 'contacts', geometry: 'pairs', numbering: 'auth', value: contactsAD() });
  const iface = b.data.computed({ ref: 'iface', provider: 'paired', pairs: [contacts], partner: null });
  const effects = b.data.annotation({
    ref: 'effects',
    geometry: 'matrix',
    target: barnase,
    states: AMINO_ACIDS,
    source: { host_file: 'barnase-effects.json', format: 'json' },
  });
  const hyd = b.data.computed({ ref: 'hyd', provider: 'org.example:hydropathy', options: { of: barnase, window: 9 } });
  for (const copies of ['one', 'each'] as const) {
    const v = b.view({ ref: `complex-${copies}`, label: copies === 'one' ? 'Complex' : 'Each chain', copies });
    v.track({ label: '' }).ruler();
    v.track({ label: 'Residues' }).letters();
    const st = v.group({ label: 'Structure (1BRS)' });
    st.track({ label: 'Secondary structure' }).secondary({ data: ss });
    // One color each: the interface has no categories, the contacts only a distance (not a category).
    const uniform = (color: string) => ({ theme: 'uniform', params: { color } });
    st.track({ label: 'Interface' }).points({ data: iface, color: uniform('#d1495b') });
    st.track({ label: 'Contacts' }).pairs({ data: contacts, color: uniform('#4a6fa5') });
    const pr = v.group({ label: 'Predictions', collapsed: true });
    pr.track({ label: 'Variant effects' }).matrix({ data: effects });
    pr.track({ label: 'Hydropathy' }).values({ data: hyd });
    // One row per feature type, collapsed at open (a ranges layer's default).
    v.track({ label: 'UniProt' }).ranges({ data: feats, split: 'category' });
  }
  return b.getDocument();
}
// #endregion

const MODEL_ID = 'AF-BARNASE-BARSTAR';
const CHAINS = [
  ['A', 'P00648'],
  ['B', 'P11540'],
] as const;
const AFDB_BANDS = {
  kind: 'discrete',
  mode: 'absolute',
  colors: [
    ['#FF7D45', 0],
    ['#FFDB13', 50],
    ['#65CBF3', 70],
    ['#0053D6', 90],
  ],
} as const;

/**
 * Contacts within 4 Å between the model's chains A (barnase) and B (barstar), in auth numbering, which for the model is
 * UniProt's (the model holds the full chains). The residue pairs are 1BRS's, standing in for the predicted interface.
 */
export function contactsAB() {
  const c = BRS.contacts as [number, number, number][];
  return {
    target: c.map(() => 'model/A'),
    beg: c.map((r) => r[0]),
    b_target: c.map(() => 'model/B'),
    b_beg: c.map((r) => r[1]),
    distance: c.map((r) => r[2]),
  };
}

/** AlphaFold DB's PAE colors (`@molstar/seq-molstar-poc`'s `PaeScale`): dark green at 0 Å to white at 30 Å. */
const PAE_SCALE = {
  theme: 'field_palette',
  params: {
    field: 'value',
    palette: {
      kind: 'continuous',
      mode: 'absolute',
      colors: [
        ['#0e5326', 0],
        ['#3b9a54', 10],
        ['#a5d8a7', 20],
        ['#f4fbf2', 30],
      ],
      overflow_color: '#f4fbf2',
    },
  },
};

// #region case8
export function case8(): SeqViewSpecDoc {
  const b = createSeqViewSpec({ title: 'AlphaFold DB complex' });
  const m = b.data.complex({
    ref: 'model',
    provider: 'molstar:structure',
    origin: 'predicted',
    source: { service: 'afdb', resource: 'structure', id: MODEL_ID },
  });
  const feats = CHAINS.map(([chain, acc]) => {
    const q = b.data.sequence({ ref: acc, source: { service: 'uniprot', resource: 'sequence', id: acc } });
    b.data.mapping({ from: q, to: m.unit(chain), method: 'identity' });
    return b.data.annotation({
      ref: `feat-${acc}`,
      geometry: 'ranges',
      target: q,
      source: { service: 'uniprot', resource: 'features', id: acc },
    });
  });
  const plddt = b.data.annotation({
    ref: 'plddt',
    geometry: 'values',
    from_source: m,
    category_name: 'ma_qa_metric_local',
    numbering: 'label',
    fields: { target: 'label_asym_id', pos: 'label_seq_id', value: 'metric_value' },
  });
  const ss = b.data.computed({ ref: 'ss', provider: 'molstar:secondary_structure', options: { of: m } });
  // The interface: the contacts between the chains (here 1BRS's residue pairs), and the residues they touch.
  const contacts = b.data.annotation({ ref: 'contacts', geometry: 'pairs', numbering: 'auth', value: contactsAB() });
  const iface = b.data.computed({ ref: 'iface', provider: 'paired', pairs: [contacts], partner: null });
  const bands = { theme: 'field_palette', params: { field: 'value', palette: AFDB_BANDS } };
  const uniform = (color: string) => ({ theme: 'uniform', params: { color } });
  // The predicted aligned error over the complex's residues, barnase then barstar: AlphaFold DB's matrix as a grid.
  const paeSource = b.data.source({ ref: 'pae-src', service: 'afdb', resource: 'pae', id: MODEL_ID });
  const pae = b.data.computed({ ref: 'pae', provider: 'molstar:pae', options: { of: m, source: paeSource } });
  const v = b.view({ ref: 'chains', label: 'Chains' });
  v.track({ label: '' }).ruler();
  v.track({ label: 'Residues' }).letters();
  v.track({ label: 'Secondary structure' }).secondary({ data: ss });
  // The model's confidence: the scores as a line in the AlphaFold DB bands, and the bands as a thin strip under it,
  // last before a group's header (no row borders: next to another track the strip would read as that track's).
  v.track({ label: 'pLDDT', height: 2 }).values({ data: plddt, encoding: 'line', domain: [0, 100], color: bands });
  v.track({ label: '', height: 0.5 }).color({ data: plddt, color: bands });
  const itf = v.group({ label: 'Interface' });
  itf.track({ label: 'Interface residues' }).points({ data: iface, color: uniform('#d1495b') });
  itf.track({ label: 'Contacts' }).pairs({ data: contacts, color: uniform('#4a6fa5') });
  // Scored residues across, the residues aligned on down: the chains' own blocks low, the interface's patch too; its
  // contacts as dots over it. In a group closed at open: the plot is tall.
  const plot = v.group({ label: 'Aligned error (PAE)', collapsed: true }).track({ label: 'PAE, contacts' });
  plot.dot_plot({ data: pae, color: PAE_SCALE, domain: [0, 30] });
  plot.dot_plot({ data: contacts, color: uniform('#d1495b') });
  // UniProt's features by kind (the model's own helices and strands stand for UniProt's).
  const uni = v.group({ label: 'UniProt' });
  uni.track({ label: 'Processing' }).ranges({
    data: feats,
    filter: { field: 'category', values: ['Signal', 'Propeptide', 'Chain', 'Initiator methionine'] },
    labels: true,
    color: uniform('#cc9933'),
  });
  uni.track({ label: 'Active sites' }).points({
    data: feats,
    filter: { field: 'category', values: ['Active site'] },
    color: uniform('#e0427a'),
  });
  uni.track({ label: 'Mutagenesis' }).points({
    data: feats,
    filter: { field: 'category', values: ['Mutagenesis'] },
    color: uniform('#7a8ff0'),
  });
  return b.getDocument();
}
// #endregion

const [ROWS, COLS, HBA_ROW, HBB_ROW] = [PF00042.rows, PF00042.columns, PF00042.hbaRow, PF00042.hbbRow];

// #region case9
export function case9(): SeqViewSpecDoc {
  const b = createSeqViewSpec({ title: 'PF00042' });
  const aln = b.data.alignment({
    ref: 'pf00042',
    row_count: ROWS,
    column_count: COLS,
    source: { uri: 'https://example.org/pfam/PF00042/full', protocol: 'seqstar-tiles/1' },
  });
  aln.tree({ source: { uri: 'https://example.org/pfam/PF00042/full.nwk', format: 'newick' }, leaves: 'order' });
  const hb = b.data.complex({
    ref: '1a3n',
    provider: 'molstar:structure',
    source: { service: 'pdbe', resource: 'structure', id: '1a3n' },
  });
  aln.bind(HBA_ROW, hb.sequence('entity:1'), { method: 'align', ref: 'pf00042-hba' });
  aln.bind(HBB_ROW, hb.sequence('entity:2'), { method: 'align', ref: 'pf00042-hbb' });
  const cons = b.data.computed({ ref: 'cons', provider: 'column_stats', alignment: aln, measure: 'conservation' });
  const ss = b.data.computed({ ref: 'ss', provider: 'molstar:secondary_structure', options: { of: hb } });
  const meta = b.data.annotation({
    ref: 'meta',
    geometry: 'rows',
    target: aln,
    source: { uri: 'https://example.org/pfam/PF00042/rows.json', format: 'json' },
  });
  const v = b.view({ ref: 'msa', label: 'Alignment' });
  v.block({ alignment: aln });
  v.track({ label: 'Conservation' }).values({ data: cons, domain: [0, 1] });
  v.track({ label: 'Helices (1A3N α, β)' }).secondary({ data: ss });
  const rows = v.track({ label: 'Rows' });
  rows.rows({
    data_columns: [
      { data: meta, field: 'species', label: 'Species', type: 'text' },
      { data: meta, field: 'identity', label: 'Id.', type: 'fraction' },
    ],
  });
  // Each bound row's own chain's helices on its row (α's on HBA's, β's on HBB's).
  rows.secondary({ data: ss });
  v.section({ ref: 'msa.chains', label: '1A3N' })
    .block({ unit: hb.unit('A') })
    .block({ unit: hb.unit('B') })
    .track({ label: 'Secondary structure' })
    .secondary({ data: ss });
  return b.getDocument();
}
// #endregion

const HINGE_COLUMN = 50;
const KD_REGIONS = [
  ['P-loop', 248, 256],
  ['αC helix', 279, 291],
  ['Gatekeeper', 315, 316],
  ['Hinge', 316, 322],
  ['DFG', 381, 384],
  ['A-loop', 384, 403],
].map(([label, s, e]) => ({ start: Number(s) - 1, stop: Number(e) - 1, label }));

// #region case10
export function case10(): SeqViewSpecDoc {
  const b = createSeqViewSpec({ title: 'Stories drafts' });
  const { kinase, sti } = abl(b);
  const gate = b.data.annotation({
    ref: 'site.gatekeeper',
    geometry: 'ranges',
    target: kinase,
    numbering: 'auth',
    value: { beg: [315] },
  });
  const gk = b.data.alignment({
    ref: 'gatekeeper-msa',
    value: KINASE_ROWS,
    labels: ['ABL1', 'ABL1 T315I', 'SRC', 'KIT', 'PDGFRA', 'EGFR'],
    parent_rows: [null, 0, null, null, null, null],
  });
  const sites = b.data.annotation({
    ref: 'msa-sites',
    geometry: 'ranges',
    target: gk,
    value: {
      start: [GK_COLUMN, HINGE_COLUMN],
      stop: [GK_COLUMN + 1, HINGE_COLUMN + 3],
      label: ['Gatekeeper', 'Hinge'],
    },
  });
  const gkCol = b.data.annotation({
    ref: 'msa-gatekeeper',
    geometry: 'ranges',
    target: gk,
    value: { start: [GK_COLUMN] },
  });
  const ident = b.data.computed({ ref: 'msa-identity', provider: 'column_stats', alignment: gk, measure: 'identity' });

  b.data.computed({
    ref: 'sti-contacts',
    provider: 'query',
    label: 'STI contacts',
    query: Q.read('(or (partners 1iep/C :via contacts) (partners 1iep/A :via contacts))'),
  });
  const chains = b
    .view({ ref: 'chains', label: 'Chains', numbering: 'auth' })
    .block({ unit: kinase })
    .block({ unit: sti });
  chains.track({ label: 'Residues' }).letters();
  // A glyph and a color per kind of contact: where several kinds meet at STI's one cell, each shows beside the others.
  chains.track({ ref: 'chains.contacts', label: 'Contacts (STI)' }).pairs({
    data: 'contacts',
    glyph: interactionGlyphs(),
    color: { theme: 'field_palette', params: { field: 'category' } },
    arcs: 'hover',
  });
  const msa = b.view({ ref: 'alignment', label: 'Alignment', style: 'msa' }).block({ alignment: gk });
  msa.track({ label: 'Sites' }).ranges({ data: sites, display: 'brackets' });
  const kinases = msa.track({ label: 'Kinases' });
  kinases.rows({ color: { theme: 'identity_to_row', params: { row: 0 } } });
  kinases.emphasis({ data: gkCol });
  msa.track({ label: 'Identity' }).values({ data: ident });
  const kd = b.data.unit({ ref: 'abl1-kd', sequence: 'abl1a', range: [241, 493], label: 'ABL1 kinase domain' });
  const kdSites = b.data.annotation({
    ref: 'kd-sites',
    geometry: 'ranges',
    target: 'abl1a',
    numbering: 'uniprot:P00519-1',
    value: { beg: [253, 255, 315, 351] },
  });
  const text = b.view({ ref: 'text', label: 'Text', layout: 'text', style: 'compact' }).block({ unit: kd });
  const kt = text.track({ label: 'Sequence' });
  kt.letters();
  kt.emphasis({ data: kdSites });
  text
    .track({ label: 'Regions' })
    .ranges({ display: 'underline', annotation: { geometry: 'ranges', target: 'abl1a', value: KD_REGIONS } }); // six regions

  const strip = b.view({ ref: 'strip', label: 'Strip', style: 'strip', height: 'content' }).block({ unit: kinase });
  strip.track({ label: 'Residues' }).letters();
  strip.track({ label: 'Gatekeeper' }).points({ data: gate });
  const both = b.view({ ref: 'both', label: 'Annotations and MSA' });
  both
    .section({ ref: 'both.ann', height: 'content', scroll_group: 'x' })
    .block({ alignment: gk })
    .track({ label: 'Sites' })
    .ranges({ data: sites });
  both
    .section({ ref: 'both.msa', height: 'box', height_rows: 12, scroll_group: 'x' })
    .block({ alignment: gk })
    .track({ label: 'Rows' })
    .rows();

  b.snapshot({ ref: 's1', title: 'Meet ABL and imatinib', ops: [] });
  b.snapshot({
    ref: 's2',
    title: 'The pocket',
    // STI's contacts, both sides: the kinase residues touching imatinib and its atoms touching the kinase (`pocket`
    // alone leaves STI out; `contacts` keeps residues only the docked poses touch, which this track draws nothing for).
    ops: [{ op: 'set', ref: 'chains.main', params: { positions: 'sti-contacts' } }],
  });
  b.snapshot({
    ref: 's3',
    title: 'The gatekeeper',
    transition_ms: 400,
    ops: [
      { op: 'set', ref: 'root', params: { default_view: 'strip' } },
      { op: 'set', ref: 'chains.contacts', params: { hidden: true } },
      {
        op: 'insert',
        parent: 'strip.main',
        node: {
          kind: 'focus',
          ref: 'strip.focus',
          params: { target: '1iep/A', beg: 305, end: 325, numbering: 'auth' },
        },
      },
    ],
  });
  return b.getDocument();
}
// #endregion

const HEME_ELEMENTS = ['FE', 'N', 'N', 'N', 'N', ...Array.from({ length: 34 }, () => 'C'), 'O', 'O', 'O', 'O'];

// #region case12
export function case12(): SeqViewSpecDoc {
  const b = createSeqViewSpec({ title: 'Design 7' });
  const design = b.data.sequence({ ref: 'design', label: 'Design 7', value: DESIGN_7 });
  const heme = b.data.ligand({ ref: 'heme', value: HEME_ELEMENTS, chem_comp_ids: ['HEM'] });
  const dimer = b.data.complex({ ref: 'dimer', label: 'Designed dimer' });
  for (const c of ['A', 'B']) dimer.unit({ ref: `dimer.${c}`, sequence: design, label: c, origin: 'designed' });
  dimer.unit({ ref: 'dimer.heme', ligand: heme });
  const hom = b.data.alignment({ ref: 'homologs', source: { value: HOMOLOGS_STO, format: 'stockholm' } });
  hom.bind(0, design, { method: 'align' });
  const cons = b.data.computed({ ref: 'cons', provider: 'column_stats', alignment: hom, measure: 'conservation' });
  const v = b.view({ ref: 'design-view', label: 'Design' }).block({ complex: dimer });
  v.track({ label: 'Residues' }).letters();
  v.track({ label: 'Conservation' }).values({ data: cons });
  return b.getDocument();
}
// #endregion

// #region case13
export function case13(): SeqViewSpecDoc {
  const b = createSeqViewSpec({ title: 'p53 features' });
  const p53 = b.data.sequence({ ref: 'p53', source: { service: 'uniprot', resource: 'sequence', id: 'P04637' } });
  const feats = b.data.annotation({
    ref: 'p53-features',
    geometry: 'ranges',
    target: p53,
    source: { service: 'uniprot', resource: 'features', id: 'P04637' },
    provenance: { source: 'UniProtKB/Swiss-Prot', method: 'curated' },
  });
  const variants = b.data.annotation({
    ref: 'p53-variants',
    geometry: 'ranges',
    target: p53,
    source: { service: 'uniprot', resource: 'variants', id: 'P04637' },
    provenance: { source: 'Proteins API', method: 'imported', note: 'polymorphisms and somatic hotspots' },
  });
  const v = b.view({ ref: 'features', label: 'Features' });
  v.track({ label: '' }).ruler();
  v.track({ label: 'Sequence' }).letters();
  // UniProt's categories as groups, a feature type per track in UniProt's colors; some open collapsed.
  const groups = [
    ['Molecule processing', [['Chain', 'ranges', '#cc9933']], false],
    [
      'PTM',
      [
        ['Modified residue', 'points', '#c33b3b'],
        ['Cross-link', 'points', '#8a4fbf'],
      ],
      true,
    ],
    [
      'Sequence information',
      [
        ['Compositional bias', 'ranges', '#9e9e9e'],
        ['Sequence conflict', 'points', '#5b5b5b'],
      ],
      true,
    ],
    [
      'Domains & sites',
      [
        ['Domain', 'ranges', '#7a8ff0'],
        ['Region', 'ranges', '#d17b3a'],
        ['Motif', 'ranges', '#3aa58f'],
        ['Binding site', 'points', '#e0427a'],
        ['Site', 'points', '#b39020'],
      ],
      false,
    ],
  ] as const;
  // Domains, regions and motifs carry names worth reading in place (`labels`, where they fit); chains and biases not.
  const named = new Set(['Domain', 'Region', 'Motif']);
  for (const [label, tracks, collapsed] of groups) {
    const g = v.group({ label, collapsed });
    for (const [type, layer, color] of tracks)
      g.track({ label: type })[layer]({
        data: feats,
        filter: { field: 'category', values: [type] },
        color: { theme: 'uniform', params: { color } },
        ...(named.has(type) ? { labels: true } : {}),
      });
  }
  // The variants as UniProt's grid: a row per residue they change to (`split: "alt"`), expanded at open.
  v.group({ label: 'Variants' }).track({ label: 'Natural variants' }).variants({ data: variants, split: 'alt' });
  return b.getDocument();
}
// #endregion

// #region case14
export function case14(): SeqViewSpecDoc {
  const b = createSeqViewSpec({ title: 'A protein over its coding DNA' });
  const GENE =
    'GAATTCGCCACC' + // 5′ UTR
    'ATGGCTTCCAAGGAACTGACCCCACGTGACGGCC' + // exon 1's coding part
    'GTAAGTCTTTCCTAACTTCTGTTTCCACAG' + // the intron
    'ATTGGTACTTCAACAGCCAGGTGGAGTAA' + // exon 2's coding part and stop
    'GCGGCCGCATAAAGT'; // 3′ UTR
  const gene = b.data.sequence({ ref: 'gene', label: 'Gene fragment', alphabet: 'dna', value: GENE });
  const protein = translate(b.data, {
    ref: 'protein',
    label: 'Protein',
    from: gene,
    cds: [
      [12, 46],
      [76, 105],
    ],
  });
  const exons = b.data.annotation({
    ref: 'exons',
    geometry: 'ranges',
    target: gene,
    value: { start: [0, 76], stop: [46, 120], label: ['Exon 1', 'Exon 2'] },
  });
  const domain = b.data.annotation({
    ref: 'domain',
    geometry: 'ranges',
    target: protein,
    numbering: 'ordinal',
    value: { beg: [5], end: [14], label: ['Domain'] },
  });
  const site = b.data.annotation({
    ref: 'site',
    geometry: 'ranges',
    target: protein,
    numbering: 'ordinal',
    value: { beg: [17], label: ['pS17'] },
  });
  const v = b.view({ ref: 'gene-view', label: 'Gene' }).block({ sequence: gene });
  v.track({ label: '' }).ruler();
  v.track({ label: 'DNA' }).letters();
  v.track({ label: 'Protein' }).letters({ of: protein }); // M, or Met where a codon has room
  v.track({ label: 'Exons' }).ranges({ data: exons });
  v.track({ label: 'Domain' }).ranges({ data: domain });
  v.track({ label: 'Phosphosite' }).points({ data: site });
  return b.getDocument();
}
// #endregion

// Case 15's designs: synthetic (not measured, not published), each kind of edit on barnase's mature chain, written as
// edits on the wild type in its numbering (1-based). The helpers below derive the rest as data: each design's
// residues and gapped pairwise rows, an alignment of all five, and the edits as table rows on the wild type.
type IndelEdit =
  | { readonly kind: 'substitution'; readonly at: number; readonly to: string }
  | { readonly kind: 'deletion'; readonly at: number; readonly end: number }
  /** Inserted after residue `after` (0: before the first; the length: after the last). */
  | { readonly kind: 'insertion'; readonly after: number; readonly letters: string };
export interface IndelDesign {
  readonly ref: string;
  readonly label: string;
  readonly edits: readonly IndelEdit[];
}
export const INDEL_DESIGNS: readonly IndelDesign[] = [
  {
    ref: 'd1',
    label: 'D1 surface',
    edits: [
      { kind: 'substitution', at: 15, to: 'E' },
      { kind: 'substitution', at: 26, to: 'S' },
      { kind: 'substitution', at: 91, to: 'A' },
    ],
  },
  {
    ref: 'd2',
    label: 'D2 loop trim',
    edits: [
      { kind: 'deletion', at: 64, end: 67 },
      { kind: 'substitution', at: 77, to: 'D' },
    ],
  },
  {
    ref: 'd3',
    label: 'D3 loop graft',
    edits: [
      { kind: 'deletion', at: 1, end: 2 },
      { kind: 'insertion', after: 93, letters: 'GNGKDG' },
      { kind: 'substitution', at: 104, to: 'K' },
    ],
  },
  {
    ref: 'd4',
    label: 'D4 dead, tagged',
    edits: [
      { kind: 'insertion', after: 0, letters: 'M' },
      { kind: 'deletion', at: 101, end: 103 },
      { kind: 'insertion', after: 110, letters: 'GSHHHHHH' },
    ],
  },
];

/** A design applied to the wild type. */
export interface AppliedDesign {
  readonly design: IndelDesign;
  readonly residues: string;
  /** The pairwise rows, the wild type's and the design's, one column per pair. */
  readonly gapped: [string, string];
  /** Per wild-type boundary (0..length): the letters inserted there. */
  readonly inserted: ReadonlyMap<number, string>;
  /** Per wild-type index: its index in the design, -1: deleted. */
  readonly toDesign: Int32Array;
  /** Per wild-type index: the substituted letter, or "". */
  readonly substituted: readonly string[];
  /** The inserted stretches on the design's own indices: [start, stop). */
  readonly insertRuns: readonly (readonly [number, number])[];
}

export function applyDesign(design: IndelDesign, wt: string): AppliedDesign {
  const n = wt.length;
  const deleted = new Uint8Array(n);
  const substituted = Array.from({ length: n }, () => '');
  const inserted = new Map<number, string>();
  for (const e of design.edits) {
    if (e.kind === 'substitution') substituted[e.at - 1] = e.to;
    else if (e.kind === 'deletion') deleted.fill(1, e.at - 1, e.end);
    else inserted.set(e.after, (inserted.get(e.after) ?? '') + e.letters);
  }
  const toDesign = new Int32Array(n).fill(-1);
  const insertRuns: [number, number][] = [];
  let residues = '';
  let wtRow = '';
  let designRow = '';
  for (let b = 0; b <= n; b++) {
    const ins = inserted.get(b);
    if (ins) {
      insertRuns.push([residues.length, residues.length + ins.length]);
      residues += ins;
      wtRow += '-'.repeat(ins.length);
      designRow += ins;
    }
    if (b === n) break;
    wtRow += wt[b];
    if (deleted[b]) designRow += '-';
    else {
      toDesign[b] = residues.length;
      residues += substituted[b] || wt[b];
      designRow += substituted[b] || wt[b];
    }
  }
  return { design, residues, gapped: [wtRow, designRow], inserted, toDesign, substituted, insertRuns };
}

/**
 * The wild type and the designs aligned: a column per wild-type residue, and at each boundary with an insertion as
 * many insert columns as its longest one, their residues lower case (a column of only lower-case residues is an insert
 * column: seqspec:model#alignments-and-trees).
 */
function indelAlignment(applied: readonly AppliedDesign[], wt: string): { rows: string[]; labels: string[] } {
  const n = wt.length;
  const width = (b: number) => Math.max(0, ...applied.map((a) => a.inserted.get(b)?.length ?? 0));
  const row = (a: AppliedDesign | null) => {
    let r = '';
    for (let b = 0; b <= n; b++) {
      const w = width(b);
      if (w > 0) r += (a?.inserted.get(b) ?? '').toLowerCase().padEnd(w, '.');
      if (b === n) break;
      r += a === null ? wt[b] : a.toDesign[b] < 0 ? '-' : a.substituted[b] || wt[b];
    }
    return r;
  };
  return { rows: [row(null), ...applied.map(row)], labels: ['WT barnase', ...applied.map((a) => a.design.label)] };
}

/** An edit's name in HGVS's one-letter protein form (a terminal extension of a mature chain has no clean one). */
function hgvs(e: IndelEdit, wt: string): string {
  const at = (i: number) => `${wt[i - 1]}${i}`;
  if (e.kind === 'substitution') return `${at(e.at)}${e.to}`;
  if (e.kind === 'deletion') return e.at === e.end ? `${at(e.at)}del` : `${at(e.at)}_${at(e.end)}del`;
  if (e.after === 0) return `${e.letters}0`;
  if (e.after === wt.length) return `C-term +${e.letters}`;
  return `${at(e.after)}_${at(e.after + 1)}ins${e.letters}`;
}

/** A design's edits as rows on the wild type: substitutions and deletions as ranges, insertions as boundaries. */
function indelEditRows(a: AppliedDesign, wt: string) {
  const design = a.design.label;
  return a.design.edits.map((e) => {
    const label = hgvs(e, wt);
    if (e.kind === 'substitution')
      return { start: e.at - 1, stop: e.at, category: 'substitution', alt: e.to, label, design };
    if (e.kind === 'deletion') {
      const tooltip = `${label}: ${wt.slice(e.at - 1, e.end)} deleted`;
      return { start: e.at - 1, stop: e.end, category: 'deletion', alt: '-', label, tooltip, design };
    }
    const tooltip = `${label}: ${e.letters.length} inserted`;
    return { boundary: e.after, category: 'insertion', alt: e.letters, label, tooltip, design };
  });
}

/** A design's edits on its own residues: its inserted stretches as ranges, its deletions as the junctions they leave. */
function indelOwnRows(a: AppliedDesign, wt: string) {
  return a.design.edits.flatMap((e): Record<string, unknown>[] => {
    const label = hgvs(e, wt);
    if (e.kind === 'deletion') {
      const junction = e.at === 1 ? 0 : a.toDesign[e.at - 2] + 1;
      return [
        { boundary: junction, category: 'deletion', label, tooltip: `${label}: ${wt.slice(e.at - 1, e.end)} deleted` },
      ];
    }
    if (e.kind !== 'insertion') return [];
    const run = a.insertRuns.find(([s, t]) => a.residues.slice(s, t) === e.letters) ?? [0, 0];
    return [{ start: run[0], stop: run[1], category: 'insertion', label }];
  });
}

/**
 * A design's residues named in the wild type's numbering, as Kabat names an antibody's: a kept residue its wild-type
 * number, an inserted one the number before it with a letter (93A–93F), one before the first 0, -1, …
 */
function wtLabels(a: AppliedDesign): string[] {
  const out: string[] = Array.from({ length: a.residues.length }, () => '');
  a.toDesign.forEach((d, i) => {
    if (d >= 0) out[d] = String(i + 1);
  });
  for (const [s, t] of a.insertRuns) {
    const after = s === 0 ? 0 : Number(out[s - 1]);
    for (let k = s; k < t; k++) out[k] = s === 0 ? String(k - t + 1) : `${after}${String.fromCharCode(65 + k - s)}`;
  }
  return out;
}

// #region case15
export function case15(): SeqViewSpecDoc {
  const b = createSeqViewSpec({ title: "Indels: barnase designs in the wild type's numbering" });
  const wt = BARNASE_MATURE; // UniProt P00648 48–157, numbered 1–110
  const applied = INDEL_DESIGNS.map((d) => applyDesign(d, wt));
  const wtSeq = b.data.sequence({ ref: 'wt', label: 'WT barnase', alphabet: 'protein', value: wt });
  // Each design's residues also named in the wild type's numbering, inserts with letters (93A–93F): a label numbering.
  const designs = applied.map((a) =>
    b.data
      .sequence({ ref: a.design.ref, label: a.design.label, alphabet: 'protein', value: a.residues })
      .numbering({ name: 'wt', method: 'labels', labels: wtLabels(a) }),
  );
  // Each design onto the wild type: a gapped pairwise mapping (a deletion's residues map nowhere, an insertion's too).
  for (const [k, a] of applied.entries())
    b.data.mapping({ ref: `${a.design.ref}.on_wt`, from: wtSeq, to: designs[k], method: 'pairwise', gapped: a.gapped });
  // The five aligned; the inserts' columns are lower case, so flagged as insert columns.
  const aligned = indelAlignment(applied, wt);
  const aln = b.data.alignment({ ref: 'designs', value: aligned.rows, labels: aligned.labels, alphabet: 'protein' });
  aln.bind(0, wtSeq, { method: 'identity' });
  for (const [k, d] of designs.entries()) aln.bind(k + 1, d, { method: 'identity' });

  // The wild type's features: UniProt's, from the 1BRS fixture (mature numbering: UniProt's − 47).
  const feats = (BRS.barnase.features as [string, number, number, string][])
    .map(([type, s, e, note]) => ({ type, s: s - 47, e: e - 47, note }))
    .filter((f) => f.s >= 1);
  const ss = feats.filter((f) => f.type === 'Helix' || f.type === 'Beta strand');
  const ranges = (ref: string, label: string, target: typeof wtSeq, value: unknown) =>
    b.data.annotation({ ref, geometry: 'ranges', target, label, value });
  const secondary = ranges('secondary', 'Secondary structure', wtSeq, {
    start: ss.map((f) => f.s - 1),
    stop: ss.map((f) => f.e),
    category: ss.map((f) => (f.type === 'Helix' ? 'helix' : 'strand')),
  });
  const strands = ss.filter((f) => f.type === 'Beta strand');
  const sheet = ranges('sheet', 'β-sheet', wtSeq, [
    { start: strands[0].s - 1, stop: strands[strands.length - 1].e, label: 'β-sheet' },
  ]);
  const sites = ranges(
    'active_site',
    'Active site',
    wtSeq,
    feats
      .filter((f) => f.type === 'Active site')
      .map((f) => ({ start: f.s - 1, stop: f.e, label: `${wt[f.s - 1]}${f.s}`, tooltip: f.note.toLowerCase() })),
  );
  // The residues within 4 Å of barstar in 1BRS (the fixture's contacts).
  const contacts = [...new Set((BRS.contacts as [number, number, number][]).map((c) => c[0] - 47))].sort(
    (x, y) => x - y,
  );
  const iface = ranges(
    'barstar_interface',
    'Barstar interface',
    wtSeq,
    contacts.map((r) => ({ start: r - 1, stop: r, label: `${wt[r - 1]}${r}` })),
  );
  // The edits as tables on the wild type: substitutions and deletions as ranges, insertions as boundaries.
  const edits = applied.map((a) => indelEditRows(a, wt));
  const subs = applied.map((a, k) =>
    ranges(
      `${a.design.ref}.substitutions`,
      `${a.design.label}: substitutions`,
      wtSeq,
      edits[k].filter((r) => r.category === 'substitution'),
    ),
  );
  const indels = applied.map((a, k) =>
    ranges(
      `${a.design.ref}.indels`,
      `${a.design.label}: indels`,
      wtSeq,
      edits[k].filter((r) => r.category !== 'substitution'),
    ),
  );
  const allEdits = ranges('edits', 'Edits', wtSeq, edits.flat());
  // On each design's own residues: its inserted stretches, and the junctions its deletions leave.
  const own = applied.map((a, k) =>
    ranges(`${a.design.ref}.own_edits`, `${a.design.label}: edits`, designs[k], indelOwnRows(a, wt)),
  );
  const inserted = applied.flatMap((a, k) =>
    a.insertRuns.length === 0
      ? []
      : [
          ranges(
            `${a.design.ref}.inserted`,
            `${a.design.label}: inserted`,
            designs[k],
            a.insertRuns.map(([s, t]) => ({ start: s, stop: t, category: 'insertion', label: a.residues.slice(s, t) })),
          ),
        ],
  );
  // The edits' colors, with a variant for a dark skin.
  const byCategory = {
    theme: 'scheme',
    params: {
      field: 'category',
      custom: {
        kind: 'categorical',
        colors: { substitution: '#4a86d4', deletion: '#d0435b', insertion: '#16955f' },
        dark: { substitution: '#7eaaf0', deletion: '#f07c90', insertion: '#46c88d' },
      },
    },
  } as const;
  const triangle = { theme: 'uniform', params: { glyph: { shape: 'triangle', label: 'insertion' } } } as const;
  const tall = { theme: 'uniform', params: { size: 9 } } as const;
  const noCellColor = { theme: 'rules', params: { rules: [] } } as const; // the cell under the letter shows

  // 1. The wild type's axis, a row per design placed through its mapping: deleted residues leave their cells empty,
  // inserted ones have none; substitutions on a colored cell, deletions as bars, insertions as glyphs on boundaries.
  const w = b.view({ ref: 'wt_frame', label: 'On the wild type' });
  w.block({ sequence: wtSeq });
  w.track({ label: '' }).ruler();
  w.track({ label: 'WT barnase' }).letters();
  w.track({ label: 'Secondary' }).secondary({ data: secondary });
  w.track({ label: 'Active site, interface' })
    .ranges({ data: iface, display: 'union', color: { theme: 'uniform', params: { color: '#9b7fd4' } } })
    .points({ data: sites });
  const g = w.group({ label: 'Designs' });
  for (const [k, a] of applied.entries())
    g.track({ label: a.design.label })
      .ranges({ data: subs[k], color: byCategory })
      .letters({ of: designs[k], color: noCellColor })
      .ranges({ data: indels[k], color: byCategory, glyph: triangle, radius: tall });
  // The same edits as variants: a deletion one marker on a line over its residues, an insertion its length.
  w.track({ label: 'Edits as variants', height: 5 }).variants({ data: allEdits, split: 'design', color: byCategory });

  // 2. The alignment: the inserts' columns between the wild type's, which its numbering skips.
  const al = b.view({ ref: 'aligned', label: 'Aligned' });
  al.block({ alignment: aln });
  al.track({ label: '' }).ruler();
  al.track({ label: 'WT numbering' }).ruler({ of: wtSeq });
  al.track({ label: 'Secondary' }).secondary({ data: secondary });
  al.track({ label: 'Inserted' }).ranges({ data: inserted, color: byCategory });
  al.track({ label: 'Designs' }).rows();

  // 3. Each design on its own axis, named in the wild type's numbering (93A–93F for the graft), the wild type's
  // features placed through its mapping; its edits on its own residues. On the wild type, the designs' inserted
  // stretches land on the boundaries they sit at, their own data placed through their mappings.
  const axes = b.view({ ref: 'own_axes', label: 'Own axes' });
  for (const k of [-1, ...applied.keys()]) {
    const label = k < 0 ? 'WT barnase' : applied[k].design.label;
    const s = axes.section({ ref: `own_axes.${k < 0 ? 'wt' : applied[k].design.ref}`, label });
    s.block({ sequence: k < 0 ? wtSeq : designs[k] });
    s.track({ label: '' }).ruler();
    if (k >= 0) s.track({ label: 'WT numbering' }).ruler({ numbering: 'wt' });
    s.track({ label }).letters();
    s.track({ label: 'Features' }).secondary({ data: secondary }).points({ data: sites });
    s.track({ label: 'β-sheet' }).ranges({ data: sheet, labels: true });
    if (k < 0)
      s.track({ label: "Designs' inserts" }).ranges({
        data: inserted,
        color: byCategory,
        glyph: triangle,
        radius: tall,
      });
    else
      s.track({ label: 'Edits' }).ranges({
        data: own[k],
        color: byCategory,
        glyph: triangle,
        radius: tall,
        labels: true,
      });
  }
  return b.getDocument();
}
// #endregion

/** A deterministic synthetic DNA of `n` nucleotides (a linear congruential walk over ACGT): case 16's letters. */
function syntheticDna(n: number, seed: number): string {
  let x = seed;
  let out = '';
  for (let i = 0; i < n; i++) {
    x = (x * 1103515245 + 12345) & 0x7fffffff;
    out += 'ACGT'[(x >> 16) & 3];
  }
  return out;
}

// #region case16
export function case16(): SeqViewSpecDoc {
  const b = createSeqViewSpec({ title: 'pSyn, a plasmid map' });
  const p = b.data.sequence({
    ref: 'pSyn',
    label: 'pSyn',
    alphabet: 'dna',
    cyclic: true,
    value: syntheticDna(2686, 19), // random letters, a pUC19-like layout
  });
  const f: [string, number, number, string][] = [
    ['lacZα', 146, 470, '-'],
    ['MCS', 396, 454, ''],
    ['lac promoter', 514, 545, '-'],
    ['rep', 800, 1100, '+'],
    ['ori', 1155, 1744, '-'],
    ['AmpR', 1626, 2487, '-'],
    ['AmpR promoter', 2487, 2592, '-'],
    ['f1 region', 2600, 120, '+'],
  ];
  const features = b.data.annotation({
    ref: 'features',
    geometry: 'ranges',
    target: p,
    value: {
      start: f.map((x) => x[1]),
      stop: f.map((x) => x[2]),
      label: f.map((x) => x[0]),
      strand: f.map((x) => x[3]),
      category: f.map((x) => (x[0].includes('promoter') ? 'promoter' : x[0] === 'MCS' ? 'site' : 'feature')),
    },
  });
  const cut: [string, number][] = [
    ['EcoRI', 396],
    ['SacI', 402],
    ['KpnI', 408],
    ['BamHI', 417],
    ['XbaI', 423],
    ['SalI', 429],
    ['PstI', 435],
    ['SphI', 441],
    ['HindIII', 447],
    ['NdeI', 183],
    ['PvuI', 2066],
    ['ScaI', 2177],
    ['SspI', 2501],
    ['AatII', 2617],
  ];
  const sites = b.data.annotation({
    ref: 'sites',
    geometry: 'ranges',
    target: p,
    label: 'Cut sites',
    value: { boundary: cut.map((x) => x[1]), label: cut.map((x) => x[0]) },
  });
  const repeats = b.data.annotation({
    ref: 'repeats',
    geometry: 'pairs',
    target: p,
    value: {
      start: [300, 950, 2650],
      stop: [312, 962, 2662],
      b_target: ['pSyn', 'pSyn', 'pSyn'],
      b_start: [1900, 2300, 700],
      b_stop: [1912, 2312, 712],
      category: ['repeat', 'repeat', 'repeat'],
    },
  });
  const v = b.view({ ref: 'pmap', label: 'Plasmid map' });
  // A ring (labels outside with leaders, arrowheads by strand, cut sites as ticks, repeats as chords) and its linear
  // detail following it.
  b.preset.plasmidMap(v, p, { ref: 'pmap.ring', label: 'pSyn', features, sites, pairs: repeats });
  return b.getDocument();
}
// #endregion

// Case 17's samples: synthetic genotypes (2,000 samples over 300 variant sites, 0, 1 or 2 alternative alleles), the
// allele frequency varying by population and in blocks of 20 sites; seeded, so every reader sees the same samples.
export const GENOTYPE_SAMPLES = 2000;
export const GENOTYPE_SITES = 300;
export const POPULATIONS = ['AFR', 'AMR', 'EAS', 'EUR', 'SAS'] as const;

/** A small deterministic generator (mulberry32). */
function mulberry(seed: number) {
  let a = seed >>> 0;
  return () => {
    a = (a + 0x6d2b79f5) >>> 0;
    let t = a;
    t = Math.imul(t ^ (t >>> 15), t | 1);
    t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

/**
 * The samples' populations (in a shuffled data order), genotypes per site, and each sample's alternative alleles; fewer
 * samples and sites for the fuzz, which clones and diffs each document it mutates.
 */
export function genotypeSamples(samples = GENOTYPE_SAMPLES, sites = GENOTYPE_SITES) {
  const rnd = mulberry(20261013);
  const population = Array.from({ length: samples }, () => Math.floor(rnd() * POPULATIONS.length));
  // Per population, an alternative allele frequency per block of 20 sites, a site's own a little off its block's.
  const block = POPULATIONS.map(() => Array.from({ length: Math.ceil(sites / 20) }, () => rnd() ** 2.5 * 0.8));
  const freq = POPULATIONS.map((_, p) =>
    Array.from({ length: sites }, (_, i) => Math.min(0.95, block[p][Math.floor(i / 20)] * (0.6 + 0.8 * rnd()))),
  );
  const genotypes = Array.from({ length: sites }, (_, i) =>
    population.map((p) => (rnd() < freq[p][i] ? 1 : 0) + (rnd() < freq[p][i] ? 1 : 0)),
  );
  const alleles = population.map((_, q) => genotypes.reduce((n, site) => n + site[q], 0));
  return { population, genotypes, alleles };
}

// #region case17
export function case17(samples = GENOTYPE_SAMPLES, sites = GENOTYPE_SITES): SeqViewSpecDoc {
  const { population, genotypes, alleles } = genotypeSamples(samples, sites);
  const b = createSeqViewSpec({ title: 'A genotype matrix' });
  const rnd = mulberry(7);
  const region = b.data.sequence({
    ref: 'region',
    label: 'Region',
    alphabet: 'dna',
    value: Array.from({ length: sites }, () => 'ACGT'[Math.floor(rnd() * 4)]).join(''),
  });
  // The samples are the matrix's states: a value (a genotype) per state and site.
  const states = population.map((_, q) => `HG${String(q + 1).padStart(5, '0')}`);
  const index = Array.from({ length: sites }, (_, i) => i);
  const matrix = b.data.annotation({
    ref: 'genotypes',
    label: 'Genotypes',
    geometry: 'matrix',
    target: region,
    states,
    value: { index, value: genotypes },
  });
  // A key table by state: what the layer's `order` and `filter` read.
  b.data.annotation({
    ref: 'samples',
    geometry: 'rows',
    value: { state: states, population: population.map((p) => POPULATIONS[p]), alt_alleles: alleles },
  });
  const v = b.view({ ref: 'genotypes-view', label: 'Genotypes' }).block({ sequence: region });
  v.track({ label: '' }).ruler();
  // Heterozygous amber, homozygous for the alternative allele red; the reference genotype not colored.
  v.track({ ref: 'genotypes-track', label: 'Samples' }).matrix({
    ref: 'genotypes.matrix',
    data: matrix,
    color: {
      theme: 'field_palette',
      params: { palette: { kind: 'categorical', colors: { '1': '#e8a33d', '2': '#c2412d' } } },
    },
    order: 'data',
  });
  return b.getDocument();
}
// #endregion

// #region case18
export function case18(): SeqViewSpecDoc {
  const b = createSeqViewSpec({ title: 'Finding things in p53' });
  const p53 = b.data.sequence({ ref: 'p53', source: { service: 'uniprot', resource: 'sequence', id: 'P04637' } });
  // Items with names to search by: case 1's domains, the zinc ligands, and the residues the core domain's structures
  // show touching the DNA (positions of UniProt P04637).
  const domains = b.data.annotation({
    ref: 'domains',
    geometry: 'ranges',
    target: p53,
    label: 'Domains',
    value: P53_REGIONS,
  });
  const residues = (ref: string, label: string, name: string, sites: [string, number][]) =>
    b.data.annotation({
      ref,
      geometry: 'ranges',
      target: p53,
      label,
      numbering: 'uniprot:P04637',
      value: { beg: sites.map(([, at]) => at), label: sites.map(([aa, at]) => `${name} ${aa}${at}`) },
    });
  const zinc = residues('zinc', 'Zinc ligands', 'Zinc ligand', [
    ['C', 176],
    ['H', 179],
    ['C', 238],
    ['C', 242],
  ]);
  const dna = residues('dna-contacts', 'DNA-binding residues', 'DNA contact', [
    ['K', 120],
    ['S', 241],
    ['R', 248],
    ['R', 273],
    ['A', 276],
    ['C', 277],
    ['R', 280],
    ['R', 283],
  ]);
  // Two queries, run at load: p53's cysteines, and the domains holding four or more of them (a query naming the other).
  const cysteines = b.data.computed({
    ref: 'cysteines',
    provider: 'query',
    label: 'Cysteines',
    query: Q.read('(and p53 (letters C))'),
  });
  const rich = b.data.computed({
    ref: 'rich',
    provider: 'query',
    label: 'Cysteine-rich domains',
    query: Q.read('(pick (items :of domains) :where (>= (count cysteines) 4))'),
  });
  // Bookmarks: an inline ranges table under the section's name, which the bookmark actions edit and step.
  const bookmarks = b.data.annotation({
    ref: 'tracks.main.bookmarks',
    geometry: 'ranges',
    label: 'Bookmarks',
    value: [
      { target: 'p53/unit', start: 19, stop: 26, id: 'b1', label: 'MDM2 binding' },
      { target: 'p53/unit', start: 119, stop: 136, id: 'b2', label: 'L1 loop' },
      { target: 'p53/unit', start: 240, stop: 250, id: 'b3', label: 'R248' },
    ],
  });
  // The section's location sets (`highlights`): the rich domains as a band behind the tracks, the cysteines lit, the
  // zinc ligands accented, the bookmarks.
  const tracks = b.view({
    ref: 'tracks',
    label: 'Tracks',
    highlights: [
      { data: rich, style: 'band' },
      { data: cysteines, style: 'highlight' },
      { data: zinc, style: 'accent' },
      { data: bookmarks },
    ],
  });
  tracks.track({ label: '' }).ruler();
  tracks.track({ label: 'Sequence' }).letters();
  tracks.track({ label: 'Domains' }).ranges({
    data: domains,
    labels: true,
    color: { theme: 'field_palette', params: { field: 'label', palette: { kind: 'categorical', colors: 'Set2' } } },
  });
  tracks.track({ label: 'Zinc ligands' }).points({ data: zinc });
  tracks.track({ label: 'DNA-binding residues' }).points({ data: dna });
  return b.getDocument();
}
// #endregion

/** Case 19's places on chr17: a title, then [start, stop) on it (17p13.1 is hg38's cytoband; the letters are synthetic). */
export const CHR17_REGIONS: readonly (readonly [string, number, number])[] = [
  ['17p13.1', 6_500_000, 10_800_000],
  ['The whole of chr17', 0, 83_257_441],
  ['200 kb of 17p13.1', 7_600_000, 7_800_000],
  ['50 positions in an exon: the letters', 7_595_000, 7_595_050],
];

// #region case19
export function case19(): SeqViewSpecDoc {
  const b = createSeqViewSpec({ title: 'chr17 by region' });
  // GRCh38 from a letters stream: its manifest names the chromosomes and their lengths; letters come by region, when
  // zoomed in far enough to draw them. Genes, variants and a phyloP-like signal are streams over its chromosomes.
  const genome = b.data.source({ ref: 'genome', stream: 'hg38_genome', data_kind: 'letters' });
  const hg38 = b.data.complex({ ref: 'hg38', source_ref: genome, label: 'GRCh38' });
  for (const k of ['genes', 'variants']) {
    b.data.source({ ref: `${k}_src`, stream: `hg38_${k}`, data_kind: 'items', coordinate_system: 'hg38' });
    b.data.annotation({ ref: k, geometry: 'ranges', source_ref: `${k}_src` });
  }
  b.data.source({ ref: 'phylop_src', stream: 'hg38_phylop', data_kind: 'values', coordinate_system: 'hg38' });
  const phylop = b.data.annotation({ ref: 'phylop', geometry: 'values', source_ref: 'phylop_src' });
  // One chromosome, its first place in focus; a snapshot per place (the ideogram is the viewer's navigator).
  const v = b.view({ ref: 'chr17', label: 'chr17' });
  v.block({ unit: hg38.unit('chr17') });
  const [first] = CHR17_REGIONS;
  v.focus({ ref: 'chr17.focus', target: 'hg38/chr17', start: first[1], stop: first[2] });
  v.track({ label: '' }).ruler();
  v.track({ label: 'Sequence' }).letters({ color: b.colors.scheme('nucleotide') });
  v.track({ label: 'Genes', height: 4 }).ranges({
    data: 'genes',
    labels: true,
    color: { theme: 'field_palette', params: { field: 'category' } },
  });
  v.track({ label: 'Variants' }).variants({
    data: 'variants',
    color: { theme: 'letter', params: { colors: 'Nucleotide' } },
  });
  v.track({ label: 'phyloP' }).values({ data: phylop });
  for (const [k, [title, start, stop]] of CHR17_REGIONS.entries())
    b.snapshot({
      ref: `region${k + 1}`,
      title,
      transition_ms: 400,
      ops: k === 0 ? [] : [{ op: 'set', ref: 'chr17.focus', params: { start, stop } }],
    });
  return b.getDocument();
}
// #endregion

// #region case20
export function case20(): SeqViewSpecDoc {
  const b = createSeqViewSpec({ title: 'P1 on its genome' });
  const genome = b.data.source({ ref: 'genome', stream: 'hg38_genome', data_kind: 'letters' });
  const hg38 = b.data.complex({ ref: 'hg38', source_ref: genome, label: 'GRCh38' });
  for (const k of ['genes', 'variants']) {
    b.data.source({ ref: `${k}_src`, stream: `hg38_${k}`, data_kind: 'items', coordinate_system: 'hg38' });
    b.data.annotation({ ref: k, geometry: 'ranges', source_ref: `${k}_src` });
  }
  b.data.source({ ref: 'phylop_src', stream: 'hg38_phylop', data_kind: 'values', coordinate_system: 'hg38' });
  b.data.annotation({ ref: 'phylop', geometry: 'values', source_ref: 'phylop_src' });
  // P1 (synthetic, 1,334 residues) is coded by three exons over 220 kb of chr2: a codon map, its codons read across the
  // introns (residues 500 and 900 each split by one).
  const p1 = b.data.sequence({ ref: 'p1', label: 'P1', alphabet: 'protein', value: P1 });
  b.data.mapping({ ref: 'p1.codon', from: p1, to: 'hg38/chr2', method: 'codon', cds: CODING_CDS });
  // The genes stream serves T1, the transcript that codes P1, with its map (a codon map over the same exons): when a
  // region holding it arrives it becomes an item coordinate system, `genes_src/T1`, mapped onto chr2.
  const t1 = `genes_src/${TRANSCRIPT}`;
  const domains = b.data.annotation({
    ref: 'domains',
    geometry: 'ranges',
    target: p1,
    label: 'P1 domains',
    value: { start: [20, 480, 880], stop: [200, 620, 1000], label: ['Alpha', 'Beta', 'Gamma'] },
  });
  const v = b.view({ ref: 'coding', label: 'P1 on chr2' });
  // P1's own axis: the genome's genes, variants and phyloP placed through the codon map.
  const p = v.section({ ref: 'coding.protein', label: "P1, with chr2's tracks through its codon map" });
  p.block({ sequence: p1 });
  p.track({ label: '' }).ruler();
  p.track({ label: 'P1' }).letters();
  p.track({ label: 'Domains' }).ranges({ data: domains, labels: true });
  p.track({ label: 'Genes', height: 3 }).ranges({ data: 'genes', labels: true });
  p.track({ label: 'Variants' }).variants({ data: 'variants' });
  p.track({ label: 'phyloP' }).values({ data: 'phylop' });
  // chr2, in focus at the gene (5 kb either side): P1's residues over their codons and its domains along the exons, the
  // same tracks directly.
  const g = v.section({ ref: 'coding.genome', label: "chr2 at P1's gene" });
  g.block({ unit: hg38.unit('chr2') });
  g.focus({ ref: 'coding.focus', target: 'hg38/chr2', start: CODING_CDS[0][0] - 5000, stop: CODING_CDS[2][1] + 5000 });
  g.track({ label: '' }).ruler();
  g.track({ label: 'Genes', height: 3 }).ranges({ data: 'genes', labels: true });
  g.track({ label: 'P1' }).letters({ of: p1 });
  g.track({ label: 'P1 domains' }).ranges({ data: domains, labels: true });
  g.track({ label: 'Variants' }).variants({ data: 'variants' });
  g.track({ label: 'phyloP' }).values({ data: 'phylop' });
  // T1's own residues: its block names the item coordinate system, the genome's data placed through T1's map.
  const t = v.section({ ref: 'coding.transcript', label: "T1, the streamed transcript's own residues" });
  t.block({ sequence: t1 });
  t.track({ label: '' }).ruler();
  t.track({ label: 'Genes', height: 3 }).ranges({ data: 'genes', labels: true });
  t.track({ label: 'Variants' }).variants({ data: 'variants' });
  t.track({ label: 'phyloP' }).values({ data: 'phylop' });
  return b.getDocument();
}
// #endregion

/** Case 21's gene, GENE2: its protein (synthetic), translated from its CDS by the codon map. */
const GENE2_PROTEIN =
  'MSEKLNPELLKAFQEAGVDTSKWSLEELRRIAEENNLPLTGHKEDLIRRILAHFGPEAQKKSSQWEEFVAARLKEDGLDVSKLTPEQLRELAEKHGLS';
/** One codon per residue (a back-translation; the codon map translates them back). */
const CODON: Readonly<Record<string, string>> = {
  A: 'GCC',
  C: 'TGC',
  D: 'GAC',
  E: 'GAG',
  F: 'TTC',
  G: 'GGC',
  H: 'CAC',
  I: 'ATC',
  K: 'AAG',
  L: 'CTG',
  M: 'ATG',
  N: 'AAC',
  P: 'CCC',
  Q: 'CAG',
  R: 'CGC',
  S: 'AGC',
  T: 'ACC',
  V: 'GTG',
  W: 'TGG',
  Y: 'TAC',
};
const COMPLEMENT: Readonly<Record<string, string>> = { A: 'T', C: 'G', G: 'C', T: 'A' };
const reverseComplement = (s: string) =>
  [...s]
    .reverse()
    .map((c) => COMPLEMENT[c] ?? c)
    .join('');

/**
 * Case 21's region: 1,400 synthetic nucleotides with GENE2 on its reverse strand. The gene is written as its mRNA (a
 * 5′ UTR, the CDS, a 3′ UTR) cut into four exons, introns (GT…AG) between them; the region holds its reverse complement
 * between two flanks. Returned on the region's own (+) indices: the gene's span, its exons in transcript order, the
 * CDS ranges (the highest first, as a minus-strand codon map reads them), the splice junctions with read counts, and
 * an RNA-seq-like coverage.
 */
export function reverseStrandGene() {
  const cds = `${[...GENE2_PROTEIN].map((r) => CODON[r]).join('')}TAA`;
  const mrna = syntheticDna(36, 5) + cds + syntheticDna(48, 7);
  const cuts = [0, 92, 205, 301, mrna.length];
  const exonsSense = cuts.slice(1).map((e, k) => [cuts[k], e] as const);
  const introns = [150, 220, 180].map((n, k) => `GT${syntheticDna(n - 4, 11 + k)}AG`);
  const sense = exonsSense.map(([a, e], k) => mrna.slice(a, e) + (introns[k] ?? '')).join('');
  // Each exon's start on the gene's own (sense) strand.
  const at: number[] = [];
  let x = 0;
  for (const [k, [a, e]] of exonsSense.entries()) {
    at.push(x);
    x += e - a + (introns[k]?.length ?? 0);
  }
  const flank5 = syntheticDna(200, 23);
  const region = flank5 + reverseComplement(sense) + syntheticDna(1400 - 200 - sense.length, 29);
  const [geneStart, geneStop] = [flank5.length, flank5.length + sense.length];
  /** A sense-strand half-open range [a, e) of the gene as the region's own [start, stop). */
  const onRegion = (a: number, e: number): [number, number] => [geneStop - e, geneStop - a];
  const exons = exonsSense.map(([a, e], k) => onRegion(at[k], at[k] + e - a));
  // The CDS: the mRNA's [36, 36 + cds) on each exon it meets, in transcript order.
  const [cdsA, cdsE] = [36, 36 + cds.length];
  const cdsRanges = exonsSense.flatMap(([a, e], k): [number, number][] => {
    const [s, t] = [Math.max(a, cdsA), Math.min(e, cdsE)];
    return s < t ? [onRegion(at[k] + s - a, at[k] + t - a)] : [];
  });
  // Junctions: the donor (an exon's last nucleotide, read 5′ → 3′) and the acceptor (the next one's first), each on
  // the region's indices; exon 2 skipped by a few reads; two weak junctions of an antisense transcript on the + strand,
  // in the gene's introns.
  const donor = (k: number) => exons[k][0];
  const acceptor = (k: number) => exons[k][1] - 1;
  const junctions = {
    start: [donor(0), donor(1), donor(2), donor(0), exons[3][1] + 30, exons[2][1] + 40],
    b_start: [acceptor(1), acceptor(2), acceptor(3), acceptor(2), exons[2][0] - 25, exons[1][0] - 35],
    count: [48, 41, 37, 6, 3, 2],
    strand: ['-', '-', '-', '-', '+', '+'],
  };
  const depth = Array.from({ length: region.length }, (_, i) => {
    const k = exons.findIndex(([s, t]) => i >= s && i < t);
    return k < 0 ? (i % 41 === 0 ? 2 : 0) : [52, 47, 40, 36][k] + ((i * 13) % 7) - 3;
  });
  return { region, gene: [geneStart, geneStop] as [number, number], exons, cds: cdsRanges, junctions, depth };
}

// #region case21
export function case21(): SeqViewSpecDoc {
  const b = createSeqViewSpec({ title: 'GENE2, a gene on the reverse strand' });
  const g = reverseStrandGene();
  const region = b.data.sequence({ ref: 'region', label: 'Region', alphabet: 'dna', value: g.region });
  // The transcript: its exons as rows of one item (one `id`), on the minus strand, numbered in transcript order.
  const transcript = b.data.annotation({
    ref: 'gene2_201',
    geometry: 'ranges',
    target: region,
    label: 'GENE2-201',
    value: {
      start: g.exons.map(([s]) => s),
      stop: g.exons.map(([, e]) => e),
      id: g.exons.map(() => 'GENE2-201'),
      label: g.exons.map((_, k) => `exon ${k + 1}`),
      strand: g.exons.map(() => '-'),
    },
  });
  // Its protein, translated now from the CDS's ranges read downward, complemented: a minus-strand codon map.
  const protein = translate(b.data, {
    ref: 'gene2_protein',
    label: 'GENE2 protein',
    from: region,
    cds: g.cds,
    strand: '-',
  });
  const junctions = b.data.annotation({
    ref: 'junctions',
    geometry: 'pairs',
    target: region,
    value: { ...g.junctions, b_target: g.junctions.start.map(() => 'region') },
  });
  const coverage = b.data.annotation({
    ref: 'coverage',
    geometry: 'values',
    target: region,
    value: { index: g.depth.map((_, i) => i), value: g.depth },
  });
  // The gene as a piece of the region in its own direction: read 5′ → 3′, the ruler counting down, its letters
  // complemented (a letter map; the data stays the + strand's).
  const piece = { sequence: region, range: g.gene, direction: -1, label: 'GENE2 (−)' } as const;
  const v = b.view({ ref: 'gene2', label: 'GENE2' });
  const own = v.section({ ref: 'gene2.own', label: 'GENE2 on its own strand, 5′ → 3′' }).block(piece);
  own.track({ label: '' }).ruler();
  own.track({ label: 'DNA (− strand)' }).letters({ letter_map: 'complement' });
  own.track({ label: 'GENE2-201' }).ranges({ data: transcript, labels: true, direction_field: 'strand' });
  own.track({ label: 'Protein' }).letters({ of: protein });
  own.track({ label: 'Coverage', height: 2 }).values({ data: coverage });
  b.preset.sashimi(own, junctions, { label: 'Junctions' });
  // The mRNA as its own track: the same piece showing only the exons (`positions`), its letters as RNA's, complemented.
  const mrna = v.section({ ref: 'gene2.mrna', label: 'The mRNA, 5′ → 3′', follow: own.handle }).block(piece);
  b.preset.exonsOnly(mrna, transcript);
  mrna.track({ label: '' }).ruler();
  mrna.track({ label: 'mRNA' }).letters({ letter_map: 'complement_rna' });
  mrna.track({ label: 'Protein' }).letters({ of: protein });
  mrna.track({ label: 'Exons' }).ranges({ data: transcript, labels: true });
  // The region as written, the + strand: the gene's arrowheads point left, its junctions hang below the line.
  const ref = v.section({ ref: 'gene2.reference', label: 'The region as written, + strand', follow: own.handle });
  ref.block({ sequence: region });
  ref.track({ label: '' }).ruler();
  ref.track({ label: 'DNA (+ strand)' }).letters();
  ref.track({ label: 'GENE2-201' }).ranges({ data: transcript, labels: true, direction_field: 'strand' });
  ref.track({ label: 'Coverage', height: 2 }).values({ data: coverage });
  b.preset.sashimi(ref, junctions, { label: 'Junctions' });
  return b.getDocument();
}
// #endregion

/**
 * Case 22's review status: Swiss-Prot names its entries by a mnemonic (HBA_HUMAN), TrEMBL by its accession
 * (A0A010R001_9PEZI); an entry name read so, `reviewed` or `unreviewed`.
 */
const ACCESSION = /^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})_/;
export const reviewStatus = (label: string) => (ACCESSION.test(label) ? 'unreviewed' : 'reviewed');
/** Case 22's rows layer, the target of its triage patches. */
export const FAMILY_MEMBERS = 'family.members';

// #region case22
export function case22(): SeqViewSpecDoc {
  const b = createSeqViewSpec({ title: 'Globins triaged' });
  // PF00042's 32 globins (Pfam's seed as committed): a guide tree computed at load (neighbor joining on p-distances).
  const aln = b.data.alignment({
    ref: 'globins',
    label: 'PF00042',
    value: GLOBINS.map((r) => r.gapped),
    labels: GLOBINS.map((r) => r.label),
  });
  aln.tree({ provider: 'nj', options: { distance: 'p' } });
  // Key tables on the rows: the review status (the document's), identity to HBA_HUMAN and coverage (`row_stats`),
  // clusters at 50% identity (`redundancy`, taken from the most identical down).
  const review = b.data.annotation({
    ref: 'review',
    geometry: 'rows',
    target: aln,
    value: { row: GLOBINS.map((_, k) => k), status: GLOBINS.map((r) => reviewStatus(r.label)) },
  });
  const stats = b.data.computed({ ref: 'stats', provider: 'row_stats', alignment: aln });
  const nr50 = b.data.computed({
    ref: 'nr50',
    provider: 'redundancy',
    alignment: aln,
    threshold: 0.5,
    order: { data: 'stats', field: 'identity', direction: 'desc' },
  });
  const cons = b.data.computed({ ref: 'cons', provider: 'column_stats', alignment: aln, measure: 'conservation' });
  const v = b.view({ ref: 'family', label: 'PF00042', style: 'msa', height: 'box', height_rows: 16 });
  v.block({ alignment: aln });
  v.track({ label: '' }).ruler();
  v.track({ label: 'Conservation' }).values({ data: cons, domain: [0, 1] });
  // In tree order (the tree drawn, its clades folding); the triage patches the layer's `order` and `filter`.
  v.track({ label: 'Members' }).rows({
    ref: FAMILY_MEMBERS,
    order: 'tree',
    data_columns: [
      { data: review, field: 'status', label: 'Status', type: 'text' },
      { data: stats, field: 'identity', label: 'Identity', type: 'fraction' },
      { data: stats, field: 'coverage', label: 'Coverage', type: 'fraction' },
      { data: nr50, field: 'cluster', label: 'Cluster', type: 'number' },
    ],
  });
  return b.getDocument();
}
// #endregion

export const CASES: Record<string, () => SeqViewSpecDoc> = {
  case1,
  case2,
  case3,
  case4,
  case5,
  case6,
  case7,
  case8,
  case9,
  case10,
  case12,
  case13,
  case14,
  case15,
  case16,
  case17,
  case18,
  case19,
  case20,
  case21,
  case22,
};
export type { SvsRegistries };
export { ABL1A };
