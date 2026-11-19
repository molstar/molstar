/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The render tests' specs: p53 as tracks (letters, domains, conservation), p53 on a toy gene of its own (the gene's
// section following the protein's), and a 3M-position system read by region from a synthetic runtime with levels.

import { categoryColumn, type Node } from '@molstar/core/repr/tree';
import {
  type Coverage,
  type FormatRuntime,
  type OpenedSource,
  type RegionAnswer,
  type RegionRequest,
} from '@molstar/seqx-core/data';
import { createSeqSpec } from '@molstar/seqx-core/lang';

/** Human p53 (UniProt P04637), 393 residues. */
export const P53 =
  'MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGPDEAPRMPEAAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGTAKSVTCTYSPALNKMFCQLAKTCPVQLWVDSTPPPGTRVRAMAIYKQSQHMTEVVRRCPHHERCSDSDGLAPPQHLIRVEGNLRVEYLDDRNTFRHSVVVPYEPPEVGSDCTTIHYNYMCNSSCMGGMNRRPILTIITLEDSSGNLLGRNSFEVRVCACPGRDRRTEEENLRKKGEPHHELPPGSTKRALPNNTSSSPQPKKKPLDGEYFTLQIRGRERFEMFRELNEALELKDAQAGKEPGGSRAHSSHLKSKKGQSTSRHKKLMFKTEGPDSD';

/** A seeded generator (a 32-bit LCG): fixtures are generated, never stored. */
export function seeded(seed: number): () => number {
  let s = seed >>> 0;
  return () => {
    s = (Math.imul(s, 1664525) + 1013904223) >>> 0;
    return s / 0x100000000;
  };
}

/** Conservation: a value a residue in [0, 1], two missing (NaN). */
export function conservation(): Float32Array {
  const r = seeded(53);
  const v = Float32Array.from({ length: P53.length }, () => Math.round(r() * 100) / 100);
  v[100] = Number.NaN;
  v[101] = Number.NaN;
  return v;
}

function domains(b: ReturnType<typeof createSeqSpec>, p53: Parameters<typeof b.data.ranges>[1]['target']) {
  return b.data.ranges('domains', {
    target: p53,
    rows: {
      start: new Int32Array([0, 63, 93, 318, 323, 363]),
      end: new Int32Array([61, 94, 292, 360, 355, 393]),
      id: ['TAD', 'PRR', 'DBD', 'NLS', 'TET', 'CTD'],
      label: [
        'Transactivation',
        'Proline-rich',
        'DNA-binding',
        'Nuclear localization',
        'Tetramerization',
        'C-terminal',
      ],
      category: ['activation', 'linker', 'binding', 'signal', 'oligomer', 'regulatory'],
    },
  });
}

/**
 * A light chain's rulers: 107 positions numbered by Kabat, 100 of the constant domain not, under a ruler marking them
 * and one labelling them by the ordinal (`unnumbered: fallback`).
 */
export function holeySpec(): Node {
  const b = createSeqSpec({ title: 'holey rulers' });
  const labels = Array.from({ length: 207 }, (_, i) => (i < 107 ? String(i + 1) : null));
  const light = b.data.sequence('light', { alphabet: 'protein', letters: 'A'.repeat(207) }, (s) => {
    s.numbering('kabat', { method: 'labels', labels, scheme: 'kabat' });
  });
  b.views.view('main', (v) => {
    const s = v.section('light', { of: light });
    s.track('mark', { label: 'Light' }).layer('r', { type: 'ruler', numbering: 'kabat' });
    s.track('fallback', { label: 'Light, ordinal' }).layer('r', {
      type: 'ruler',
      numbering: 'kabat',
      unnumbered: 'fallback',
    });
    s.track('letters', { label: 'Sequence' }).layer('letters', { type: 'letters' });
  });
  return b.build();
}

/** p53 as tracks: its letters, its domains (in lanes), its conservation. */
export function p53TracksSpec(): Node {
  const b = createSeqSpec({ title: 'p53 tracks' });
  const p53 = b.data.sequence('p53', { alphabet: 'protein', letters: P53 });
  const d = domains(b, p53);
  const c = b.data.values('conservation', { target: p53, values: conservation() });
  b.views.view('main', (v) => {
    const s = v.section('protein', { of: p53 });
    s.track('letters', { label: 'p53' }).layer('letters', { type: 'letters' });
    s.track('domains', { label: 'Domains' }).layer('ranges', { type: 'ranges', data: d });
    s.track('conservation', { label: 'Conservation' }).layer('values', { type: 'values', data: c });
  });
  return b.build();
}

/** The toy gene's exons (minus strand, genomic order): 300 + 420 + 462 = 1182 nucleotides, 393 codons and a stop. */
export const GENE_EXONS = { start: [100, 480, 960], end: [400, 900, 1422] };

export const GENE_LENGTH = 1500;

/** The toy gene's letters: seeded. */
export function geneLetters(): string {
  const r = seeded(17);
  let s = '';
  for (let i = 0; i < GENE_LENGTH; i++) s += 'ACGT'[Math.floor(r() * 4)];
  return s;
}

/**
 * p53 on its gene: the protein's section (letters, domains), the gene's following it (letters, the domains on it), read
 * along it (a minus strand: reversed); with `coding`, on its coding pieces (its introns collapsed), the protein's
 * letters a codon wide on it too.
 */
export function p53GeneSpec(options: { readonly coding?: boolean } = {}): Node {
  const b = createSeqSpec({ title: 'p53 on its gene' });
  const p53 = b.data.sequence('p53', { alphabet: 'protein', letters: P53 });
  const gene = b.data.sequence('tp53', { alphabet: 'dna', letters: geneLetters() });
  b.data.map('p53-on-gene', { from: p53, to: gene, method: 'codon', exons: GENE_EXONS, strand: '-' });
  const d = domains(b, p53);
  b.views.view('main', (v) => {
    const protein = v.section('protein', { of: p53 });
    protein.track('letters', { label: 'p53' }).layer('letters', { type: 'letters' });
    protein.track('domains', { label: 'Domains' }).layer('ranges', { type: 'ranges', data: d });
    const g = v.section('gene', { of: gene, follow: protein, ...(options.coding ? { pieces: 'coding' } : {}) });
    if (options.coding) g.track('residues', { label: 'p53' }).layer('letters', { type: 'letters', data: p53 });
    g.track('letters', { label: 'TP53' }).layer('letters', { type: 'letters' });
    g.track('domains', { label: 'Domains on the gene' }).layer('ranges', { type: 'ranges', data: d });
  });
  return b.build();
}

// ---- a 3M-position system read by region ----

export const BIG = 3_000_000;

/** Bin widths the synthetic source serves. */
export const BigLevels = [1, 64, 4096];

const BigCoverage: Coverage = {
  systems: [{ name: 'chr', length: BIG, alphabet: 'dna' }],
  shapes: [
    { system: 'chr', shape: 'letters' },
    { system: 'chr', shape: 'ranges' },
    { system: 'chr', shape: 'values' },
  ],
  levels: BigLevels,
};

const hash = (p: number) => (Math.imul(p ^ 0x9e3779b9, 0x85ebca6b) >>> 13) & 0xffff;

/** The big system's value at a position (level 1). */
export const bigValue = (p: number) => hash(p) / 0xffff;

/**
 * Items every 50 positions, 70 long (two lanes), each `r<start>`; a region's answer holds the items overlapping it
 * (the contract), so one starting in the tile before is in it too; a level above 1 answers their density.
 */
function bigAnswer(r: RegionRequest): RegionAnswer {
  const tables: { letters?: unknown; ranges?: unknown; values?: unknown } = {};
  const L = r.level;
  if (r.shapes.includes('letters') && L === 1) {
    let s = '';
    for (let p = r.start; p < r.end; p++) s += 'ACGT'[hash(p) & 3];
    tables.letters = { form: 'letters', start: r.start, letters: s };
  }
  if (r.shapes.includes('ranges')) {
    if (L === 1) {
      const starts: number[] = [];
      for (let p = Math.max(0, Math.ceil((r.start - 69) / 50) * 50); p < r.end; p += 50) starts.push(p);
      tables.ranges = {
        form: 'list',
        start: Int32Array.from(starts),
        end: Int32Array.from(starts, (s) => Math.min(s + 70, BIG)),
        id: categoryColumn(starts.map((s) => `r${s}`)),
      };
    } else {
      const lo = Math.floor(r.start / L);
      const hi = Math.ceil(r.end / L);
      const bins = Int32Array.from({ length: hi - lo }, (_, i) => lo + i);
      tables.ranges = { form: 'density', bin: bins, count: Float32Array.from(bins, (b) => (L / 50) * (1 + (b % 3))) };
    }
  }
  if (r.shapes.includes('values')) {
    if (L === 1)
      tables.values = {
        form: 'positions',
        start: r.start,
        value: Float32Array.from({ length: r.end - r.start }, (_, i) => hash(r.start + i) / 0xffff),
      };
    else {
      const start = Math.floor(r.start / L) * L;
      const n = Math.ceil(r.end / L) - start / L;
      const v = Float32Array.from({ length: n }, (_, i) => hash(start + i * L) / 0xffff);
      tables.values = { form: 'bins', start, min: v, max: v, mean: v, count: new Float32Array(n).fill(L) };
    }
  }
  return { ok: true, tables: tables as never };
}

export const BigRuntime: FormatRuntime = {
  name: 'synthetic',
  formats: ['synthetic'],
  async open(): Promise<OpenedSource> {
    return {
      coverage: BigCoverage,
      read: async (r) => bigAnswer(r),
      dispose() {},
    };
  },
};

/**
 * A 3M-position system from the synthetic runtime: letters, ranges and values, each read by region; with `arrow`, the
 * ranges draw an arrowhead at their end.
 */
export function bigSpec(options: { arrow?: boolean } = {}): Node {
  const b = createSeqSpec({ title: '3M positions' });
  const src = b.data.source('synthetic', { format: 'synthetic', uri: 'synthetic:chr' });
  const chr = b.data.sequence('chr', { alphabet: 'dna', source: src });
  const r = b.data.ranges('items', { target: chr, source: src });
  const v = b.data.values('signal', { target: chr, source: src });
  b.views.view('main', (view) => {
    const s = view.section('chr', { of: chr });
    s.track('letters').layer('letters', { type: 'letters' });
    s.track('items').layer('ranges', { type: 'ranges', data: r, ...(options.arrow ? { arrow: 'end' } : {}) });
    s.track('signal').layer('values', { type: 'values', data: v });
  });
  return b.build();
}

/** Single-position pieces on a reversed axis (EC-33): every other position of a 300,000-nucleotide sequence. */
export const PIECES = 150_000;

/**
 * A sequence of `2 · PIECES` nucleotides shown reversed on its every other position (a ranges node of `PIECES`
 * single-position rows as its `pieces`): 150,000 pieces and as many breaks, with its letters and a signal. The pieces
 * node is not drawn as a ranges layer: a ranges layout walks every row of the entries it reads (TODO), which is not the
 * axis's cost.
 */
export function piecesSpec(): Node {
  const b = createSeqSpec({ title: 'pieces' });
  const r = seeded(33);
  let letters = '';
  for (let i = 0; i < 2 * PIECES; i++) letters += 'ACGT'[Math.floor(r() * 4)];
  const s = b.data.sequence('s', { alphabet: 'dna', letters });
  const start = new Int32Array(PIECES);
  const end = new Int32Array(PIECES);
  for (let i = 0; i < PIECES; i++) {
    start[i] = 2 * i;
    end[i] = 2 * i + 1;
  }
  const pieces = b.data.ranges('pieces', { target: s, rows: { start, end } });
  const signal = b.data.values('signal', {
    target: s,
    values: Float32Array.from({ length: 2 * PIECES }, (_, i) => Math.sin(i / 50)),
  });
  b.views.view('main', (v) => {
    const sec = v.section('s', { of: s, direction: 'reversed', pieces });
    sec.track('letters', { label: 'letters' }).layer('letters', { type: 'letters' });
    sec.track('signal', { label: 'signal' }).layer('values', { type: 'values', data: signal });
  });
  return b.build();
}

/**
 * p53 with the small layers: a ruler, secondary structure (helices, strands, turns), variants (three at one site), the
 * conservation as a line and as an area, the domains dashed with an arrowhead and as highlights.
 */
export function smallLayersSpec(): Node {
  const b = createSeqSpec({ title: 'p53 small layers' });
  const p53 = b.data.sequence('p53', { alphabet: 'protein', letters: P53 });
  const d = domains(b, p53);
  const c = b.data.values('conservation', { target: p53, values: conservation() });
  const ss = b.data.ranges('ss', { target: p53, rows: SecondaryRows });
  const vs = b.data.ranges('variants', { target: p53, rows: VariantRows });
  b.views.view('main', (v) => {
    const s = v.section('protein', { of: p53 });
    s.track('ruler', { label: 'Position' }).layer('ruler', { type: 'ruler' });
    s.track('ss', { label: 'Secondary' }).layer('ss', { type: 'secondary', data: ss });
    s.track('variants', { label: 'Variants' }).layer('variants', { type: 'variants', data: vs });
    s.track('line', { label: 'Line' }).layer('line', { type: 'values', data: c, encoding: 'line' });
    s.track('area', { label: 'Area' }).layer('area', { type: 'values', data: c, encoding: 'area' });
    s.track('dashed', { label: 'Dashed' }).layer('dashed', { type: 'ranges', data: d, line: 'dashed', arrow: 'end' });
    s.track('highlight', { label: 'Highlight' }).layer('hl', { type: 'ranges', data: d, line: 'highlight' });
  });
  return b.build();
}

/** Secondary structure on p53: illustrative bounds, `category` the class. */
export const SecondaryRows = {
  start: new Int32Array([20, 40, 60, 100, 120, 160, 180, 250, 300, 330]),
  end: new Int32Array([30, 55, 66, 115, 140, 172, 200, 262, 320, 345]),
  id: ['h1', 's1', 'h2', 's2', 'h3', 's3', 's4', 'h4', 's5', 'h5'],
  category: ['helix', 'strand', 'helix', 'strand', 'helix', 'strand', 'strand', 'helix', 'strand', 'helix'],
};

/** Variants on p53: three at one site (R175), the others alone. */
export const VariantRows = {
  start: new Int32Array([10, 72, 130, 174, 174, 174, 244, 247, 272, 281, 306, 341]),
  end: new Int32Array([11, 73, 131, 175, 175, 175, 245, 248, 273, 282, 307, 342]),
  id: ['v1', 'v2', 'v3', 'v4', 'v5', 'v6', 'v7', 'v8', 'v9', 'v10', 'v11', 'v12'],
  category: [
    'missense',
    'missense',
    'nonsense',
    'missense',
    'missense',
    'synonymous',
    'missense',
    'missense',
    'frameshift',
    'missense',
    'nonsense',
    'missense',
  ],
};

/** p53 as text: its letters, its two domains under them (a band that collapses where no domain is). */
export function p53TextSpec(section: Record<string, unknown> = {}): Node {
  const b = createSeqSpec({ title: 'p53 text' });
  const p53 = b.data.sequence('p53', { alphabet: 'protein', letters: P53 });
  const d = b.data.ranges('domains', {
    target: p53,
    rows: {
      start: Int32Array.from([0, 94]),
      end: Int32Array.from([61, 292]),
      id: categoryColumn(['TAD', 'DBD']),
    },
  });
  b.views.view('main', { layout: 'text' }, (v) => {
    const s = v.section('protein', { of: p53, ...section });
    s.track('letters', { label: 'p53' }).layer('letters', { type: 'letters' });
    s.track('domains', { label: 'Domains' }).layer('ranges', { type: 'ranges', data: d });
  });
  return b.build();
}
