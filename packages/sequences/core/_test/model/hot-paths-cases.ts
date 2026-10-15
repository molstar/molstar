/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The model's hot-path cases, for their heap guard and its quick twin (hot-paths.heap.ts, hot-paths.quick.ts) (M1, M2a,
// M2b; MD2's coordinate maps; MD3's indexes): each is a step of the paths a hover or a draw takes through the model, over
// fixtures built here. `Loci.forEachLocation` is a step of a whole walk (a million cells), the others of a few lookups.

import {
  AlignedLocation,
  type AlignedLoci,
  type AlignmentId,
  type Annotation,
  type AnnotationId,
  Annotations,
  type ColumnIndex,
  ColumnLocation,
  CoordinateSystemKind,
  Coords,
  Fasta,
  GeneticCode,
  ItemLocation,
  KeyKind,
  type KeyTarget,
  Loci,
  type MapHit,
  type Mapping,
  Mappings,
  type NameId,
  type NameNumbering,
  Newick,
  Numbering,
  type NumberingInput,
  OrderedSet,
  Pairs,
  type PathStep,
  Paths,
  type PositionLabel,
  type RowIndex,
  RowLocation,
  type SeqIndex,
  SeqLocation,
  type SeqLoci,
  SeqSets,
  SortedRanges,
  StringTable,
  Trees,
  type UnitId,
} from '../../src/model/index.js';
import { pf00042, SH2_A3M } from './alignment-fixtures.js';
import { contacts, iep1 } from './annotation-fixtures.js';
import { CDS, GENE, PROTEIN } from './codon-fixtures.js';
import { a3n1, brs1 } from './fixtures.js';
import { gene as dmd } from './support/generators.js';

/** Steps in the shorter run; the longer runs four times as many. */
const N = 250_000;
/** Accumulators stay small integers (Smis): a double is boxed when stored, and in V8's lower tiers. */
const MASK = 0xfffff;

/** The steps' results, shared by every path and kept small integers, so a step boxes no numbers. */
let acc = 0;
const ctx = { n: 0 };

const strings = new StringTable();
const build = (input: NumberingInput, length: number) => Numbering.build(input, length, strings, [], '') as Numbering;
const L = 400;
/**
 * Suffixed labels of both kinds: a letter (`3A`) and a multi-digit `.n` (`7.12`), whose rank the pre-CM-09 code made
 * from a string a position (LA-14: a letter or a one-digit `.n` made none, and the case was blind).
 */
const labels: (number | string)[] = Array.from({ length: L }, (_, i) =>
  i % 10 === 3 ? `${i}A` : i % 10 === 7 ? `${i}.12` : i + 1,
);
const offset = build({ name: 'o', method: 'offset', first: -47, skipZero: true }, L);
const segments = build(
  {
    name: 's',
    method: 'segments',
    segments: [
      [0, 100, 150],
      [160, 400, 240],
    ],
  },
  L,
);
const labelNumbering = build({ name: 'l', method: 'labels', labels }, L);
const names = build({ name: 'n', method: 'names', names: labels.map(String) }, L) as NameNumbering;
const suffixA = strings.find('A');
const nameIds = Int32Array.from(labels, (x) => strings.find(String(x)));
const out: PositionLabel = { num: 0, suffix: -1 as NameId, name: -1 as NameId };
const range = new Int32Array(2);
const beg: PositionLabel = { num: 10, suffix: -1 as NameId, name: -1 as NameId };
const end: PositionLabel = { num: 300, suffix: -1 as NameId, name: -1 as NameId };
/** An endpoint the numbering lacks (index 3 is "3A", so no 4), placed by order (CM-09: no string a position). */
const absent: PositionLabel = { num: 4, suffix: -1 as NameId, name: -1 as NameId };

const sparse = OrderedSet.ofSortedArray(Int32Array.from({ length: 5000 }, (_, i) => i * 3 + (i % 2)));
const dense = OrderedSet.ofRange(1000, 900_000);
const observed = Int32Array.from({ length: 200 }, (_, i) => i * 10 + (i % 2) * 5);

const { set } = SeqSets.build(brs1());
const seqLoci: SeqLoci = {
  kind: 'seq-loci',
  set,
  elements: [0, 2, 4].map((u) => ({
    unit: u as UnitId,
    indices: OrderedSet.ofSortedArray(Int32Array.of(3, 9, 27, 81)),
  })),
};
const rect: AlignedLoci = {
  kind: 'aligned-loci',
  set,
  alignment: 0 as AlignmentId,
  elements: {
    rows: OrderedSet.ofRange(0 as RowIndex, 100 as RowIndex),
    columns: OrderedSet.ofRange(0 as ColumnIndex, 10_000 as ColumnIndex),
  },
};
const seq = SeqLocation.create(set);
const cell = AlignedLocation.create(set);
const target: KeyTarget = { kind: -1, id: -1, row: -1 };
const unitKey = set.units.key[4] as NameId;

const globin = pf00042();
const rowName = 'globin/17';
const tree = Newick.parse(`(${globin.map((r, i) => (i % 2 ? `${r.label})` : `(${r.label}`)).join(',')});`).tree;
const msa = SeqSets.build({
  molecules: [],
  alignments: [
    {
      key: 'globin',
      flavor: 'aligned',
      gapped: globin.map((r) => r.gapped),
      labels: globin.map((r) => r.label),
      tree: tree
        ? {
            ...tree,
            leafRow: Trees.matchLeaves(
              tree,
              globin.map((r) => r.label),
              'labels',
            ),
          }
        : undefined,
    },
    Fasta.alignment(Fasta.parse(SH2_A3M).records, 'sh2', 'a3m').input,
    {
      key: 'lazy',
      flavor: 'aligned',
      gapped: [],
      lazy: {
        rowCount: 100_000,
        columnCount: globin[0].gapped.length,
        rows: new Map(globin.map((r, i) => [i * 3001, { label: r.label, gapped: r.gapped }])),
        labelOf: () => null,
      },
    },
  ],
}).set.alignments;
const [aligned, a3m, lazy] = msa.map((a) => a as NonNullable<typeof a>);
const rowMap = new Int32Array(1000);

// M2b: PF00042's HBA_HUMAN row bound to 1A3N's α entity (copies A and C), and β's entity bound to row 1 (the row on the
// map's `to` side); a path from a UniProt-like sequence's own unit to unit C (A to C breaks the copy rule).
const HBA = a3n1().molecules[0].residues as string;
const bound = SeqSets.build({
  ...a3n1(),
  molecules: [...a3n1().molecules, { key: 'P69905', residues: HBA }],
  units: [...(a3n1().units ?? []), { key: 'P69905/u', molecule: 'P69905' }],
  alignments: [{ key: 'globin', flavor: 'aligned', gapped: globin.map((r) => r.gapped) }],
  mappings: [
    { key: 'hba', from: 'globin/0', to: '1a3n/entity:1', method: 'segments', segments: [[0, 25, null]] },
    { key: 'hbb', from: '1a3n/entity:2', to: 'globin/1', method: 'segments', segments: [[0, 0, 60]] },
    { key: 'uniprot', from: 'P69905', to: '1a3n/entity:1', method: 'identity' },
  ],
}).set;
const hba = bound.mappings[0] as Mapping;
const steps: PathStep[] = [];
const unitC = { kind: CoordinateSystemKind.unit, id: 1, row: -1 };
const uniprot = { kind: CoordinateSystemKind.unit, id: 8, row: -1 };
Mappings.path(bound, uniprot, unitC, steps);
const AtoC = Mappings.compose(bound, uniprot, steps);
const mapHit: MapHit = { element: -1, phase: -1 };
function hitSum(h: MapHit, ctx: { n: number }) {
  ctx.n = (ctx.n + h.element + h.phase) & MASK;
}
const boundColumns = (bound.alignments[0] as NonNullable<(typeof bound.alignments)[number]>).columnCount;
const hover = {
  seq: SeqLocation.create(bound),
  seqOut: SeqLocation.create(bound),
  column: ColumnLocation.create(bound),
  row: RowLocation.create(bound),
  rowOut: RowLocation.create(bound),
  cell: AlignedLocation.create(bound),
};

function count(_: unknown, ctx: { n: number }) {
  ctx.n = (ctx.n + 1) & MASK;
}

// M3: 1IEP's contacts (five ligand poses) as pairs; ranges (with values) and dense values over 1BRS.
const iep = SeqSets.build({ ...iep1(), annotations: [contacts()] }).set;
const pairs = iep.annotations[0] as Annotation;
const ranges = Array.from({ length: 400 }, (_, k) => ({ t: k % 3 ? '1brs/A' : '1brs/entity:1', s: (k * 7) % 100 }));
const brs = SeqSets.build({
  ...brs1(),
  annotations: [
    {
      key: 'r',
      geometry: 'ranges',
      columns: {
        target: ranges.map((r) => r.t),
        start: ranges.map((r) => r.s),
        stop: ranges.map((r) => r.s + 1 + (r.s % 9)),
        id: ranges.map((_, k) => `i${k % 50}`),
        value: ranges.map((_, k) => k),
      },
    },
    { key: 'v', geometry: 'values', target: '1brs/D', columns: { value: Array.from({ length: 89 }, (_, i) => i) } },
  ],
}).set;
const [rangesAnn, valuesAnn] = brs.annotations as Annotation[];
const itemQuery = { seq: SeqLocation.create(brs), item: ItemLocation.create(brs) };
function rowSum(row: number, ctx: { n: number }) {
  ctx.n = (ctx.n + row) & MASK;
}
function pairSum(row: number, side: number, ctx: { n: number }) {
  ctx.n = (ctx.n + row + side) & MASK;
}
function partnerSum(unit: number, index: number, row: number, ctx: { n: number }) {
  ctx.n = (ctx.n + unit + index + row) & MASK;
}

// C1: case 14's protein over its gene; one codon step each way, composed.
const coding = SeqSets.build({
  molecules: [
    { key: 'gene', alphabet: 'dna', residues: GENE, ownUnit: true },
    { key: 'protein', alphabet: 'protein', residues: PROTEIN, ownUnit: true },
  ],
  mappings: [{ key: 'c', from: 'protein', to: 'gene', method: 'codon', cds: CDS }],
}).set;
// Own units in molecule order: gene/unit 0, protein/unit 1.
const [gene, protein] = [0, 1].map((u) => ({ kind: CoordinateSystemKind.unit, id: u, row: -1 }));
const codonSteps: PathStep[] = [];
if (Mappings.path(coding, protein, gene, codonSteps) !== 3) throw new Error('no codon path');
const down = Mappings.compose(coding, protein, codonSteps);
Mappings.path(coding, gene, protein, codonSteps);
const up = Mappings.compose(coding, gene, codonSteps);
const codonHover = { l: SeqLocation.create(coding), out: SeqLocation.create(coding) };
function runSum(start: number, stop: number, ctx: { n: number }) {
  ctx.n = (ctx.n + start + stop) & MASK;
}

// MD2: the DMD-sized gene (79 exons over 2.2M nucleotides) and its protein: lookups both ways and range walks.
const dmdGene = dmd();
const dmdSet = SeqSets.build({
  molecules: [
    { key: 'dmd', alphabet: 'dna', residues: dmdGene.dna },
    { key: 'dystrophin', alphabet: 'protein', residues: dmdGene.protein },
  ],
  mappings: [{ key: 'c', from: 'dystrophin', to: 'dmd', method: 'codon', cds: dmdGene.exons.map(([s, t]) => [s, t]) }],
}).set;
const dmdMap = (dmdSet.mappings[0] as Mapping).map;
const span = dmdMap.toStop - dmdMap.toStart;

// MD3: the indexes on the hover path. A pairs annotation of 64k contacts on one 20k unit, some whole-unit endpoints
// (touching's binary searches; endUnits), and 1,000 chains each bound to its own row.
const BIG = 20_000;
const bigPairs = SeqSets.build({
  molecules: [
    { key: 'u', residues: 'A'.repeat(BIG), ownUnit: true },
    { key: 'w', residues: 'AAAA', ownUnit: true },
  ],
  annotations: [
    {
      key: 'c',
      geometry: 'pairs',
      target: 'u/unit',
      columns: {
        start: Array.from({ length: 64_000 }, (_, i) => (i * 7919) % BIG),
        stop: Array.from({ length: 64_000 }, (_, i) => ((i * 7919) % BIG) + 1 + (i % 9)),
        b_target: Array.from({ length: 64_000 }, (_, i) => (i % 5 ? 'w/unit' : 'u/unit')),
        b_start: Array.from({ length: 64_000 }, (_, i) => (i % 101 === 0 ? null : i % 4)),
      },
    },
  ],
}).set.annotations[0] as Annotation;
const CHAIN = 'MKVLAAGIVGLPNVGKSTLFNALTKAGIEAANYPFCTIEPNTGVVPMPDP';
const chainSet = SeqSets.build({
  molecules: [{ key: 'M', residues: CHAIN }],
  units: Array.from({ length: 1000 }, (_, k) => ({ key: `c${k}`, molecule: 'M' })),
  alignments: [{ key: 'aln', flavor: 'aligned', gapped: Array.from({ length: 1000 }, () => `-${CHAIN}`) }],
  mappings: Array.from({ length: 1000 }, (_, k) => ({
    key: `m${k}`,
    from: `aln/${k}`,
    to: `c${k}`,
    method: 'identity' as const,
  })),
}).set;
const chainHover = { seq: SeqLocation.create(chainSet), row: RowLocation.create(chainSet) };
function unitSum(unit: number, ctx: { n: number }) {
  ctx.n = (ctx.n + unit) & MASK;
}

/** A guarded step: its name, the steps in the guard's shorter run and the step. */
export type HotPath = [name: string, n: number, step: (i: number) => void];

export const PATHS: HotPath[] = [
  [
    'MD3: touching over 64k pairs with whole-unit endpoints, and endUnits',
    N,
    (i) => {
      Pairs.touching(bigPairs, 0 as UnitId, ((i * 31) % BIG) as SeqIndex, pairSum, ctx);
      if ((i & 3) === 0) Pairs.endUnits(bigPairs, i % 64_000, i & 1, unitSum, ctx);
    },
  ],
  [
    'MD3: toRowLocations over 1,000 bound chains',
    N,
    (i) => {
      SeqLocation.set(chainHover.seq, i % 1000, i % 50);
      SeqLocation.toRowLocations(chainHover.seq, 0 as AlignmentId, chainHover.row, count, ctx);
    },
  ],
  [
    'coordinate maps: toAt, forEachFromAt, forEachTo (runs and boundaries) and forEachFrom on a 79-exon codon map',
    N,
    (i) => {
      const r = (i * 7919) % 3685;
      const j = Coords.toAt(dmdMap, r, i % 3);
      acc = (acc + j) & MASK;
      Coords.forEachFromAt(dmdMap, j, mapHit, hitSum, ctx);
      Coords.forEachFromAt(dmdMap, dmdMap.toStart + ((i * 104729) % span), mapHit, hitSum, ctx);
      if ((i & 3) === 0) Coords.forEachTo(dmdMap, r, r + (i % 40), runSum, ctx);
      if ((i & 3) === 1)
        Coords.forEachFrom(
          dmdMap,
          dmdMap.toStart + ((i * 31) % span),
          dmdMap.toStart + ((i * 31) % span) + 5000,
          runSum,
          ctx,
        );
    },
  ],
  [
    'codons: SeqLocation.toSystemEach both ways, Mappings.forEachRun, GeneticCode.translate and isStart',
    N,
    (i) => {
      SeqLocation.set(codonHover.l, protein.id, i % 20);
      SeqLocation.toSystemEach(codonHover.l, down, codonHover.out, count, ctx);
      SeqLocation.set(codonHover.l, gene.id, i % 120);
      SeqLocation.toSystemEach(codonHover.l, up, codonHover.out, count, ctx);
      if (i % 4 === 0) Coords.forEachTo(down.map, i % 16, (i % 16) + 4, runSum, ctx);
      if (i % 4 === 1) Coords.forEachFrom(up.map, i % 110, (i % 110) + 9, runSum, ctx);
      if (i % 4 === 2) Coords.forEachTo(down.map, i % 21, i % 21, runSum, ctx);
      const k = i % 118;
      const a = GENE.charCodeAt(k);
      const b = GENE.charCodeAt(k + 1);
      const c = GENE.charCodeAt(k + 2);
      acc = (acc + GeneticCode.translate(1 + (i & 1), a, b, c) + (GeneticCode.isStart(11, a, b, c) ? 1 : 0)) & MASK;
    },
  ],
  [
    'offset indexOf and labelOf',
    N,
    (i) => {
      acc = (acc + Numbering.indexOf(offset, (i % 300) - 50, -1 as NameId)) & MASK;
      if (Numbering.labelOf(offset, (i % L) as SeqIndex, out)) acc = (acc + out.num) & MASK;
    },
  ],
  [
    'segments indexOf and labelOf',
    N,
    (i) => {
      acc = (acc + Numbering.indexOf(segments, 90 + (i % 600), -1 as NameId)) & MASK;
      if (Numbering.labelOf(segments, (i % L) as SeqIndex, out)) acc = (acc + out.num) & MASK;
    },
  ],
  [
    'labels indexOf with suffixes, labelOf and rangeOf',
    N,
    (i) => {
      const k = i % L;
      acc = (acc + Numbering.indexOf(labelNumbering, k, (k % 10 === 3 ? suffixA : -1) as NameId)) & MASK;
      if (Numbering.labelOf(labelNumbering, k as SeqIndex, out)) acc = (acc + out.num + out.suffix) & MASK;
      if (i % 16 === 0 && Numbering.rangeOf(labelNumbering, beg, end, range)) acc = (acc + range[1]) & MASK;
      if (i % 16 === 8 && Numbering.rangeOf(labelNumbering, absent, end, range)) acc = (acc + range[0]) & MASK;
    },
  ],
  [
    'names indexOfName',
    N,
    (i) => {
      acc = (acc + Numbering.indexOfName(names, nameIds[i % L] as NameId)) & MASK;
    },
  ],
  [
    'OrderedSet has, indexOf, lowerBound and getAt on both forms',
    N,
    (i) => {
      const x = (i * 7) % 16000;
      acc = (acc + (OrderedSet.has(sparse, x) ? 1 : 0)) & MASK;
      acc =
        (acc + OrderedSet.indexOf(sparse, x) + OrderedSet.lowerBound(sparse, x) + OrderedSet.getAt(sparse, i % 5000)) &
        MASK;
      acc = (acc + (OrderedSet.has(dense, x) ? 1 : 0)) & MASK;
      acc =
        (acc + OrderedSet.indexOf(dense, x) + OrderedSet.lowerBound(dense, x) + OrderedSet.getAt(dense, i % 5000)) &
        MASK;
      acc = (acc + (SortedRanges.has(observed, x % 2000) ? 1 : 0)) & MASK;
    },
  ],
  [
    'locations: set, Loci.has, SeqSet.byKey, byName',
    N,
    (i) => {
      SeqLocation.set(seq, i % 6, i % 100);
      AlignedLocation.set(cell, 0, i % 120, i % 12_000);
      acc = (acc + (Loci.has(seqLoci, seq) ? 1 : 0) + (Loci.has(rect, cell) ? 1 : 0)) & MASK;
      if (set.byKey(unitKey, target) && target.kind === KeyKind.unit) acc = (acc + target.id) & MASK;
      if (bound.byName(rowName, target)) acc = (acc + target.row) & MASK;
    },
  ],
  [
    'Loci.forEachLocation over a million cells',
    2,
    () => {
      // A step is a whole walk: runs of 2 and 8 walks.
      Loci.forEachLocation(rect, cell, count, ctx);
      Loci.forEachLocation(seqLoci, seq, count, ctx);
      acc = (acc + ctx.n) & MASK;
    },
  ],
  [
    'colToPos and colOf on aligned, A3M and lazy rows',
    N,
    (i) => {
      const C = aligned.columnCount;
      const r = (i % 32) as RowIndex;
      acc = (acc + aligned.rows.colToPos(r, ((i * 7) % C) as ColumnIndex)) & MASK;
      acc = (acc + a3m.rows.colToPos((i & 3) as RowIndex, (i % 30) as ColumnIndex)) & MASK;
      acc = (acc + lazy.rows.colToPos(((i % 64) * 1500) as RowIndex, (i % C) as ColumnIndex)) & MASK;
      acc = (acc + aligned.rows.colOf(r, (i * 5) % 60) + a3m.rows.colOf((i & 3) as RowIndex, i % 30)) & MASK;
    },
  ],
  [
    'posToCol into a reused array, posMap hits and cladeRows',
    N,
    (i) => {
      const t = aligned.tree;
      if (i % 8 === 0 && a3m.rows.posToCol((i & 3) as RowIndex, rowMap)) acc = (acc + rowMap[i % 28]) & MASK;
      const map = aligned.rows.posMap((i % 32) as RowIndex);
      if (map) acc = (acc + map[i % 50]) & MASK;
      if (t) Trees.cladeRows(t, i % t.count, count, ctx);
    },
  ],
  [
    'mapping lookups and SeqLocation.toSystem along a composed path',
    N,
    (i) => {
      const k = i % 141;
      acc = (acc + Coords.forEachFromAt(hba.map, k, mapHit, hitSum, ctx) + hba.status[i % 111]) & MASK;
      acc = (acc + Coords.toAt(hba.map, i % 120, 0)) & MASK;
      SeqLocation.set(hover.seq, uniprot.id, k);
      if (SeqLocation.toSystem(hover.seq, AtoC, hover.seqOut)) acc = (acc + hover.seqOut.index) & MASK;
    },
  ],
  [
    'hover conversions: column, cell and row to units, and back',
    N,
    (i) => {
      ColumnLocation.set(hover.column, 0, i % boundColumns);
      ColumnLocation.toSeqLocations(hover.column, hover.seqOut, count, ctx);
      AlignedLocation.set(hover.cell, 0, i & 7, i % boundColumns);
      if (AlignedLocation.toRow(hover.cell, hover.row)) {
        RowLocation.toSeqLocations(hover.row, hover.seqOut, count, ctx);
        if (RowLocation.toAligned(hover.row, hover.cell)) acc = (acc + hover.cell.column) & MASK;
      }
      SeqLocation.set(hover.seq, i & 1, i % 141);
      SeqLocation.toRowLocations(hover.seq, 0 as AlignmentId, hover.rowOut, count, ctx);
      SeqLocation.set(hover.seq, 2 + (i & 1), i % 146);
      SeqLocation.toRowLocations(hover.seq, 0 as AlignmentId, hover.rowOut, count, ctx);
    },
  ],
  [
    'annotation interval queries and valueAt, sparse and dense',
    N,
    (i) => {
      const at = i % 110;
      Annotations.itemsOver(rangesAnn, i & 1, at, at + (i & 7) + 1, rowSum, ctx);
      acc =
        (acc + (Annotations.valueAt(rangesAnn, i & 1, at) | 0) + (Annotations.valueAt(valuesAnn, 0, i % 89) | 0)) &
        MASK;
    },
  ],
  [
    'pair lookups: touching and partnersOf',
    N,
    (i) => {
      Pairs.touching(pairs, 0 as UnitId, (20 + (i % 160)) as SeqIndex, pairSum, ctx);
      Pairs.touching(pairs, (1 + (i % 5)) as UnitId, (i % 37) as SeqIndex, pairSum, ctx);
      if ((i & 15) === 0) Pairs.partnersOf(pairs, (1 + (i % 5)) as UnitId, partnerSum, ctx);
    },
  ],
  [
    'SeqLocation.itemsAt over a unit and its molecule',
    N,
    (i) => {
      SeqLocation.set(itemQuery.seq, i % 3, i % 110);
      SeqLocation.itemsAt(itemQuery.seq, 0 as AnnotationId, itemQuery.item, count, ctx);
    },
  ],
];

/** The steps' results, small integers: a guard over steps that read nothing would pass blind. */
export const sums = () => acc + ctx.n;

/** The fixtures' own checks: what each shared fixture is (name, what it reads, what it must be). */
export const FIXTURES: [name: string, actual: () => unknown, expected: unknown][] = [
  [
    'runs over a bound row and a composed path',
    () => {
      const b = Mappings.bindingsOf(bound, 0);
      return [steps.length, Paths.at(AtoC, 50), b?.unit.length, Array.from(b?.rowIsFrom ?? [])];
    },
    [3, 50, 4, [1, 1, 0, 0]],
  ],
  ['runs over a DMD-sized codon map', () => [dmdMap.count, dmdMap.toStop - dmdMap.toStart > 2_000_000], [79, true]],
  [
    'runs over built annotations',
    () => [pairs.count, rangesAnn.count, rangesAnn.slots.unit.length, valuesAnn.dense?.[0].length],
    [83, 400, 2, 89],
  ],
  [
    'runs over built alignments and a tree',
    () => [aligned.tree?.count, a3m.columnCount, lazy.rows.count],
    [49, 30, 32],
  ],
  [
    'places a label endpoint the numbering lacks by order',
    () => [Numbering.rangeOf(labelNumbering, absent, end, range), range[0] > 0],
    [true, true],
  ],
];
