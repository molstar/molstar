/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Locations and loci (seqspec:model#locations-and-loci). A location is a few numbers in a mutable object the caller
// owns and reuses; loci are many locations of one kind, stored compactly. Conversions follow alignments and mappings:
// location-level ones are the hover path and allocate nothing; loci-level ones are edges and report into a reused
// `Conversion`.

import { type Annotation, Annotations, Pairs } from './annotation.js';
import type {
  AlignmentId,
  AnnotationId,
  BoundaryIndex,
  ColumnIndex,
  ItemIndex,
  RowIndex,
  SeqIndex,
  UnitId,
} from './common.js';
import { CoordinateSystemKind, Coords, type MapHit } from './coordinates.js';
import { type ComposedPath, Mappings, Paths } from './mapping.js';
import { OrderedSet, SortedArray } from './ordered-set.js';
import { KeyKind, type KeyTarget, type SeqSet } from './set.js';

export interface SeqLocation {
  readonly kind: 'seq';
  set: SeqSet;
  unit: UnitId;
  index: SeqIndex;
}
export interface BoundaryLocation {
  readonly kind: 'boundary';
  set: SeqSet;
  unit: UnitId;
  at: BoundaryIndex;
}
/** A residue of an alignment row (unaligned). */
export interface RowLocation {
  readonly kind: 'row';
  set: SeqSet;
  alignment: AlignmentId;
  row: RowIndex;
  index: SeqIndex;
}
/** A cell: a row at a column. */
export interface AlignedLocation {
  readonly kind: 'aligned';
  set: SeqSet;
  alignment: AlignmentId;
  row: RowIndex;
  column: ColumnIndex;
}
export interface ColumnLocation {
  readonly kind: 'column';
  set: SeqSet;
  alignment: AlignmentId;
  column: ColumnIndex;
}
export interface ItemLocation {
  readonly kind: 'item';
  set: SeqSet;
  annotation: AnnotationId;
  item: ItemIndex;
}
export type Location = SeqLocation | BoundaryLocation | RowLocation | AlignedLocation | ColumnLocation | ItemLocation;

export interface SeqLoci {
  readonly kind: 'seq-loci';
  readonly set: SeqSet;
  readonly elements: readonly { readonly unit: UnitId; readonly indices: OrderedSet<SeqIndex> }[];
}
export interface BoundaryLoci {
  readonly kind: 'boundary-loci';
  readonly set: SeqSet;
  readonly elements: readonly { readonly unit: UnitId; readonly indices: OrderedSet<BoundaryIndex> }[];
}
export interface RowLoci {
  readonly kind: 'row-loci';
  readonly set: SeqSet;
  readonly alignment: AlignmentId;
  readonly elements: readonly { readonly row: RowIndex; readonly indices: OrderedSet<SeqIndex> }[];
}
/** A row band or a range over many rows: one rectangle. */
export interface AlignedRect {
  readonly rows: OrderedSet<RowIndex>;
  readonly columns: OrderedSet<ColumnIndex>;
}
export interface AlignedLoci {
  readonly kind: 'aligned-loci';
  readonly set: SeqSet;
  readonly alignment: AlignmentId;
  readonly elements: readonly { readonly row: RowIndex; readonly columns: OrderedSet<ColumnIndex> }[] | AlignedRect;
}
export interface ColumnLoci {
  readonly kind: 'column-loci';
  readonly set: SeqSet;
  readonly alignment: AlignmentId;
  readonly columns: OrderedSet<ColumnIndex>;
}
export interface ItemLoci {
  readonly kind: 'item-loci';
  readonly set: SeqSet;
  readonly elements: readonly { readonly annotation: AnnotationId; readonly items: OrderedSet<ItemIndex> }[];
}
export type Loci = SeqLoci | BoundaryLoci | RowLoci | AlignedLoci | ColumnLoci | ItemLoci;

type Rec = Record<string, number>;
type El = Readonly<Record<string, OrderedSet>>;
type Grouped = Exclude<Loci, ColumnLoci>;

/** Per location kind: its fields after `set`; the loci kind; the element's owner and set fields. */
const FIELDS = {
  seq: ['unit', 'index'],
  boundary: ['unit', 'at'],
  row: ['alignment', 'row', 'index'],
  aligned: ['alignment', 'row', 'column'],
  column: ['alignment', 'column'],
  item: ['annotation', 'item'],
} as const;
const LOCI_OF = {
  seq: 'seq-loci',
  boundary: 'boundary-loci',
  row: 'row-loci',
  aligned: 'aligned-loci',
  column: 'column-loci',
  item: 'item-loci',
} as const;
/** [owner, set, the location's value field]. */
const SHAPE = {
  'seq-loci': ['unit', 'indices', 'index'],
  'boundary-loci': ['unit', 'indices', 'at'],
  'row-loci': ['row', 'indices', 'index'],
  'aligned-loci': ['row', 'columns', 'column'],
  'item-loci': ['annotation', 'items', 'item'],
} as const;

function locationKind<L extends Location>(kind: L['kind']) {
  const f: readonly string[] = FIELDS[kind];
  return {
    create(set: SeqSet): L {
      const l: Record<string, unknown> = { kind, set };
      for (const k of f) l[k] = -1;
      return l as L;
    },
    /** The fields in `FIELDS` order: (unit, index), (alignment, row, column), … */
    set(out: L, a: number, b: number, c = -1): L {
      const o = out as unknown as Rec;
      o[f[0]] = a;
      o[f[1]] = b;
      if (f.length > 2) o[f[2]] = c;
      return out;
    },
    areEqual(x: L, y: L): boolean {
      if (x.set.lineage !== y.set.lineage) return false;
      for (let i = 0; i < f.length; i++) if ((x as unknown as Rec)[f[i]] !== (y as unknown as Rec)[f[i]]) return false;
      return true;
    },
  };
}

export const SeqLocation = { ...locationKind<SeqLocation>('seq'), toSystem, toSystemEach, toRowLocations, itemsAt };
export const BoundaryLocation = locationKind<BoundaryLocation>('boundary');
export const RowLocation = { ...locationKind<RowLocation>('row'), toAligned, toSeqLocations: rowToSeq };
export const AlignedLocation = { ...locationKind<AlignedLocation>('aligned'), toRow };
export const ColumnLocation = { ...locationKind<ColumnLocation>('column'), toSeqLocations: columnToSeq };
export const ItemLocation = locationKind<ItemLocation>('item');

// Location-level conversions: the hover path, allocation-free.

/** False at a gap, or a row not held. */
function toRow(l: AlignedLocation, out: RowLocation): boolean {
  const p = l.set.alignments[l.alignment]?.rows.colToPos(l.row, l.column) ?? (-1 as SeqIndex);
  if (p < 0) return false;
  out.set = l.set;
  out.alignment = l.alignment;
  out.row = l.row;
  out.index = p;
  return true;
}

/** False at an A3M insert, or a row not held. */
function toAligned(l: RowLocation, out: AlignedLocation): boolean {
  const map = l.set.alignments[l.alignment]?.rows.posMap(l.row);
  const c = map && l.index >= 0 && l.index < map.length ? map[l.index] : -1;
  if (c < 0) return false;
  out.set = l.set;
  out.alignment = l.alignment;
  out.row = l.row;
  out.column = c as ColumnIndex;
  return true;
}

/** Visits a unit's position when the unit is alive and the index within its range: 1, else 0. */
function emitSeq<C>(
  set: SeqSet,
  unit: number,
  index: number,
  out: SeqLocation,
  cb: (l: SeqLocation, ctx: C) => void,
  ctx: C,
) {
  const u = set.units;
  if (!u.alive[unit] || !(index >= u.rangeStart[unit] && index < u.rangeStop[unit])) return 0;
  out.set = set;
  out.unit = unit as UnitId;
  out.index = index as SeqIndex;
  cb(out, ctx);
  return 1;
}

/** Every bound unit's position of row position `pos`, through the row's bindings. */
function rowPosToSeq<C>(
  set: SeqSet,
  alignment: number,
  row: number,
  pos: number,
  out: SeqLocation,
  cb: (l: SeqLocation, ctx: C) => void,
  ctx: C,
): number {
  const b = Mappings.bindingsOf(set, alignment);
  const k = b ? Mappings.boundRow(b, row) : -1;
  if (!b || k < 0 || pos < 0) return 0;
  let n = 0;
  for (let e = b.start[k]; e < b.start[k + 1]; e++) {
    const map = b.map[e];
    if (b.rowIsFrom[e]) n += emitSeq(set, b.unit[e], Coords.toAt(map, pos, 0), out, cb, ctx);
    else {
      // The row on the map's to side: every element holding the residue (a source used twice gives two).
      const v = begin(set, b.unit[e], out, cb, ctx);
      Coords.forEachFromAt(map, pos, hit, emitHit, v);
      n += end(v);
    }
  }
  return n;
}

/** The state of one map walk on the hover path, by depth: a callback may convert again. */
type Via = { set: SeqSet; unit: number; out: Location; cb: (l: Location, c: unknown) => void; ctx: unknown; n: number };
const vias: Via[] = [];
let viaDepth = 0;
const hit: MapHit = { element: -1, phase: -1 };
function begin<L extends Location, C>(set: SeqSet, unit: number, out: L, cb: (l: L, c: C) => void, ctx: C): Via {
  if (!vias[viaDepth]) vias[viaDepth] = { set, unit, out, cb: noop, ctx, n: 0 };
  const v = vias[viaDepth++];
  v.set = set;
  v.unit = unit;
  v.out = out;
  v.cb = cb as Via['cb'];
  v.ctx = ctx;
  v.n = 0;
  return v;
}
function end(v: Via): number {
  viaDepth--;
  return v.n;
}
function emitHit(h: MapHit, v: Via) {
  v.n += emitSeq(v.set, v.unit, h.element, v.out as SeqLocation, v.cb, v.ctx);
}
function emitRow(h: MapHit, v: Via) {
  (v.out as RowLocation).index = h.element as SeqIndex;
  v.cb(v.out, v.ctx);
  v.n++;
}

/** Every bound unit's position of a row residue; the count visited. */
function rowToSeq<C>(l: RowLocation, out: SeqLocation, cb: (l: SeqLocation, ctx: C) => void, ctx: C): number {
  return rowPosToSeq(l.set, l.alignment, l.row, l.index, out, cb, ctx);
}

/** Every bound unit's position at this column, through each row's binding; gaps give none. */
function columnToSeq<C>(l: ColumnLocation, out: SeqLocation, cb: (l: SeqLocation, ctx: C) => void, ctx: C): number {
  const b = Mappings.bindingsOf(l.set, l.alignment);
  const rows = l.set.alignments[l.alignment]?.rows;
  if (!b || !rows) return 0;
  let n = 0;
  for (let k = 0; k < b.rows.length; k++) {
    const r = b.rows[k];
    n += rowPosToSeq(l.set, l.alignment, r, rows.colToPos(r as RowIndex, l.column), out, cb, ctx);
  }
  return n;
}

/** Whether a location can go along a path: it starts on the location's unit or molecule and ends on a unit. */
function along(l: SeqLocation, path: ComposedPath): boolean {
  const { from, to } = path;
  const own = from.kind === CoordinateSystemKind.unit ? from.id === l.unit : from.id === l.set.units.molecule[l.unit];
  return (
    own &&
    from.kind <= CoordinateSystemKind.unit &&
    to.kind === CoordinateSystemKind.unit &&
    l.index >= 0 &&
    l.index < path.length
  );
}

/** A unit's position along a composed path (from the unit or its molecule) to a unit, a span's first; false when unmapped. */
function toSystem(l: SeqLocation, path: ComposedPath, out: SeqLocation): boolean {
  return along(l, path) && emitSeq(l.set, path.to.id, Paths.at(path, l.index), out, noop, null) === 1;
}

/**
 * Every position along a path: a residue's three nucleotides through a codon step (seqspec:model#codon-maps); the
 * count.
 */
function toSystemEach<C>(
  l: SeqLocation,
  path: ComposedPath,
  out: SeqLocation,
  cb: (l: SeqLocation, ctx: C) => void,
  ctx: C,
): number {
  if (!along(l, path)) return 0;
  if (!path.heldSpans) return emitSeq(l.set, path.to.id, Paths.at(path, l.index), out, cb, ctx);
  const set = l.set;
  const unit = path.to.id;
  const i = l.index;
  return (
    emitSeq(set, unit, Paths.spanAt(path, i, 0), out, cb, ctx) +
    emitSeq(set, unit, Paths.spanAt(path, i, 1), out, cb, ctx) +
    emitSeq(set, unit, Paths.spanAt(path, i, 2), out, cb, ctx)
  );
}
function noop() {}

/** The residues of an alignment's rows bound to this unit's position; the count visited. */
function toRowLocations<C>(
  l: SeqLocation,
  alignment: AlignmentId,
  out: RowLocation,
  cb: (l: RowLocation, ctx: C) => void,
  ctx: C,
): number {
  const b = Mappings.bindingsOf(l.set, alignment);
  const k = b ? Mappings.boundUnit(b, l.unit) : -1;
  if (!b || k < 0 || l.index < 0) return 0;
  let n = 0;
  out.set = l.set;
  out.alignment = alignment;
  // The unit's entries through the unit index, in row order.
  for (let j = b.unitStart[k], jn = b.unitStart[k + 1]; j < jn; j++) {
    const e = b.byUnit[j];
    out.row = b.rows[b.entryRow[e]] as RowIndex;
    if (b.rowIsFrom[e]) {
      // The unit's position on the map's to side: every row residue holding it.
      const v = begin(l.set, l.unit, out, cb, ctx);
      Coords.forEachFromAt(b.map[e], l.index, hit, emitRow, v);
      n += end(v);
      continue;
    }
    const p = Coords.toAt(b.map[e], l.index, 0);
    if (p < 0) continue;
    out.index = p as SeqIndex;
    cb(out, ctx);
    n++;
  }
  return n;
}

/** The state of one item query (not reentrant for one annotation: a callback must not query it again). */
const items = {
  ann: null as Annotation | null,
  out: null as ItemLocation | null,
  cb: null as ((l: ItemLocation, ctx: unknown) => void) | null,
  ctx: null as unknown,
  n: 0,
};
type ItemQuery = typeof items;
function emitItem(q: ItemQuery, item: number) {
  if (!Annotations.firstVisit(q.ann as Annotation, item)) return;
  (q.out as ItemLocation).item = item as ItemIndex;
  q.cb?.(q.out as ItemLocation, q.ctx);
  q.n++;
}
const visitRow = (row: number, q: ItemQuery) => emitItem(q, (q.ann as Annotation).item[row]);
const visitPair = (row: number, _side: number, q: ItemQuery) => emitItem(q, (q.ann as Annotation).item[row]);

/** Starts an item query: each item is reported once until the next. */
function beginQuery<C>(
  ann: Annotation,
  id: number,
  set: SeqSet,
  out: ItemLocation,
  cb: (l: ItemLocation, ctx: C) => void,
  ctx: C,
) {
  items.ann = ann;
  items.out = out;
  items.cb = cb as ItemQuery['cb'];
  items.ctx = ctx;
  items.n = 0;
  out.set = set;
  out.annotation = id as AnnotationId;
  Annotations.beginItems(ann);
  return items;
}

/**
 * The items of an annotation at a unit's position, each once: rows over the unit or its molecule; for pairs,
 * either endpoint. The count visited. Allocation-free.
 */
function itemsAt<C>(
  l: SeqLocation,
  annotation: AnnotationId,
  out: ItemLocation,
  cb: (l: ItemLocation, ctx: C) => void,
  ctx: C,
): number {
  const ann = l.set.annotations[annotation];
  if (!ann || !l.set.units.alive[l.unit]) return 0;
  const q = beginQuery(ann, annotation, l.set, out, cb, ctx);
  if (ann.pairs) Pairs.touching(ann, l.unit, l.index, visitPair, q);
  else {
    const t = ann.slots;
    const mol = l.set.units.molecule[l.unit];
    for (let s = 0; s < t.unit.length; s++) {
      if (t.unit[s] === l.unit || (t.unit[s] < 0 && t.alignment[s] < 0 && t.molecule[s] === mol)) {
        Annotations.itemsOver(ann, s, l.index, l.index + 1, visitRow, q);
      }
    }
  }
  return q.n;
}

function isRect(e: AlignedLoci['elements']): e is AlignedRect {
  return !Array.isArray(e);
}

/** Per-owner elements; a rectangle expanded per row (edges). */
function elementsOf(loci: Grouped): readonly El[] {
  if (loci.kind !== 'aligned-loci' || !isRect(loci.elements)) return loci.elements as unknown as readonly El[];
  const { rows, columns } = loci.elements;
  const out: El[] = [];
  for (let i = 0; i < OrderedSet.size(rows); i++)
    out.push({ row: OrderedSet.getAt(rows, i), columns } as unknown as El);
  return out;
}

function fromLocation(l: Location): Loci {
  const one = OrderedSet.ofSingleton;
  switch (l.kind) {
    case 'column':
      return { kind: 'column-loci', set: l.set, alignment: l.alignment, columns: one(l.column) };
    case 'aligned':
      return {
        kind: 'aligned-loci',
        set: l.set,
        alignment: l.alignment,
        elements: [{ row: l.row, columns: one(l.column) }],
      };
    case 'row':
      return {
        kind: 'row-loci',
        set: l.set,
        alignment: l.alignment,
        elements: [{ row: l.row, indices: one(l.index) }],
      };
    case 'seq':
      return { kind: 'seq-loci', set: l.set, elements: [{ unit: l.unit, indices: one(l.index) }] };
    case 'boundary':
      return { kind: 'boundary-loci', set: l.set, elements: [{ unit: l.unit, indices: one(l.at) }] };
    case 'item':
      return { kind: 'item-loci', set: l.set, elements: [{ annotation: l.annotation, items: one(l.item) }] };
  }
}

/** Allocation-free. */
function has(loci: Loci, l: Location): boolean {
  if (LOCI_OF[l.kind] !== loci.kind) return false;
  const o = l as unknown as Rec;
  if ('alignment' in loci && loci.alignment !== o.alignment) return false;
  if (loci.kind === 'column-loci') return OrderedSet.has(loci.columns, o.column);
  if (loci.kind === 'aligned-loci' && isRect(loci.elements)) {
    return OrderedSet.has(loci.elements.rows, o.row) && OrderedSet.has(loci.elements.columns, o.column);
  }
  const s = SHAPE[loci.kind];
  const owner = o[s[0]];
  const els = loci.elements as unknown as readonly El[];
  for (let i = 0; i < els.length; i++) if (els[i][s[0]] === owner) return OrderedSet.has(els[i][s[1]], o[s[2]]);
  return false;
}

/** Visits every location through `out` (its kind's); allocation-free. Returns the count. */
function forEachLocation<L extends Location, C>(loci: Loci, out: L, cb: (l: L, ctx: C) => void, ctx: C): number {
  const o = out as unknown as Rec;
  let n = 0;
  (out as { set: SeqSet }).set = loci.set; // Never a stale set's bindings
  if ('alignment' in loci) o.alignment = loci.alignment;
  if (loci.kind === 'column-loci') {
    for (let i = 0, k = OrderedSet.size(loci.columns); i < k; i++, n++) {
      o.column = OrderedSet.getAt(loci.columns, i);
      cb(out, ctx);
    }
    return n;
  }
  if (loci.kind === 'aligned-loci' && isRect(loci.elements)) {
    const { rows, columns } = loci.elements;
    for (let r = 0, nr = OrderedSet.size(rows), nc = OrderedSet.size(columns); r < nr; r++) {
      o.row = OrderedSet.getAt(rows, r);
      for (let c = 0; c < nc; c++, n++) {
        o.column = OrderedSet.getAt(columns, c);
        cb(out, ctx);
      }
    }
    return n;
  }
  const s = SHAPE[loci.kind];
  const els = loci.elements as unknown as readonly El[];
  for (let e = 0; e < els.length; e++) {
    o[s[0]] = els[e][s[0]] as number;
    const v = els[e][s[1]];
    for (let i = 0, k = OrderedSet.size(v); i < k; i++, n++) {
      o[s[2]] = OrderedSet.getAt(v, i);
      cb(out, ctx);
    }
  }
  return n;
}

/** Per-owner elements sorted by owner, empty sets dropped; an aligned result of equal columns is a rectangle when `rect`. */
function grouped<L extends Grouped>(like: L, pairs: Map<number, OrderedSet>, rect: boolean, set: SeqSet): L {
  const s = SHAPE[like.kind];
  const owners = [...pairs.keys()].filter((k) => !OrderedSet.isEmpty(pairs.get(k) as OrderedSet)).sort((a, b) => a - b);
  const els = owners.map((k) => ({ [s[0]]: k, [s[1]]: pairs.get(k) }));
  const first = pairs.get(owners[0]) as OrderedSet;
  const isRectangle =
    rect && owners.length > 0 && owners.every((k) => OrderedSet.areEqual(pairs.get(k) as OrderedSet, first));
  const elements = isRectangle ? { rows: OrderedSet.ofSortedArray(Int32Array.from(owners)), columns: first } : els;
  return { ...like, set, elements };
}

function collect(loci: Grouped, into: Map<number, OrderedSet>) {
  const s = SHAPE[loci.kind];
  for (const e of elementsOf(loci)) {
    const k = e[s[0]] as number;
    const prev = into.get(k);
    into.set(k, prev === undefined ? e[s[1]] : OrderedSet.union(prev, e[s[1]]));
  }
  return into;
}

/** Same kind (and alignment); edges. A rectangle stays one when the result is one. */
function union<L extends Loci>(a: L, b: L): L {
  if ((a as Loci).kind !== b.kind || ('alignment' in a && 'alignment' in b && a.alignment !== b.alignment)) {
    throw new Error(`Loci.union: ${a.kind} and ${b.kind} over different elements.`);
  }
  if (a.kind === 'column-loci') return { ...a, columns: OrderedSet.union(a.columns, (b as ColumnLoci).columns) };
  const ga = a as Grouped;
  const gb = b as Grouped;
  const rect = ga.kind === 'aligned-loci' && (isRect(ga.elements) || isRect((gb as AlignedLoci).elements));
  return grouped(ga, collect(gb, collect(ga, new Map())), rect, a.set) as L;
}

function isEmpty(loci: Loci): boolean {
  if (loci.kind === 'column-loci') return OrderedSet.isEmpty(loci.columns);
  return elementsOf(loci).every((e) => OrderedSet.isEmpty(e[SHAPE[loci.kind][1]]));
}

function areEqual(a: Loci, b: Loci): boolean {
  if (a.kind !== b.kind || a.set.lineage !== b.set.lineage) return false;
  if ('alignment' in a && 'alignment' in b && a.alignment !== b.alignment) return false;
  if (a.kind === 'column-loci') return OrderedSet.areEqual(a.columns, (b as ColumnLoci).columns);
  const ma = collect(a, new Map());
  const mb = collect(b as Grouped, new Map());
  const keys = (m: Map<number, OrderedSet>) => [...m.keys()].filter((k) => !OrderedSet.isEmpty(m.get(k) as OrderedSet));
  const ka = keys(ma);
  return (
    ka.length === keys(mb).length &&
    ka.every((k) => mb.has(k) && OrderedSet.areEqual(ma.get(k) as OrderedSet, mb.get(k) as OrderedSet))
  );
}

/**
 * Loci of `from` in `to`, by key: units, alignments and annotations by their keys, items by their `id`
 * (by index when they have none), rows by index; what is gone or out of bounds is dropped.
 */
function remap<L extends Loci>(loci: L, from: SeqSet, to: SeqSet): L {
  const target: KeyTarget = { kind: 0, id: -1, row: -1 };
  const resolve = (key: number, kind: number) =>
    to.byName(from.strings.get(key), target) && target.kind === kind ? target.id : -1;
  const below = (s: OrderedSet, n: number) => OrderedSet.intersect(s, OrderedSet.ofRange(0, Math.max(0, n)));
  let alignment = -1;
  if ('alignment' in loci) {
    const fa = from.alignments[loci.alignment];
    alignment = fa ? resolve(fa.key, KeyKind.alignment) : -1;
  }
  const ta = to.alignments[alignment];
  if (loci.kind === 'column-loci') {
    return { ...loci, set: to, alignment, columns: below(loci.columns, ta?.columnCount ?? 0) };
  }
  const pairs = new Map<number, OrderedSet>();
  const s = SHAPE[loci.kind];
  for (const e of elementsOf(loci)) {
    const owner = e[s[0]] as number;
    const v = e[s[1]];
    let id = -1;
    let next: OrderedSet = v;
    if (loci.kind === 'seq-loci' || loci.kind === 'boundary-loci') {
      id = resolve(from.units.key[owner], KeyKind.unit);
      const length = id < 0 ? 0 : (to.molecules[to.units.molecule[id]]?.length ?? 0);
      next = below(v, length + (loci.kind === 'boundary-loci' ? 1 : 0));
    } else if (loci.kind === 'item-loci') {
      id = resolve(from.annotations[owner]?.key ?? -1, KeyKind.annotation);
      next = id < 0 ? OrderedSet.Empty : remapItems(v, from, to, owner, id);
    } else if (ta && owner < ta.rowCount) {
      id = owner;
      if (loci.kind === 'aligned-loci') next = below(v, ta.columnCount);
    }
    if (id >= 0) pairs.set(id, pairs.has(id) ? OrderedSet.union(pairs.get(id) as OrderedSet, next) : next);
  }
  const rect = loci.kind === 'aligned-loci' && isRect(loci.elements);
  const out = grouped(loci, pairs, rect, to);
  return ('alignment' in out ? { ...out, alignment } : out) as L;
}

function remapItems(items: OrderedSet, from: SeqSet, to: SeqSet, fromId: number, toId: number): OrderedSet {
  const a = from.annotations[fromId];
  const b = to.annotations[toId];
  if (!a || !b) return OrderedSet.Empty;
  const out: number[] = [];
  OrderedSet.forEach(
    items,
    (item) => {
      const key = a.itemKey[item] ?? -1;
      if (key < 0) {
        if (item < b.itemCount) out.push(item);
        return;
      }
      const k = to.strings.find(from.strings.get(key));
      const j = k < 0 ? -1 : Annotations.itemOf(b, k);
      if (j >= 0) out.push(j);
    },
    null,
  );
  return OrderedSet.ofSortedArray(SortedArray.ofUnsortedArray(Int32Array.from(out)));
}

function lociKind<L extends Loci, Loc extends Location>() {
  return {
    /** Edges. */
    fromLocation: fromLocation as unknown as (l: Loc) => L,
    has: has as (loci: L, l: Loc) => boolean,
    forEachLocation: forEachLocation as <C>(loci: L, out: Loc, cb: (l: Loc, ctx: C) => void, ctx: C) => number,
    union: union as (a: L, b: L) => L,
  };
}

// Loci-level conversions: edges, each reporting into a reused Conversion.

export const ConversionStatus = { exact: 0, partial: 1, ambiguous: 2, unmapped: 3 } as const;
/**
 * `dropped`: source locations with no counterpart (gaps going to residues, inserts going to columns,
 * unbound rows); `multiplied`: source locations with several (copies, several bindings).
 */
export interface Conversion {
  status: number;
  dropped: number;
  multiplied: number;
}
export const Conversion = { create: (): Conversion => ({ status: ConversionStatus.exact, dropped: 0, multiplied: 0 }) };

type Visit<R> = (r: R, acc: Acc) => void;
type Acc = { pairs: Map<number, number[]>; owner: string; value: string };
const collectOne: Visit<Location> = (r, acc) => {
  const o = r as unknown as Rec;
  const list = acc.pairs.get(o[acc.owner]);
  if (list) list.push(o[acc.value]);
  else acc.pairs.set(o[acc.owner], [o[acc.value]]);
};
const single =
  <L, R>(f: (l: L, out: R) => boolean) =>
  (l: L, out: R, cb: Visit<R>, acc: Acc) => {
    if (!f(l, out)) return 0;
    cb(out, acc);
    return 1;
  };

/** Every location of `loci` through `one`, collected by `owner` into sets of `value`; counts into `c`. */
function convert<L extends Location, R extends Location>(
  loci: Loci,
  from: L,
  to: R,
  one: (l: L, out: R, cb: Visit<R>, acc: Acc) => number,
  owner: string,
  value: string,
  c: Conversion,
): Map<number, OrderedSet> {
  const acc: Acc = { pairs: new Map(), owner, value };
  let inputs = 0;
  c.dropped = 0;
  c.multiplied = 0;
  forEachLocation(
    loci,
    from,
    (l) => {
      inputs++;
      const n = one(l, to, collectOne as Visit<R>, acc);
      if (n === 0) c.dropped++;
      else if (n > 1) c.multiplied++;
    },
    null,
  );
  c.status =
    inputs > 0 && acc.pairs.size === 0
      ? ConversionStatus.unmapped
      : c.dropped
        ? ConversionStatus.partial
        : ConversionStatus.exact;
  const sets = new Map<number, OrderedSet>();
  for (const [k, v] of acc.pairs)
    sets.set(k, OrderedSet.ofSortedArray(SortedArray.ofUnsortedArray(Int32Array.from(v))));
  return sets;
}

/** Ambiguous when a bound row the result passes through has several bindings (their union is given). */
function markAmbiguous(c: Conversion, set: SeqSet, alignment: number, rows: Iterable<number>) {
  const b = Mappings.bindingsOf(set, alignment);
  if (!b || c.status === ConversionStatus.unmapped) return;
  for (const r of rows) {
    const k = Mappings.boundRow(b, r);
    if (k >= 0 && b.bindings[k] > 1) c.status = ConversionStatus.ambiguous;
  }
}

/** Cells to residues; gaps dropped. */
function alignedToRowLoci(loci: AlignedLoci, c: Conversion): RowLoci {
  const pairs = convert(
    loci,
    AlignedLocation.create(loci.set),
    RowLocation.create(loci.set),
    single(toRow),
    'row',
    'index',
    c,
  );
  return grouped({ kind: 'row-loci', set: loci.set, alignment: loci.alignment, elements: [] }, pairs, false, loci.set);
}

/** Residues to cells; A3M inserts and rows not held dropped. */
function rowToAlignedLoci(loci: RowLoci, c: Conversion): AlignedLoci {
  const to = AlignedLocation.create(loci.set);
  const pairs = convert(loci, RowLocation.create(loci.set), to, single(toAligned), 'row', 'column', c);
  return grouped(
    { kind: 'aligned-loci', set: loci.set, alignment: loci.alignment, elements: [] },
    pairs,
    false,
    loci.set,
  );
}

/** The columns any cell is in. */
function alignedToColumnLoci(loci: AlignedLoci, c: Conversion): ColumnLoci {
  const all = [...collect(loci, new Map()).values()];
  const columns = all.reduce((a, b) => OrderedSet.union(a, b), OrderedSet.Empty) as OrderedSet<ColumnIndex>;
  Object.assign(c, Conversion.create());
  return { kind: 'column-loci', set: loci.set, alignment: loci.alignment, columns };
}

/** The columns' cells in `rows`: one rectangle. */
function columnToAlignedLoci(loci: ColumnLoci, rows: OrderedSet<RowIndex>, c: Conversion): AlignedLoci {
  Object.assign(c, Conversion.create());
  return { kind: 'aligned-loci', set: loci.set, alignment: loci.alignment, elements: { rows, columns: loci.columns } };
}

/** Row residues to the units their rows are bound to; a molecule-level binding gives every unit over it. */
function rowToSeqLoci(loci: RowLoci, c: Conversion): SeqLoci {
  const pairs = convert(loci, RowLocation.create(loci.set), SeqLocation.create(loci.set), rowToSeq, 'unit', 'index', c);
  markAmbiguous(
    c,
    loci.set,
    loci.alignment,
    Array.from(loci.elements, (e) => e.row as number),
  );
  return grouped({ kind: 'seq-loci', set: loci.set, elements: [] }, pairs, false, loci.set);
}

/** Unit positions to the residues of an alignment's rows bound to them. */
function seqToRowLoci(loci: SeqLoci, alignment: AlignmentId, c: Conversion): RowLoci {
  const one = (l: SeqLocation, out: RowLocation, cb: Visit<RowLocation>, acc: Acc) =>
    toRowLocations(l, alignment, out, cb, acc);
  const pairs = convert(loci, SeqLocation.create(loci.set), RowLocation.create(loci.set), one, 'row', 'index', c);
  markAmbiguous(c, loci.set, alignment, pairs.keys());
  return grouped({ kind: 'row-loci', set: loci.set, alignment, elements: [] }, pairs, false, loci.set);
}

/** Unit positions along a composed path to a unit (a linked section); a residue's three nucleotides are `multiplied`. */
function seqAlongPath(loci: SeqLoci, path: ComposedPath, c: Conversion): SeqLoci {
  const one = (l: SeqLocation, out: SeqLocation, cb: Visit<SeqLocation>, acc: Acc) =>
    toSystemEach(l, path, out, cb, acc);
  const pairs = convert(loci, SeqLocation.create(loci.set), SeqLocation.create(loci.set), one, 'unit', 'index', c);
  return grouped({ kind: 'seq-loci', set: loci.set, elements: [] }, pairs, false, loci.set);
}

export const SeqLoci = { ...lociKind<SeqLoci, SeqLocation>(), toRowLoci: seqToRowLoci, toSystem: seqAlongPath };
export const BoundaryLoci = lociKind<BoundaryLoci, BoundaryLocation>();
export const RowLoci = {
  ...lociKind<RowLoci, RowLocation>(),
  toSeqLoci: rowToSeqLoci,
  toAlignedLoci: rowToAlignedLoci,
};
export const AlignedLoci = {
  ...lociKind<AlignedLoci, AlignedLocation>(),
  toRowLoci: alignedToRowLoci,
  toColumnLoci: alignedToColumnLoci,
};
export const ColumnLoci = { ...lociKind<ColumnLoci, ColumnLocation>(), toAlignedLoci: columnToAlignedLoci };
export const ItemLoci = lociKind<ItemLoci, ItemLocation>();

/**
 * The annotation's items at the loci's locations: unit positions (as `itemsAt`), row residues, cells (the
 * row's residue and the column) and columns, or the annotation's part of item loci. Edges.
 */
function toItemLoci(loci: Loci, annotation: AnnotationId, c: Conversion): ItemLoci {
  const { set } = loci;
  const ann = set.annotations[annotation];
  const found: number[] = [];
  let inputs = 0;
  c.dropped = 0;
  c.multiplied = 0;
  const out = ItemLocation.create(set);
  const push = (l: ItemLocation) => found.push(l.item);
  const count = (n: number) => {
    inputs++;
    if (n === 0) c.dropped++;
    else if (n > 1) c.multiplied++;
  };
  const space = { kind: CoordinateSystemKind.columns, id: -1, row: -1 };
  const over = (kind: number, id: number, row: number, index: number) => {
    Object.assign(space, { kind, id, row });
    const slot = Annotations.slotOf(ann as Annotation, space);
    if (slot >= 0 && index >= 0) Annotations.itemsOver(ann as Annotation, slot, index, index + 1, visitRow, items);
  };
  if (ann && loci.kind === 'item-loci') {
    for (const e of loci.elements)
      if (e.annotation === annotation) OrderedSet.forEach(e.items, (i) => found.push(i), null);
    inputs = found.length;
  } else if (ann && loci.kind === 'seq-loci') {
    forEachLocation(loci, SeqLocation.create(set), (l) => count(itemsAt(l, annotation, out, push, null)), null);
  } else if (ann && (loci.kind === 'row-loci' || loci.kind === 'aligned-loci' || loci.kind === 'column-loci')) {
    const at = loci.kind === 'row-loci' ? RowLocation : loci.kind === 'aligned-loci' ? AlignedLocation : ColumnLocation;
    const rows = set.alignments[loci.alignment]?.rows;
    forEachLocation(
      loci,
      at.create(set) as Location,
      (l) => {
        const q = beginQuery(ann, annotation, set, out, push, null);
        const o = l as unknown as Rec;
        if (l.kind === 'row') over(CoordinateSystemKind.row, l.alignment, l.row, l.index);
        else over(CoordinateSystemKind.columns, o.alignment, -1, o.column);
        if (l.kind === 'aligned')
          over(CoordinateSystemKind.row, l.alignment, l.row, rows?.colToPos(l.row, l.column) ?? -1);
        count(q.n);
      },
      null,
    );
  } else forEachLocation(loci, {} as Location, () => count(0), null);
  c.status =
    inputs > 0 && found.length === 0
      ? ConversionStatus.unmapped
      : c.dropped
        ? ConversionStatus.partial
        : ConversionStatus.exact;
  const elements = found.length
    ? [{ annotation, items: OrderedSet.ofSortedArray(SortedArray.ofUnsortedArray(Int32Array.from(found))) }]
    : [];
  return { kind: 'item-loci', set, elements: elements as ItemLoci['elements'] };
}

export const Loci = { fromLocation, has, forEachLocation, union, isEmpty, areEqual, remap, toItemLoci };
