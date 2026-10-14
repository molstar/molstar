/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Hit testing: from (x, y) on a section's surface to its region, entry, block, shown and model
// index, and the locations for the hit's kind, written into a reused `Hit`. Headless and allocation-free: entries
// by a scan of the few pinned ones or a binary search, blocks by a binary search over offsets, the shown index to the
// block's through its layout map (a reversed or cyclic piece alike), items through the annotation's interval index, so
// hits are logarithmic in features. `Hit`, its regions and kinds live here (hover reads them).

import { Lanes } from '@molstar/seq-core/format';
import {
  AlignedLocation,
  Annotations,
  ColumnLocation,
  type ComposedPath,
  type CoordinateSystem,
  CoordinateSystemKind,
  Coords,
  ItemLocation,
  type MapHit,
  RowLocation,
  SeqLocation,
  type SeqSet,
} from '@molstar/seq-core/model';
import { AxisGeometry, type SectionAxisView, type ShownIndex } from './axis.js';
import { placedRingLabels } from './draw/ring-labels.js';
import { hitLayer } from './draw/sizes.js';
import { hitStreamed } from './draw/streamed.js';
import {
  EntryKind,
  entryBorder,
  facingOf,
  hasScrollbar,
  type SectionEntry,
  type SectionLayout,
  SectionLayouts,
  type SplitFeatures,
} from './layout.js';
import type { SeqViewerRegistry } from './registry.js';
import { hitItem } from './renderer-helpers.js';
import { Ring } from './ring.js';
import { linePositionAt, linePositionDrawnNear, runXInto } from './runs.js';
import type { StreamView } from './streams.js';

/** A layer's regions by its ref (`ctx.streams.get`): a hit reads a streamed layer's items through `v.stream`. */
export type StreamLookup = (layer: string) => StreamView | null;

export const HitRegion = {
  none: 0,
  data: 1,
  column: 2,
  header: 3,
  gutter: 4,
  dataColumn: 5,
  dataColumnHeader: 6,
  navigator: 7,
  vnav: 8,
  scrollbar: 9,
} as const;
export const HoverKind = { none: 0, seq: 1, row: 2, aligned: 3, column: 4, item: 5, header: 6 } as const;

/** What is under the pointer, reused; -1: none. S4's hit testing fills it. */
export interface Hit {
  /** The point tested, surface px. */
  x: number;
  y: number;
  /**
   * The same point as doubles, for renderers' `hit` and the header's: read from it in a loop, a field of this record
   * boxed each read (16 bytes a move, measured on the tree's hit).
   */
  readonly point: Float64Array;
  section: number;
  /** `HitRegion`. */
  region: number;
  entry: number;
  track: number;
  layer: number;
  lane: number;
  displayRow: number;
  dataRow: number;
  block: number;
  shown: number;
  /** The model index on the block's axis. */
  index: number;
  annotation: number;
  item: number;
  side: number;
  treeNode: number;
  dataColumn: number;
  /** A matrix's state (its band) under the point: a fingerprint's ligand; -1: none. */
  state: number;
  /** A split track's row: its value's, -2 its label's entry; -1: none. */
  split: number;
  /** A group's header entry: the group; -1: none (its entry is in the key). */
  group: number;
  /** A renderer's part of the element under the point: a fingerprint's kind (its stripe); -1: none. */
  part: number;
  /** A dot plot's partner, the position under y (seqspec:interaction#hover): its block, unit and index; -1: none. */
  pairBlock: number;
  pairUnit: number;
  pairIndex: number;
  /** A streamed item's [start, stop) on the hit's unit (its `part` its number in the layer; ST4b); -1, -1: none. */
  readonly span: Int32Array;
  /** `HoverKind`. */
  kind: number;
  /**
   * The element's coordinate system (seqspec:interaction#channels: the element): the unit's, a row's off a gap, the
   * columns', an item's.
   */
  system: CoordinateSystem;
  seq: SeqLocation;
  row: RowLocation;
  aligned: AlignedLocation;
  column: ColumnLocation;
  itemLoc: ItemLocation;
}
type NumberKeys<T> = { [K in keyof T]: T[K] extends number ? K : never }[keyof T];
/**
 * What makes a hover "the same element" (seqspec:interaction#input-events): the hit's numbers, never a move inside one;
 * its `system` (kind, id, row) is the key's last three slots. One list for the dispatcher and the hover manager,
 * type-checked against `Hit` and compared by `HitKey` (no loop over the names: a keyed load per field is megamorphic).
 */
export const HIT_ELEMENT = [
  'section',
  'region',
  'kind',
  'entry',
  'layer',
  'displayRow',
  'block',
  'shown',
  'annotation',
  'item',
  'side',
  'treeNode',
  'dataColumn',
  'state',
  'split',
  'part',
  'pairBlock',
  'pairUnit',
  'pairIndex',
] as const satisfies readonly NumberKeys<Hit>[];
/** `HitKey`'s record: the `HIT_ELEMENT` names in order, then the system's kind, id and row. */
const SIZE = HIT_ELEMENT.length + 3;
export const HitKey = {
  SIZE,
  /** Straight-line on purpose: `_test/hit-key.test.ts` changes every field of the list and of `system` in turn. */
  write(h: Hit, k: Int32Array) {
    k[0] = h.section;
    k[1] = h.region;
    k[2] = h.kind;
    k[3] = h.entry;
    k[4] = h.layer;
    k[5] = h.displayRow;
    k[6] = h.block;
    k[7] = h.shown;
    k[8] = h.annotation;
    k[9] = h.item;
    k[10] = h.side;
    k[11] = h.treeNode;
    k[12] = h.dataColumn;
    k[13] = h.state;
    k[14] = h.split;
    k[15] = h.part;
    k[16] = h.pairBlock;
    k[17] = h.pairUnit;
    k[18] = h.pairIndex;
    k[19] = h.system.kind;
    k[20] = h.system.id;
    k[21] = h.system.row;
  },
  same(h: Hit, k: Int32Array): boolean {
    return (
      k[0] === h.section &&
      k[1] === h.region &&
      k[2] === h.kind &&
      k[3] === h.entry &&
      k[4] === h.layer &&
      k[5] === h.displayRow &&
      k[6] === h.block &&
      k[7] === h.shown &&
      k[8] === h.annotation &&
      k[9] === h.item &&
      k[10] === h.side &&
      k[11] === h.treeNode &&
      k[12] === h.dataColumn &&
      k[13] === h.state &&
      k[14] === h.split &&
      k[15] === h.part &&
      k[16] === h.pairBlock &&
      k[17] === h.pairUnit &&
      k[18] === h.pairIndex &&
      k[19] === h.system.kind &&
      k[20] === h.system.id &&
      k[21] === h.system.row
    );
  },
};
const NO_SET = null as unknown as SeqSet;
const ids = {
  track: -1,
  layer: -1,
  lane: -1,
  treeNode: -1,
  dataColumn: -1,
  state: -1,
  split: -1,
  group: -1,
  part: -1,
  pairBlock: -1,
  pairUnit: -1,
  pairIndex: -1,
};
export const Hit = {
  /** Everything but the point and the section: no element. */
  reset: (out: Hit) => reset(out),
  create: (): Hit => ({
    ...ids,
    x: 0.5,
    y: 0.5,
    point: new Float64Array(2),
    span: new Int32Array([-1, -1]),
    section: -1,
    region: HitRegion.none,
    entry: -1,
    displayRow: -1,
    dataRow: -1,
    block: -1,
    shown: -1,
    index: -1,
    annotation: -1,
    item: -1,
    side: -1,
    kind: HoverKind.none,
    system: { kind: -1, id: -1, row: -1 },
    seq: SeqLocation.create(NO_SET),
    row: RowLocation.create(NO_SET),
    aligned: AlignedLocation.create(NO_SET),
    column: ColumnLocation.create(NO_SET),
    itemLoc: ItemLocation.create(NO_SET),
  }),
};

/** A block index's data index through a placement's path (to → from; a source used twice: its first). -1: none. */
export function sourceAt(p: ComposedPath, i: number): number {
  if (!Coords.sameSystem(p.map.from, p.from)) return Coords.toAt(p.map, i, 0);
  mapHit.first = -1;
  Coords.forEachFromAt(p.map, i, mapHit, takeFirst, mapHit);
  return mapHit.first;
}
/**
 * Block index i's copy of its source residue through `p`: its rank, in block index order, among the block indices that
 * map to that residue; 1 for the first or only one, and for a reversed piece the same as forward.
 */
export function copyAt(p: ComposedPath | null | undefined, i: number): number {
  if (!p || Coords.sameSystem(p.map.from, p.from)) return 1;
  const j = Coords.toAt(p.map, i, 0);
  copies.at = i;
  copies.rank = 1;
  return j >= 0 && Coords.forEachFromAt(p.map, j, mapHit, below, copies) > 1 ? copies.rank : 1;
}
const copies = { at: 0, rank: 1 };
const below = (h: MapHit, c: typeof copies) => {
  if (h.element < c.at) c.rank++;
};
const mapHit = { element: -1, phase: 0, first: -1 };
const takeFirst = (h: MapHit, c: typeof mapHit) => {
  if (c.first < 0) c.first = h.element;
};
const system = (out: Hit, kind: number, id: number, row: number) => {
  out.system.kind = kind;
  out.system.id = id;
  out.system.row = row;
};
/** The point (x from the line's start), the line's left and the DPR, for the wrapped position search. */
const point = new Float64Array(3);
/** How far past a block's last cell (CSS px) its layers are still asked: the largest glyph's radius. */
const END_REACH = 12;

/** Px at the right of a `box` section's rows that hold its drawn scrollbar when the rows overflow. */
export const SCROLLBAR = 8;

/** The best row so far, the rows the filter keeps, and on a split row its value's rows (`split[row] === want`). */
const job = { best: -1, kept: null as Uint8Array | null, split: null as Int32Array | null, want: -1 };
type Job = typeof job;
const visit = (row: number, j: Job) => {
  if ((j.kept === null || j.kept[row]) && (j.split === null || j.split[row] === j.want) && row > j.best) j.best = row;
};
const noop = (_l: SeqLocation, _c: null) => {};

function reset(out: Hit) {
  out.region = HitRegion.none;
  out.kind = HoverKind.none;
  out.entry = out.track = out.layer = out.lane = out.displayRow = out.dataRow = -1;
  out.block = out.shown = out.index = out.annotation = out.item = out.side = out.treeNode = out.dataColumn = -1;
  out.state = out.split = out.group = out.part = out.span[0] = out.span[1] = -1;
  out.pairBlock = out.pairUnit = out.pairIndex = -1;
  system(out, -1, -1, -1);
}

/**
 * Fills `out` for its point (`out.x`, `out.y`, surface px); false when outside the surface (region none). The point
 * is read from the record: a fractional double passed to a call V8 does not inline is boxed.
 */
export function hitTest(
  e: SectionEntry,
  v: SectionAxisView,
  registry: SeqViewerRegistry,
  out: Hit,
  streams: StreamLookup | null = null,
  ownResidue = false,
): boolean {
  const l = e.layout;
  if (l.geometry !== AxisGeometry.circular || !(out.x >= 0 && out.y >= 0 && out.x < l.width && out.y < l.height))
    return hitAt(e, v, registry, out, streams, ownResidue);
  // Polar (seqspec:viewer#the-render-backend): the point into the unrolled axis the rings are drawn from, hit there as
  // a linear one.
  const x = out.x;
  const y = out.y;
  ring.from(l, v, 1).toUnrolled(x, y, polar);
  // Outside the outermost ring: an item's label (CG3) is hit as its item, at a point of its track; else nothing.
  const labels = polar[1] < l.scrolledTop ? placedRingLabels.get(l) : undefined;
  const k = labels && labels.at === e.live.version.value ? labels.labelAt(out) : -1;
  out.x = polar[0];
  if (k >= 0 && labels) {
    out.x = labels.ux[3 * k];
    polar[1] = labels.ux[3 * k + 1];
  }
  out.y = polar[1] >= l.scrolledTop ? polar[1] : -1;
  v.geometry = AxisGeometry.linear;
  try {
    const hit = hitAt(e, v, registry, out, streams, ownResidue);
    if (k < 0 || !labels || !e.model || out.region !== HitRegion.data) return hit;
    out.layer = labels.item[3 * k];
    return hitItem(out, e.model.set, labels.item[3 * k + 1], labels.item[3 * k + 2]);
  } finally {
    v.geometry = AxisGeometry.circular;
    out.x = x;
    out.y = y;
  }
}
const ring = new Ring();
const polar = new Float64Array(2);

function hitAt(
  e: SectionEntry,
  v: SectionAxisView,
  registry: SeqViewerRegistry,
  out: Hit,
  streams: StreamLookup | null,
  ownResidue: boolean,
): boolean {
  const x = out.x;
  const y = out.y;
  out.point[0] = x;
  out.point[1] = y;
  reset(out);
  out.section = e.index;
  const l = e.layout;
  const ax = e.axis;
  const model = e.model;
  // A section with no blocks (a source failed) hits nothing.
  const inside = l.geometry === AxisGeometry.circular ? y >= 0 : x >= 0 && y >= 0 && x < l.width && y < l.height;
  if (!ax || !model || ax.count === 0 || !inside) return false;
  if (y < l.navHeight) return region(out, HitRegion.navigator);
  // The rulerRail: the pinned ruler's row (its lower part on columns), over the data area, is the navigator.
  const re = l.railEntry;
  if (re >= 0 && x >= l.dataLeft && x < l.dataRight) {
    if (y >= l.railTop && y < l.entryTop[re] + l.entryHeight[re]) return region(out, HitRegion.navigator);
  }
  if (l.vnavWidth > 0 && x >= l.width - l.vnavWidth) return region(out, HitRegion.vnav);
  // The entry: a scan of the pinned ones, else a binary search in the scrolled content.
  let entry = -1;
  let yc = y;
  if (y < l.scrolledTop) {
    for (let k = 0; k < l.pinnedCount && entry < 0; k++)
      if (y >= l.entryTop[k] && y < l.entryTop[k] + l.entryHeight[k]) entry = k;
  } else {
    // `SectionLayouts.entryAt`, inline: a fractional double passed to a call V8 does not inline is boxed.
    yc = y - l.scrolledTop + v.scrollY;
    let lo = l.pinnedCount;
    for (let hi = l.entryCount - 1; lo < hi; ) {
      const mid = (lo + hi + 1) >> 1;
      if (l.entryTop[mid] <= yc) lo = mid;
      else hi = mid - 1;
    }
    if (lo < l.entryCount && yc >= l.entryTop[lo] && yc < l.entryTop[lo] + l.entryHeight[lo]) entry = lo;
  }
  // The region from x against the layout's edges.
  const p = e.props;
  const nc = l.columnCount;
  const columnsEnd = nc > 0 ? l.columnX[nc - 1] + l.columnWidth[nc - 1] : l.headerWidth + p.controlWidth;
  if (x < l.headerWidth) out.region = HitRegion.header;
  else if (x < columnsEnd && nc > 0) {
    out.region = y < l.scrolledTop ? HitRegion.dataColumnHeader : HitRegion.dataColumn;
    for (let k = 0; k < nc; k++) if (x >= l.columnX[k] && x < l.columnX[k] + l.columnWidth[k]) out.dataColumn = k;
  } else if (x < l.dataLeft) out.region = HitRegion.gutter;
  else if (x >= l.dataRight) {
    const over = hasScrollbar(p) && l.vscroll;
    out.region = over && x >= l.width - l.vnavWidth - SCROLLBAR ? HitRegion.scrollbar : HitRegion.gutter;
  } else out.region = HitRegion.data;
  if (entry < 0) return true;
  out.entry = entry;
  // A group's header: the group, whose hover lights its tracks (header-hover.ts).
  if (l.entryKind[entry] === EntryKind.group && out.region === HitRegion.header) {
    out.kind = HoverKind.header;
    out.group = l.entryGroup[entry];
    return true;
  }
  if (l.entryKind[entry] !== EntryKind.track) return true;
  out.track = l.entryTrack[entry];
  out.layer = l.entryLayer[entry];
  if (l.entryRow[entry] === -2 || SectionLayouts.splitRow(l, entry) >= 0) out.split = l.entryRow[entry];
  const rv = model.resolved;
  const set = model.set;
  const rows = l.rows;
  const onRows =
    rows.layer >= 0 && l.geometry === AxisGeometry.linear && out.track === rows.track && entry >= rows.entry;
  if (onRows) {
    const d = rows.lazy ? Math.floor((yc - l.entryTop[entry]) / rows.height) : l.entryRow[entry];
    out.displayRow = d;
    out.lane = SectionLayouts.laneOf(rows, d);
    out.dataRow = out.lane >= 0 ? Lanes.rowOf(rv.lanes, out.lane) : -1;
    if (rows.node && rows.node[d] >= 0) out.treeNode = rows.node[d];
  }
  if (out.region === HitRegion.header) {
    out.kind = HoverKind.header;
    const r = onRows ? registry.get('renderers', set.strings.get(rv.layers.type[rows.layer])) : undefined;
    const h = r?.header?.hit;
    if (h) {
      // The header's view: display row 0's y and a row's height (the tree's rows; restored for the passes).
      const top = v.top;
      const height = v.height;
      v.top = l.scrolledTop + l.entryTop[rows.entry] - v.scrollY;
      v.height = rows.height;
      h(v, model.prepared[rows.layer]?.state, rv.layers.params[rows.layer], out);
      v.top = top;
      v.height = height;
    }
    return true;
  }
  if (out.region !== HitRegion.data) return true;

  // The block and the shown index: wrapped by the line, else by a binary search over offsets.
  let b: number;
  let shown: number;
  // Below a device px a cell, the first of the pixel's cells (`drawnFirst`): an item drawn in its column is hit.
  let shownLo = -1;
  const wrapped0 = l.geometry === AxisGeometry.wrapped;
  if (wrapped0) {
    const line = l.entryRow[entry];
    b = l.lines.block[line];
    const run = e.runs?.[b];
    point[0] = x - l.dataLeft;
    point[1] = l.dataLeft;
    point[2] = v.dpr;
    shown = run ? linePositionAt(run, l.lines.start[line], l.lines.stop[line], p.wrappedGroupGap, point) : -1;
  } else {
    // A binary search over block starts, inline: a double passed to a call V8 does not inline is boxed.
    b = ax.stacked ? l.entryBlock[entry] : 0;
    for (let hi = ax.stacked ? 0 : ax.count - 1; b < hi; ) {
      const mid = (b + hi + 1) >> 1;
      if (v.blockX[mid] <= x) b = mid;
      else hi = mid - 1;
    }
    const pd = Math.floor(x * v.dpr) | 0;
    // By the drawn edges (GG-01): the next block owns the px its first edge rounds onto, and a block does not own the
    // px before its first edge (LC-26: the search compares unsnapped starts).
    if (!ax.stacked && b + 1 < ax.count && pd >= Math.round(v.blockX0[b + 1] * v.dpr) - v.scrollDev) b++;
    else if (!ax.stacked && b > 0 && pd < Math.round(v.blockX0[b] * v.dpr) - v.scrollDev) b--;
    // The cell drawn on the pointer's device px, not the one its unsnapped x falls in.
    shown = v.drawnShown(b, Math.floor((x - v.blockX[b]) / (ax.span[b] * v.cell)), pd);
    shownLo = v.drawnFirst(b, shown, pd);
  }
  // In the gap past a linear block's last cell or before the next block's first, within an end glyph's reach (the
  // nearer edge where the gap is narrower than two): that block with no position, for the layers whose glyphs overhang
  // a block's edge (`hitsPastEnd`: a boundary's) to hit by drawn geometry, on either side (GG-02).
  let past = false;
  // Before the block's first edge (else past its last).
  let first = false;
  if (!wrapped0 && b >= 0 && shown >= ax.shown[b]) {
    const after = x - v.blockX[b] - ax.shown[b] * ax.span[b] * v.cell;
    const before = !ax.stacked && b + 1 < ax.count ? v.blockX[b + 1] - x : END_REACH;
    if (before < END_REACH && before < after && ax.shown[b + 1] > 0) {
      b++;
      shown = -1;
      past = first = true;
    } else past = ax.shown[b] > 0 && after < END_REACH;
  }
  if (!past && (b < 0 || shown < 0 || shown >= ax.shown[b])) return true;
  out.block = b;
  out.shown = past ? -1 : shown;
  out.index = past ? -1 : ax.indexOf(b, shown as ShownIndex);

  // Locations for the hit's kind; just past a block's end, none: only a layer's own hit (an end glyph's) counts.
  if (!past) {
    const a = ax.alignment[b];
    if (a >= 0 && !onRows) {
      out.region = entry < l.pinnedCount ? HitRegion.column : out.region;
      out.kind = HoverKind.column;
      const c = out.column;
      c.set = set;
      c.alignment = a as never;
      c.column = out.index as never;
      system(out, CoordinateSystemKind.columns, a, -1);
    } else if (a >= 0) {
      if (out.dataRow < 0) return true;
      out.kind = HoverKind.aligned;
      out.row.row = out.seq.unit = -1 as never;
      const al = out.aligned;
      al.set = set;
      al.alignment = a as never;
      al.row = out.dataRow as never;
      al.column = out.index as never;
      // Off a gap: the row's residue, and the bound unit's through the row's binding (else their row and unit stay -1).
      const residue = AlignedLocation.toRow(al, out.row);
      if (residue) RowLocation.toSeqLocations(out.row, out.seq, noop, null);
      system(out, residue ? CoordinateSystemKind.row : CoordinateSystemKind.columns, a, residue ? out.dataRow : -1);
    } else {
      out.kind = HoverKind.seq;
      const s = out.seq;
      s.set = set;
      s.unit = ax.unit[b] as never;
      s.index = out.index as never;
      system(out, CoordinateSystemKind.unit, ax.unit[b], -1);
      if (out.split >= 0) ownUnitAt(e, b, out);
    }
  }

  // Items: the topmost layer of the track first; a renderer's own hit (on the entry's band), else ranges through the
  // interval index.
  const t = out.track;
  let last = out.layer;
  while (last + 1 < rv.layers.count && rv.layers.track[last + 1] === t) last++;
  const top = v.top;
  const height = v.height;
  const only = v.onlyBlock;
  const row = v.row;
  // Wrapped: the line's block and shown range, as its draw has them (a renderer hits its underlines on the line).
  const wrapped = l.geometry === WRAPPED;
  const line = v.line;
  const lineStart = v.lineStart;
  const lineStop = v.lineStop;
  if (wrapped) {
    v.line = l.entryRow[entry];
    v.lineStart = l.lines.start[v.line];
    v.lineStop = l.lines.stop[v.line];
  }
  // On a rows track's row its lane offset (-2: a collapsed clade's wedge): its features are found by lane.
  v.row = onRows ? (out.lane >= 0 ? out.lane - rows.laneStart : -2) : -1;
  v.top = entry < l.pinnedCount ? l.entryTop[entry] : l.scrolledTop + l.entryTop[entry] - v.scrollY;
  // The row's content stands above its border, as the composer draws it (`SectionLayouts.entryBorder`; the border's
  // device px inline: a double passed to a call is boxed).
  v.height = l.entryHeight[entry];
  if (entryBorder(l, entry, p.rowBorders)) v.height -= (Math.max(1, Math.round(v.dpr)) | 0) / v.dpr;
  v.onlyBlock = ax.stacked || wrapped ? l.entryBlock[entry] : -1;
  let found = false;
  // A split row: its layer's row alone; an expanded split track's label entry (row -2): not the split layer.
  const at = SectionLayouts.splitRow(l, entry);
  for (let k = last; k >= 0 && rv.layers.track[k] === t && !found; k--) {
    if ((at >= 0 && k !== out.layer) || (l.entryRow[entry] === -2 && k === rv.tracks.split[t])) continue;
    const type = rv.layers.type[k];
    const r = type >= 0 ? registry.get('renderers', set.strings.get(type)) : undefined;
    if (past && r?.hitsPastEnd !== true) continue;
    const s = model.prepared[k]?.state;
    // Its themes, so a renderer's `hit` sizes its glyphs and bands as its draw does (hits by drawn geometry).
    hitLayer.axis = ax;
    hitLayer.layer = k;
    hitLayer.themes = model.prepared[k]?.themes ?? NO_THEMES;
    v.facing = facingOf(l, ax, k, t, entry);
    // A streamed layer's items are its regions' (ST3b: the hover's `part` names one).
    v.stream = streams ? streams(set.strings.get(rv.layers.ref[k])) : null;
    found =
      r?.hit?.(v, at >= 0 ? ((s as SplitFeatures).split?.rows[at] ?? s) : s, rv.layers.params[k], out) === true ||
      (!past && itemAt(e, k, b, shownLo < 0 ? shown : shownLo, shown, onRows ? out.lane : -1, out, at)) ||
      (!past && hitStreamed(v, out, false));
    if (found) out.layer = k;
  }
  v.top = top;
  v.height = height;
  v.onlyBlock = only;
  v.row = row;
  v.line = line;
  v.lineStart = lineStart;
  v.lineStop = lineStop;
  v.stream = null;
  if (past && !found) {
    out.block = -1;
    return true;
  }
  // An end glyph's item, at the block's first position before it (LC-28), else at its last.
  if (past) {
    out.shown = first ? 0 : ax.shown[b] - 1;
    out.index = ax.indexOf(b, out.shown as ShownIndex);
  }
  // `UI.OwnRowHover` "residue": a bar on an own unit's row hovers that unit's residue, the item kept for its card.
  if (ownResidue && !past && found && at >= 0 && out.kind === HoverKind.item && ax.alignment[b] < 0) {
    const s = out.seq;
    s.set = set;
    s.unit = ax.unit[b] as never;
    s.index = out.index as never;
    if (ownUnitAt(e, b, out)) out.kind = HoverKind.seq;
  }
  return true;
}

/**
 * Wrapped: the line and the shown position nearest a point (`at[0]`, `at[1]`, surface px) as a text selection takes
 * them (seqspec:interaction#selection): the last line starting at or above y (above the first line the first, below
 * the last the last), then on it the last position drawn starting at or before x (`linePositionDrawnNear`), by x on
 * every line.
 * Into `out`: line, block, shown, and 1 when the point is past that position's cell (a gap, past the line's end, an
 * inline label), else 0. False: not wrapped, or no lines. Allocation-free.
 */
export function wrappedPositionNear(l: SectionLayout, v: SectionAxisView, at: Float64Array, out: Int32Array): boolean {
  const t = l.lines;
  if (l.geometry !== WRAPPED || t.count === 0) return false;
  const yc = at[1] - l.scrolledTop + v.scrollY;
  let lo = 0;
  for (let hi = t.count - 1; lo < hi; ) {
    const mid = (lo + hi + 1) >> 1;
    if (t.top[mid] <= yc) lo = mid;
    else hi = mid - 1;
  }
  const run = v.runs[t.block[lo]];
  if (!run) return false;
  out[0] = lo;
  out[1] = t.block[lo];
  near[0] = at[0] - l.dataLeft;
  near[1] = l.dataLeft;
  near[2] = v.dpr;
  // By the drawn cells, as a hover (LC-27): the last position whose snapped left edge is at or before the pointer's
  // device px, and past it where the px is past its snapped right edge.
  out[2] = linePositionDrawnNear(run, t.start[lo], t.stop[lo], v.gap, near);
  runXInto(run, t.start[lo], out[2], v.gap, edge, 0);
  out[3] = Math.floor(at[0] * v.dpr) >= Math.round((l.dataLeft + edge[0] + run.advance) * v.dpr) ? 1 : 0;
  return true;
}
const near = new Float64Array(3);
const edge = new Float64Array(1);
// Read once: a member read per pointer event allocated under the test runner's module wrapper.
const WRAPPED = AxisGeometry.wrapped;

const UNIT = CoordinateSystemKind.unit;
const NO_THEMES: readonly never[] = [];
/**
 * On a split row whose data sits on one unit (`SplitRows.unit`), the residue of that unit there, through its data's
 * placement on the block: the hover is that unit's; off its residues, the block's position as it is.
 */
function ownUnitAt(e: SectionEntry, b: number, out: Hit): boolean {
  const { resolved: rv, set, prepared } = e.model as NonNullable<SectionEntry['model']>;
  const k = rv.tracks.split[out.track];
  const unit = k >= 0 ? (prepared[k]?.state as SplitFeatures | undefined)?.split?.unit : undefined;
  const u = unit !== undefined && out.split < unit.length ? unit[out.split] : -1;
  const pl = rv.placements;
  const p0 = rv.layers.placementStart[k];
  for (let p = p0; u >= 0 && u !== out.seq.unit && p < p0 + rv.layers.placementCount[k]; p++) {
    if (pl.block[p] !== b + rv.sections.blockStart[e.index]) continue;
    if (set.annotations[pl.annotation[p]]?.slots.unit[pl.slot[p]] !== u) continue;
    const j = pl.path[p] >= 0 ? sourceAt(rv.paths[pl.path[p]], out.index) : out.index;
    if (j < 0) continue;
    out.seq.unit = u as never;
    out.seq.index = j as never;
    system(out, UNIT, u, -1);
    return true;
  }
  return false;
}

function region(out: Hit, r: number): boolean {
  out.region = r;
  return true;
}

const pixelRuns = new Int32Array(64);
/**
 * The last row (drawn on top) of a layer's ranges covering block b's shown [k0, k1] (the cells counted on the pointer's
 * device px: one, or several below a device px a cell, as a density column holds them); on rows, its lane's or the
 * columns'.
 */
function itemAt(e: SectionEntry, k: number, b: number, k0: number, k1: number, lane: number, out: Hit, at: number) {
  const { resolved: rv, set, prepared } = e.model as NonNullable<SectionEntry['model']>;
  const ax = e.axis;
  if (!ax) return false;
  const block = b + rv.sections.blockStart[e.index];
  const index = ax.indexOf(b, k1 as ShownIndex);
  // The pixel's cells as block-system runs (one run where it holds one cell).
  const runs = k0 < k1 ? Math.min(pixelRuns.length >> 1, ax.shownToRuns(b, k0, k1 + 1, pixelRuns)) : 0;
  const split = at >= 0 ? ((prepared[k]?.state as SplitFeatures | undefined)?.split ?? null) : null;
  const pl = rv.placements;
  const p0 = rv.layers.placementStart[k];
  for (let p = p0; p < p0 + rv.layers.placementCount[k]; p++) {
    if (pl.block[p] !== block || (pl.lane[p] >= 0 && pl.lane[p] !== lane)) continue;
    const ann = set.annotations[pl.annotation[p]];
    if (ann?.geometry !== 'ranges') continue;
    const path = pl.path[p] >= 0 ? rv.paths[pl.path[p]] : null;
    const j = path ? sourceAt(path, index) : index;
    if (j < 0 && runs === 0) continue;
    job.best = -1;
    job.kept = job.split = null;
    job.want = at;
    for (let d = rv.layers.dataStart[k]; d < rv.layers.dataStart[k] + rv.layers.dataCount[k]; d++)
      if (rv.data[d] === ann.id) {
        job.kept = rv.kept[d];
        job.split = split ? split.rowOf[d - rv.layers.dataStart[k]] : null;
      }
    if (j >= 0) Annotations.itemsOver(ann, pl.slot[p], j, j + 1, visit, job);
    for (let r = 0; r < runs; r++) {
      // Through a path, a run's mapped ends' hull (a run's ends off the path: its own index alone, above).
      const a = path ? sourceAt(path, pixelRuns[2 * r]) : pixelRuns[2 * r];
      const z = path ? sourceAt(path, pixelRuns[2 * r + 1] - 1) : pixelRuns[2 * r + 1] - 1;
      if (a >= 0 && z >= 0) Annotations.itemsOver(ann, pl.slot[p], Math.min(a, z), Math.max(a, z) + 1, visit, job);
    }
    if (job.best < 0) continue;
    out.kind = HoverKind.item;
    out.annotation = ann.id;
    out.item = ann.item[job.best];
    const loc = out.itemLoc;
    loc.set = set;
    loc.annotation = ann.id as never;
    loc.item = ann.item[job.best] as never;
    system(out, CoordinateSystemKind.item, ann.id, ann.item[job.best]);
    return true;
  }
  return false;
}
