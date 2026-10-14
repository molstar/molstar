/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// What a header's hover stands for (seqspec:interaction#hover), as the `headerHover` behavior describes it: a track's
// header the items its feature layers show, a split row's header that row's, a rows track's row header the items placed
// on its row, a group's header its shown tracks'. Continuous layers (letters, values, matrices) light nothing positional:
// the row band only. The description (`HeaderSet`, reused in the hover) is counted when the hover lands, allocation-free:
// its size (items, covered positions) is known before anything is resolved. Its runs for the overlay are built then only
// when it lights (at most the behavior's `maxItems`); its item loci (what `hover.toLoci()` gives a host or a 3D link)
// only when read.

import {
  type Annotation,
  Extents,
  type ItemLoci,
  CoordinateSystemKind as K,
  type Loci,
  OrderedSet,
  type SeqSet,
} from '@molstar/seq-core/model';
import type { ViewModel } from './graph.js';
import type { SplitFeatures } from './layout.js';

/** What a layer's header lights: its shown items, or the row band alone (which the hover draws for any header). */
export type HeaderLight = 'items' | 'rowBand';
/** The `headerHover` behavior's params (seqspec:viewer#behaviors). */
export interface HeaderHoverParams {
  /** A feature layer's (over item annotations: ranges, pairs): its shown items, or the row band alone. */
  readonly features: HeaderLight;
  /** A continuous layer's (letters, values dense or sparse, a matrix): the row band alone. */
  readonly continuous: 'rowBand';
  /** Per layer type (its registry name), over the two above. */
  readonly types: Readonly<Record<string, HeaderLight>>;
  /** More shown items than this: the row band only, and nothing forwarded. */
  readonly maxItems: number;
  /** Whether linked views get the set: `hover.toLoci()` (a 3D link's reading) gives it, else null. */
  readonly forward: boolean;
}
export const DefaultHeaderHoverParams: HeaderHoverParams = {
  features: 'items',
  continuous: 'rowBand',
  types: {},
  maxItems: 5_000,
  forward: true,
};

/** A loci's size without its positions (seqspec:interaction#hover): items, and covered positions (an estimate). */
export interface LociSize {
  readonly items: number;
  readonly length: number;
}

/** A header hover's description, reused (`HoverState.header`). */
export interface HeaderSet {
  /** Described by the behavior; false: a header with nothing of its own (the behavior removed, or another hover). */
  on: boolean;
  /** The track (or -1) and its split row (≥ 0) or row lane (≥ 0), else the group (or -1): resolved indices. */
  track: number;
  split: number;
  lane: number;
  group: number;
  /** Shown items, and the positions their rows cover (overlaps counted twice: an estimate). */
  items: number;
  length: number;
  /** Lit on the overlay (its runs in `extents`), and given to linked views (`hover.toLoci()`). */
  light: boolean;
  forward: boolean;
  /** Its rows' runs when it lights (reused). */
  readonly extents: Extents;
}
export const headerSet = (): HeaderSet => ({
  on: false,
  track: -1,
  split: -1,
  lane: -1,
  group: -1,
  items: 0,
  length: 0,
  light: false,
  forward: false,
  extents: Extents.create(16),
});

/** The params a description was made with (`items` per layer reads them). */
let P: HeaderHoverParams = DefaultHeaderHoverParams;
/**
 * What a counting pass has seen, over every annotation of the set at once (an annotation may recur on several tracks
 * of a group: counted once): its items and its rows, at per-annotation offsets made once per set; the annotations a
 * pass touched are cleared at the next. The set is held weakly: a module keeps no viewer's set past its dispose (the
 * browser retention guard, canvas/retention.heap.ts).
 */
let seenSet: WeakRef<SeqSet> | null = null;
let itemAt = new Int32Array(0);
let rowAt = new Int32Array(0);
let seenItem = new Uint8Array(0);
let seenRow = new Uint8Array(0);
let touched = new Int32Array(16);
let touchedN = 0;
let touchedOf = new Uint8Array(0);
function seenFor(set: SeqSet) {
  if (set !== seenSet?.deref()) {
    const n = set.annotations.length;
    itemAt = new Int32Array(n + 1);
    rowAt = new Int32Array(n + 1);
    for (let a = 0; a < n; a++) {
      itemAt[a + 1] = itemAt[a] + (set.annotations[a]?.itemCount ?? 0);
      rowAt[a + 1] = rowAt[a] + (set.annotations[a]?.count ?? 0);
    }
    seenItem = new Uint8Array(itemAt[n]);
    seenRow = new Uint8Array(rowAt[n]);
    touchedOf = new Uint8Array(n);
    touchedN = 0;
    seenSet = new WeakRef(set);
    return;
  }
  for (let k = 0; k < touchedN; k++) {
    const a = touched[k];
    seenItem.fill(0, itemAt[a], itemAt[a + 1]);
    seenRow.fill(0, rowAt[a], rowAt[a + 1]);
    touchedOf[a] = 0;
  }
  touchedN = 0;
}
const ITEMS = ['ranges', 'pairs'];

/** Whether layer `i` lights its items under the current params. */
function lightsItems(model: ViewModel, i: number): boolean {
  const { resolved: rv, set } = model;
  const own = P.types[set.strings.get(rv.layers.type[i])];
  if (own) return own === 'items';
  for (let d = rv.layers.dataStart[i]; d < rv.layers.dataStart[i] + rv.layers.dataCount[i]; d++) {
    const ann = set.annotations[rv.data[d]];
    if (!ann || !ITEMS.includes(ann.geometry)) return false;
  }
  return rv.layers.dataCount[i] > 0 && P.features === 'items';
}

/** Whether row `r` of `ann` is placed by layer `i` (on a block, not left out of a row subset; on `lane` if ≥ 0). */
function placed(model: ViewModel, i: number, ann: Annotation, r: number, lane: number): boolean {
  const p = model.resolved.placements;
  const p0 = model.resolved.layers.placementStart[i];
  for (let k = p0; k < p0 + model.resolved.layers.placementCount[i]; k++)
    if (p.annotation[k] === ann.id && p.slot[k] === ann.slot[r] && p.lane[k] !== -2 && (lane < 0 || p.lane[k] === lane))
      return true;
  return false;
}

/** Calls `row(ann, r)` for every shown row of track `t`'s item layers, on `split`'s value or `lane` when ≥ 0. */
function eachTrackRow(
  model: ViewModel,
  t: number,
  split: number,
  lane: number,
  row: (ann: Annotation, r: number) => void,
) {
  const { resolved: rv, set, prepared } = model;
  if (rv.tracks.hidden[t]) return;
  const k = rv.tracks.split[t];
  for (let i = 0; i < rv.layers.count; i++) {
    if (rv.layers.track[i] !== t || !lightsItems(model, i)) continue;
    const rowOf = split >= 0 && i === k ? (prepared[i]?.state as SplitFeatures | undefined)?.split?.rowOf : undefined;
    if (split >= 0 && i === k && !rowOf) continue;
    for (let d = rv.layers.dataStart[i]; d < rv.layers.dataStart[i] + rv.layers.dataCount[i]; d++) {
      const ann = set.annotations[rv.data[d]];
      if (!ann) continue;
      const kept = rv.kept[d];
      const of = rowOf?.[d - rv.layers.dataStart[i]];
      for (let r = 0; r < ann.count; r++)
        if ((!kept || kept[r]) && (!of || of[r] === split) && ann.item[r] >= 0 && placed(model, i, ann, r, lane))
          row(ann, r);
    }
  }
}

/** Whether group `g` is `parent` or under it, none of them hidden. */
function under(model: ViewModel, g: number, parent: number): boolean {
  const { parent: up, hidden } = model.resolved.groups;
  for (let at = g; at >= 0; at = up[at]) {
    if (hidden[at]) return false;
    if (at === parent) return true;
  }
  return false;
}

/** Every shown row the description stands for. */
function eachRow(model: ViewModel, s: HeaderSet, row: (ann: Annotation, r: number) => void) {
  if (s.track >= 0) return eachTrackRow(model, s.track, s.split, s.lane, row);
  const tracks = model.resolved.tracks;
  for (let t = 0; t < tracks.count; t++)
    if (tracks.group[t] >= 0 && under(model, tracks.group[t], s.group)) eachTrackRow(model, t, -1, -1, row);
}

// The counting pass's state, module-level so the callback is one function (no closure per landing).
let count = 0;
let length = 0;
const tally = (ann: Annotation, r: number) => {
  const a = ann.id;
  if (!touchedOf[a]) {
    touchedOf[a] = 1;
    if (touchedN === touched.length) {
      const grown = new Int32Array(2 * touched.length);
      grown.set(touched);
      touched = grown;
    }
    touched[touchedN++] = a;
  }
  const row = rowAt[a] + r;
  if (seenRow[row]) return;
  seenRow[row] = 1;
  length += Math.max(1, ann.stop[r] - ann.start[r]);
  const item = itemAt[a] + ann.item[r];
  if (!seenItem[item]) {
    seenItem[item] = 1;
    count++;
  }
};
let into: Extents | null = null;
const sys = { kind: 0, id: -1, row: -1 };
// Read once: a member read in a hover-rate path allocates under the test runner's module wrapper.
const [UNIT, MOLECULE, ROW, COLUMNS] = [K.unit, K.molecule, K.row, K.columns];
/** A row's target (a slot of `t`, or a pair's side b) into `sys`; false: none. */
function target(t: { unit: Int32Array; molecule: Int32Array; alignment: Int32Array; row: Int32Array }, i: number) {
  sys.row = -1;
  if (t.unit[i] >= 0) {
    sys.kind = UNIT;
    sys.id = t.unit[i];
  } else if (t.molecule[i] >= 0) {
    sys.kind = MOLECULE;
    sys.id = t.molecule[i];
  } else if (t.alignment[i] >= 0) {
    sys.kind = t.row[i] >= 0 ? ROW : COLUMNS;
    sys.id = t.alignment[i];
    sys.row = t.row[i];
  } else return false;
  return true;
}
const run = (ann: Annotation, r: number) => {
  const out = into as Extents;
  if (target(ann.slots, ann.slot[r])) Extents.add(out, sys, ann.start[r], ann.stop[r]);
  if (ann.pairs && target(ann.pairs.bTarget, r)) {
    const whole = ann.pairs.bStart[r] < 0;
    Extents.add(out, sys, whole ? 0 : ann.pairs.bStart[r], whole ? 0x7fffffff : ann.pairs.bStop[r]);
  }
};

/**
 * Describes a header hover into `out` (`track`, `split`, `lane`, `group` set by the caller): counts its shown items and
 * covered positions, decides whether it lights and forwards under `params`, and when it lights writes its rows' runs.
 * Allocation-free once the buffers have grown: a header crossed while sweeping is hover-rate.
 */
export function describeHeader(model: ViewModel | null | undefined, params: HeaderHoverParams, out: HeaderSet): void {
  out.extents.count = 0;
  out.items = out.length = 0;
  out.light = out.forward = false;
  out.on = !!model && (out.track >= 0 || out.group >= 0);
  if (!model || !out.on) return;
  P = params;
  count = 0;
  length = 0;
  seenFor(model.set as SeqSet);
  eachRow(model, out, tally);
  out.items = count;
  out.length = length;
  const within = count > 0 && count <= params.maxItems;
  out.light = within;
  out.forward = within && params.forward;
  if (!within) return;
  into = out.extents;
  eachRow(model, out, run);
  into = null;
}

/** The description's items as item loci (per annotation, its shown items), what linked views read. Allocates. */
export function headerLoci(model: ViewModel | null | undefined, params: HeaderHoverParams, s: HeaderSet): Loci | null {
  if (!model || !s.on || s.items === 0) return null;
  P = params;
  const by = new Map<number, number[]>();
  eachRow(model, s, (ann, r) => {
    const list = by.get(ann.id) ?? [];
    if (list.length === 0) by.set(ann.id, list);
    list.push(ann.item[r]);
  });
  const elements = [...by].map(([annotation, items]) => ({
    annotation: annotation as never,
    items: OrderedSet.ofSortedArray(Int32Array.from(new Set(items)).sort()),
  }));
  return elements.length > 0 ? ({ kind: 'item-loci', set: model.set as SeqSet, elements } as ItemLoci) : null;
}

/** Any loci's size: its items (an item loci's, else 0) and the positions its runs cover. Edge-rate: allocates. */
export function lociSize(loci: Loci | null): LociSize | null {
  if (!loci) return null;
  let items = 0;
  if (loci.kind === 'item-loci') for (const e of loci.elements) items += OrderedSet.size(e.items);
  const runs = Extents.create(16);
  Extents.addLoci(runs, loci);
  let length = 0;
  for (let i = 0; i < runs.count; i++) length += Math.max(0, Math.min(runs.stop[i], 0x7fffffff) - runs.start[i]);
  return { items, length };
}
