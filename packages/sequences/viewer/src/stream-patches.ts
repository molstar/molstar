/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A region's patch as it arrives (seqspec:sources#data-patches): an alignment's tile, a 1-D region's values or items, and
// the bytes each holds. Arrival rate.

import { type DataPatch, TILE, type Tile } from '@molstar/seq-core/format';
import { Paths } from '@molstar/seq-core/model';
import type { Features } from './prepare.js';
import { prefixReach } from './pyramids.js';
import type { StreamItems } from './stream-contract.js';
import type { Region } from './stream-layer.js';
import * as Store from './stream-layer.js';
import * as Levels from './stream-levels.js';

// Read once into module constants: a read through an import is a property read of the module per use, which
// the test runner's module wrapper turns into an allocation on per-frame paths.
const { FINE } = Store;
const { D1, LV } = Levels;

const num = (x: unknown) => (x == null ? Number.NaN : Number(x));

/** An alignment's `bins` patch as a 256 × 256 tile: as sent when whole, else its samples placed row-major. */
export function tileOf(p: DataPatch, lx: number): Tile {
  const b = p.kind === 'bins' ? p.bins : null;
  const cls = b?.class;
  const share = b?.share;
  if (
    cls instanceof Uint8Array &&
    share instanceof Uint8Array &&
    cls.length === TILE * TILE &&
    share.length === cls.length
  )
    return { cls, share };
  const t: Tile = { cls: new Uint8Array(TILE * TILE), share: new Uint8Array(TILE * TILE) };
  const across = Math.ceil((p.stop - p.start) / 2 ** lx);
  for (let k = 0; cls && share && across > 0 && k < cls.length; k++) {
    const at = Math.floor(k / across) * TILE + (k % across);
    t.cls[at] = Number(cls[k]);
    t.share[at] = Number(share[k]);
  }
  return t;
}

/**
 * A 1-D patch held by region r: letters as upper-case codes; `values` as per-bin `min, max, mean` triplets (`bins`) or
 * per-position values (`items`: a dense `value` column, or rows with `start`, `stop` and `value`); an items layer's
 * through `heldItems`. Its bytes; -1: loaded with nothing (empty; a kind not asked failed before, `take`).
 */
export function held1(r: Region, p: DataPatch, values: boolean): number {
  if (r.owner.isItems) return heldItems(r, p);
  if (p.kind === 'letters') {
    // Upper case by ASCII rules (a 2bit file's soft-masked repeats are lower case), as a molecule's code.
    const t = p.letters;
    const codes = new Uint8Array(t.length);
    for (let j = 0; j < t.length; j++) {
      const c = t.charCodeAt(j) & 127;
      codes[j] = c >= 97 && c <= 122 ? c - 32 : c;
    }
    r.codes = codes;
    return codes.length > 0 ? codes.length : -1;
  }
  if (values && p.kind === 'bins') {
    const { min, max, mean } = p.bins;
    const n = mean?.length ?? 0;
    const v = new Float64Array(3 * n);
    for (let j = 0; j < n; j++) {
      v[3 * j] = num(min?.[j] ?? mean?.[j]);
      v[3 * j + 1] = num(max?.[j] ?? mean?.[j]);
      v[3 * j + 2] = num(mean?.[j]);
    }
    r.vals = v;
    r.stride = 3;
    return n > 0 ? v.byteLength : -1;
  }
  if (values && p.kind === 'items') {
    const t = p.items;
    const v = new Float64Array(Math.max(0, p.stop - p.start)).fill(Number.NaN);
    const value = t.value;
    if (value && !t.start) for (let j = 0; j < v.length && j < value.length; j++) v[j] = num(value[j]);
    for (let k = 0; t.start && t.stop && value && k < t.start.length; k++)
      v.fill(num(value[k]), Math.max(0, Number(t.start[k]) - p.start), Math.max(0, Number(t.stop[k]) - p.start));
    r.vals = v;
    r.stride = 1;
    return value && value.length > 0 ? v.byteLength : -1;
  }
  return -1;
}

/**
 * An items region's patch held (arrival-rate): `count` bins as they came, or its rows as a feature table over the
 * blocks that show its system (`itemsOf`, and again folded where the layer has a `collapse` field). Its change redraws
 * what its items reach. Its bytes; -1: empty.
 */
function heldItems(r: Region, p: DataPatch): number {
  const c = r.owner;
  const lvl = ((r.level - D1) % LV) - 1;
  if (p.kind === 'bins') {
    const n = p.bins.count?.length ?? 0;
    const bins = Float64Array.from({ length: n }, (_, j) => num((p.bins.count as ArrayLike<unknown>)[j]) || 0);
    // The zoom gate (seqspec:sources#streamed-sources): a data request answered with a listed level's bins, held as
    // that level's.
    if (lvl < 0 && p.level >= 0 && p.level < c.widths.length) r.gate = p.level;
    const width = r.gate >= 0 ? c.widths[r.gate] : lvl < 0 || lvl === FINE ? 1 : c.widths[lvl];
    const hint = r.gate >= 0 ? 2 : p.estimate === true ? 1 : 0;
    r.items = { features: null, first: new Int32Array(0), bins, width, table: null, item: null, ends: null, hint };
    return n > 0 ? bins.byteLength : -1;
  }
  if (p.kind !== 'items') return -1;
  const t = p.items;
  const own = itemsOf(r, t, null);
  // A layer with a `collapse` field holds its items folded too, beside these (one lane a group; a draw picks).
  const fold = c.source.collapse ? itemsOf(r, t, c.source.collapse) : null;
  r.items = own.items;
  r.folded = fold?.items ?? null;
  r.lo = own.lo;
  r.hi = own.hi;
  const m = own.count + (fold?.count ?? 0);
  return m > 0 ? 48 * m : -1;
}

/**
 * One held view of a region's rows: as items (`group` null: rows of one `id` are one item's fragments, joined), or
 * folded (rows sharing the `group` field's value are one item, every row of every member a fragment). An item keeps the
 * layer's lane for its key (its `id`, else its label and extent; a group, its value), a new one takes the first lane
 * none of this region's items overlap.
 */
function itemsOf(r: Region, t: Record<string, ArrayLike<unknown>>, group: string | null) {
  const c = r.owner;
  const s = ((r.level - D1) / LV) | 0;
  const lvl = ((r.level - D1) % LV) - 1;
  const n = t.start?.length ?? 0;
  const col = (name: string, j: number) => {
    const x = t[name];
    return x && j < x.length && x[j] != null ? String(x[j]) : null;
  };
  const laneOf = group === null ? c.laneOf : c.foldedLaneOf;
  // Items by key: their rows and extent; the region holds each key's lane (`discard` lets go).
  const by = new Map<string, { rows: number[]; s: number; e: number; lane: number; id: number }>();
  // A folded row's member (its own item): where that starts, for the rule that a region draws what starts in it.
  const member = new Map<string, number>();
  const memberOf: string[] = [];
  for (let j = 0; j < n; j++) {
    const s0 = Number(t.start[j]);
    const e0 = Math.max(s0, Number(t.stop?.[j] ?? s0 + 1));
    const item = col('id', j) ?? `${col('label', j)}|${s0}|${e0}`;
    let key = item;
    if (group !== null) {
      memberOf.push(item);
      member.set(item, Math.min(member.get(item) ?? s0, s0));
      const g = col(group, j);
      if (g !== null) key = `group|${g}`;
    }
    const x = by.get(key);
    if (x) {
      x.rows.push(j);
      x.s = Math.min(x.s, s0);
      x.e = Math.max(x.e, e0);
      continue;
    }
    const held = laneOf.get(key);
    if (held) held.held++;
    by.set(key, { rows: [j], s: s0, e: e0, lane: held?.lane ?? -1, id: held?.id ?? c.ids++ });
  }
  const keys = [...by.keys()];
  if (group === null) r.keys = keys;
  else r.foldedKeys = keys;
  const list = [...by.entries()].sort((a, b) => a[1].s - b[1].s);
  const taken: [number, number][][] = [];
  const place = (x: { s: number; e: number; lane: number }) => {
    taken[x.lane] ??= [];
    taken[x.lane].push([x.s, Math.max(x.e, x.s + 1)]);
  };
  for (const [, x] of list) if (x.lane >= 0) place(x);
  for (const [key, x] of list) {
    if (x.lane >= 0) continue;
    let l = 0;
    while (taken[l]?.some(([a, e]) => a < Math.max(x.e, x.s + 1) && x.s < e)) l++;
    x.lane = l;
    laneOf.set(key, { lane: l, id: x.id, held: 1 });
    place(x);
  }
  // Fragments on each block that shows the system, by (block, lane, start); a variant's place in its position's stack.
  const sys = c.source.systems as Int32Array;
  const at = new Map<string, number>();
  const frags: {
    b: number;
    s: number;
    e: number;
    lane: number;
    row: number;
    first: number;
    last: number;
    kind: number;
    id: number;
    own: number;
  }[] = [];
  const frag = (b: number, s0: number, e0: number, x: (typeof list)[number][1], j: number) => {
    const kind = at.get(`${b}|${s0}`) ?? 0;
    at.set(`${b}|${s0}`, kind + 1);
    frags.push({
      b,
      s: s0,
      e: e0,
      lane: x.lane,
      row: j,
      first: x.s,
      last: x.e,
      kind,
      id: x.id,
      own: group === null ? x.s : (member.get(memberOf[j]) ?? x.s),
    });
  };
  // Through a block's path (ST7): a row's parts that reach the block, in its indices (an exon's coding part, a
  // variant's residue); one that reaches nothing is not on it.
  const parts: number[] = [];
  const part = (a: number, e: number, acc: number[]) => {
    acc.push(a, e);
  };
  for (let b = 0; b < sys.length; b++) {
    if (sys[b] !== s) continue;
    const path = c.pathOf(b);
    for (const [, x] of list)
      for (const j of x.rows) {
        const s0 = Number(t.start[j]);
        const e0 = Math.max(s0, Number(t.stop?.[j] ?? s0 + 1));
        if (!path) {
          frag(b, s0, e0, x, j);
          continue;
        }
        parts.length = 0;
        Paths.forEachTo(path, s0, Math.max(e0, s0 + 1), part, parts);
        for (let q = 0; q < parts.length; q += 2) frag(b, parts[q], e0 === s0 ? parts[q] : parts[q + 1], x, j);
      }
  }
  frags.sort((x, y) => x.b - y.b || x.lane - y.lane || x.s - y.s);
  const m = frags.length;
  const lanes = taken.length || 1;
  if (group === null) c.lanes = Math.max(c.lanes, lanes);
  else c.foldedLanes = Math.max(c.foldedLanes, lanes);
  const I = (f: (x: (typeof frags)[number]) => number) => Int32Array.from(frags, f);
  const alt = (j: number) => (col('alt', j) ?? '').toUpperCase().charCodeAt(0) || 0;
  // An item's label on one fragment (a BED12's blocks each carry it): drawn once, not per exon (ST4b). An item's, its
  // first; a group's, its widest (the first of its members' exons is often too narrow for a name).
  const carrier = new Map<string, number>();
  frags.forEach((x, i) => {
    const key = `${x.b}|${x.id}`;
    const q = carrier.get(key);
    // Through a path an item's first fragment on the block is its lowest (its start is a system position).
    const first = c.pathOf(x.b) === null ? q === undefined && x.s === x.first : q === undefined || x.s < frags[q].s;
    if (group === null ? first : q === undefined || x.e - x.s > frags[q].e - frags[q].s) carrier.set(key, i);
  });
  const f: Features = {
    count: m,
    block: I((x) => x.b),
    start: I((x) => x.s),
    stop: I((x) => x.e),
    annotation: new Int32Array(m).fill(-1),
    row: I((x) => x.row),
    kind: I((x) => x.kind),
    code: I((x) => alt(x.row)),
    lane: I((x) => x.lane),
    lanes,
    on: new Int32Array(m).fill(-1),
    hasRows: false,
    seg: new Int32Array(sys.length * lanes + 1),
    reach: new Int32Array(m),
    joinFrom: new Int32Array(m).fill(-1),
    label: frags.map((x, i) => (carrier.get(`${x.b}|${x.id}`) === i ? col('label', x.row) : null)),
    split: null,
    fields: { table: t, theme: null, keys: null, glyphTheme: null, glyphs: null },
  };
  for (let g = 0, k = 0; g < f.seg.length; g++) {
    while (k < m && f.block[k] * lanes + f.lane[k] < g) k++;
    f.seg[g] = k;
  }
  for (let g = 0; g + 1 < f.seg.length; g++) prefixReach(f.start, f.stop, f.seg[g], f.seg[g + 1], f.reach);
  // A hairline from an item's fragment reaching furthest so far, across each gap (exons of a gene).
  const furthest = new Map<string, number>();
  for (let k = 0; k < m; k++) {
    const key = `${f.block[k]}|${frags[k].first}|${f.lane[k]}`;
    const q = furthest.get(key);
    if (q !== undefined && f.stop[q] < f.start[k]) f.joinFrom[k] = q;
    if (q === undefined || f.stop[k] >= f.stop[q]) furthest.set(key, k);
  }
  const first = I((x) => x.first);
  const items: StreamItems = {
    features: f,
    first,
    ...(group === null ? {} : { own: I((x) => x.own) }),
    bins: null,
    width: 1,
    table: t,
    item: I((x) => x.id),
    ends: I((x) => x.last),
    hint: 0,
  };
  let lo = r.tx * c.spanOf(lvl);
  let hi = (r.tx + 1) * c.spanOf(lvl);
  // What its items reach, in the system's positions (a fragment through a path is in block indices).
  for (const x of frags) {
    lo = Math.min(lo, x.first);
    hi = Math.max(hi, c.pathOf(x.b) === null ? x.e : x.last);
  }
  return { items, count: m, lo, hi };
}
