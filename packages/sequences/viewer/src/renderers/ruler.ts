/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The `ruler` renderer (seqspec:viewer#extension-points): ticks on round numbers of each block's numbering, thinned by
// the labels' width; with several blocks, each named at its start, kept in view; no band when wrapped (its fallback
// omits it). Labels come from the layer's sources (`of` through each placement's path, a list prefixed with each
// system's label,
// a `numbering` list's first that labels a position), else from the section's numbering and `numberingPath`, read live
// from the view (a numbering switch prepares nothing). Positions are walked by the breaks' runs (`shownToRuns`): a
// break is labeled and marked; numbers that run by one (an offset or ordinal numbering on the block itself) step
// from round number to round number; others come from tick lists by power of ten, built once per numbering, path and
// prepared state in a draw-side cache, or below a step of 10 position by position. Through a codon step ticks fall at
// residue starts (a number differing from the previous position's). Side by side, labels never overlap across blocks: a
// block's name and numbers stay left of the next block's start. Where a source residue has several copies on the block
// (a construct's segments over one source), a label past the first says which: "Ub 48 (copy 2)".
import { Coords, NO_LABEL, Numbering, type PositionLabel } from "@molstar/seq-core/model";
import {
  copyAt,
  type DrawView,
  kernelEntry,
  numberText,
  RectBatch,
  type RulerState,
  SkinToken,
  shownHull,
  sourceAt,
  TextAlign,
  TextBaseline,
  TextBatch,
} from "../draw/api.ts";

const rects = new RectBatch(128);
const nums = new TextBatch(64);
const names = new TextBatch(16);
let runs = new Int32Array(64);
const label: PositionLabel = { num: 0, suffix: -1 as never, name: -1 as never };
/** The last label read: its number and its source (-1: the section's numbering). */
const at = new Int32Array(2);

/** Numbering n's label at its index j into `at[0]` (none: ordinals); false: unlabeled, a suffix or a name. */
function numbered(n: Numbering | null, j: number): boolean {
  if (j < 0) return false;
  if (!n) {
    at[0] = j + 1;
    return true;
  }
  if (!Numbering.labelOf(n, j as never, label) || label.suffix >= 0 || label.num === NO_LABEL) return false;
  at[0] = label.num;
  return true;
}
/** Block b's label at block index i into `at`: the layer's sources in order, else the section's numbering (live). */
function labelAt(v: DrawView, s: RulerState, b: number, i: number): boolean {
  const srcs = s.sources[b];
  at[1] = -1;
  if (!srcs) {
    const path = v.numberingPath[b] ?? null;
    return numbered(v.numbering[b] ?? null, path ? sourceAt(path, i) : i);
  }
  for (let q = 0; q < srcs.length; q++) {
    const src = srcs[q];
    const j = src.path ? sourceAt(src.path, i) : i;
    at[1] = q;
    if (src.numberings.length === 0 ? numbered(null, j) : false) return true;
    for (let r = 0; r < src.numberings.length; r++) if (numbered(src.numberings[r], j)) return true;
  }
  return false;
}
/** A residue start labeled: a label differing from the previous block index's (a codon's first nucleotide). */
function startAt(v: DrawView, s: RulerState, b: number, i: number): boolean {
  if (!labelAt(v, s, b, i)) return false;
  const n = at[0];
  const q = at[1];
  const same = i > 0 && labelAt(v, s, b, i - 1) && at[0] === n;
  at[0] = n;
  at[1] = q;
  return !same;
}

/**
 * Per numbering, path or prepared sources: block indices whose label is a multiple of 10^(level + 1), with it, from
 * block indices `range` [lo, hi) (what the blocks drawn show and their sources reach: never a chromosome's domain for
 * a protein's residue numbers, GM-03; seqspec:viewer#performance).
 */
interface Ticks {
  readonly numbering: Numbering | null;
  readonly lists: Int32Array[];
  readonly nums: Int32Array[];
  readonly srcs: Int8Array[];
  readonly range: Int32Array;
}
/** Per key, per block: two pieces of one unit share their numbering but never a hull across the gap between them (LC-33). */
const cache = new WeakMap<object, (Ticks | undefined)[]>();
function ticksOf(v: DrawView, s: RulerState, b: number): Ticks {
  const key = (s.sources[b] ?? v.numberingPath[b] ?? v.numbering[b]) as object;
  const numbering = s.sources[b] ? null : (v.numbering[b] ?? null);
  let blocks = cache.get(key);
  if (!blocks) {
    blocks = [];
    cache.set(key, blocks);
  }
  let t = blocks[b];
  if (t?.numbering !== numbering) {
    t = { numbering, lists: [], nums: [], srcs: [], range: new Int32Array(2) };
    blocks[b] = t;
  }
  return t;
}
const need = new Int32Array(2);
/**
 * Into `need` the block indices block b's ticks come from: its shown piece's hull, within its sources' paths' reach on
 * the block (a source on the block itself: no bound). Integers only (each frame checks the lists still cover it).
 */
function needOf(v: DrawView, s: RulerState, b: number): void {
  const ax = v.axis;
  if (!ax || ax.shown[b] === 0) {
    need[0] = need[1] = 0;
    return;
  }
  shownHull(v, b, 0, ax.shown[b], need);
  const srcs = s.sources[b];
  if (!srcs) return;
  let lo = 0x7fffffff;
  let hi = -0x80000000;
  for (let q = 0; q < srcs.length; q++) {
    const p = srcs[q].path;
    if (!p) return;
    // The map's side on the block: its `to` where it runs from the source, else its `from`.
    const toBlock = Coords.sameSystem(p.map.from, p.from);
    lo = Math.min(lo, toBlock ? p.map.toStart : p.map.fromStart);
    hi = Math.max(hi, toBlock ? p.map.toStop : p.map.fromStop);
  }
  need[0] = Math.max(need[0], lo);
  need[1] = Math.max(need[0], Math.min(need[1], hi));
}
/**
 * The tick list of 10^(level + 1), built on first use from the positions block b needs (level 0) or the list below;
 * all built again over the union where the block needs positions they do not cover. Edge-rate.
 */
function listOf(v: DrawView, s: RulerState, b: number, t: Ticks, level: number): number {
  needOf(v, s, b);
  if (need[0] < t.range[0] || need[1] > t.range[1]) widen(t);
  while (v.axis && t.lists.length <= level) build(v, s, b, t);
  return level;
}
/** The lists dropped, their range now covering `need` too (edge-rate). */
function widen(t: Ticks): void {
  const empty = t.range[1] <= t.range[0];
  t.range[0] = empty ? need[0] : Math.min(t.range[0], need[0]);
  t.range[1] = empty ? need[1] : Math.max(t.range[1], need[1]);
  t.lists.length = t.nums.length = t.srcs.length = 0;
}
/** The next tick list (edge-rate, apart: a closure in `listOf` would make its context on every frame). */
function build(v: DrawView, s: RulerState, b: number, t: Ticks) {
  const ax = v.axis;
  if (ax) {
    const p = 10 ** (t.lists.length + 1);
    const is: number[] = [];
    const ns: number[] = [];
    const qs: number[] = [];
    const keep = (i: number, n: number, q: number) => {
      if (n % p !== 0) return;
      is.push(i);
      ns.push(n);
      qs.push(q);
    };
    const below = t.lists.length - 1;
    if (below < 0)
      for (let i = t.range[0]; i < t.range[1]; i++) {
        if (startAt(v, s, b, i)) keep(i, at[0], at[1]);
      }
    else for (let x = 0; x < t.lists[below].length; x++) keep(t.lists[below][x], t.nums[below][x], t.srcs[below][x]);
    t.lists.push(Int32Array.from(is));
    t.nums.push(Int32Array.from(ns));
    t.srcs.push(Int8Array.from(qs));
  }
}

/** Px a position (a double crosses no call: V8 boxes it). */
const perPx = new Float64Array(1);
/**
 * The least round step (1, 2, 2.5 or 5 × 10^k positions, 2.5 only from 25) whose ticks are `room` px apart at
 * `perPx[0]` px a position: any length, so a long block's ticks thin out rather than crowd.
 */
function stepFor(room: number): number {
  const per = perPx[0];
  for (let p = 1; p < 1e12; p *= 10) {
    if (p * per >= room) return p;
    if (2 * p * per >= room) return 2 * p;
    if (p >= 10 && 2.5 * p * per >= room) return 2.5 * p;
    if (5 * p * per >= room) return 5 * p;
  }
  return 1e12;
}

/**
 * What a block's ticks share while they are placed: the step, the label font, its y, the clear edge, the limit; on rings
 * the minor step (0: none).
 */
const frame = new Int32Array(8);
const F = { step: 0, font: 1, y: 2, clear: 3, limit: 4, bottom: 5, broke: 6, minor: 7 } as const;

/** A tick at block b's shown k with the label in `at`, its label where it fits (break: labeled whatever its number). */
function tick(v: DrawView, s: RulerState, b: number, k: number) {
  if (k === frame[F.broke]) return;
  const dpr = v.dpr;
  const x = Math.round(v.xd(b, k) + (v.cell * v.span[b] * dpr) / 2) | 0;
  const th = Math.round(4 * dpr) | 0;
  rects.begin(SkinToken.inkQuiet);
  rects.push(x, frame[F.bottom] - th, 1, th);
  const srcs = s.sources[b];
  const prefix = srcs && at[1] >= 0 ? srcs[at[1]].prefix : "";
  const copy = copyAt(srcs && at[1] >= 0 ? srcs[at[1]].path : v.numberingPath[b], v.indexOf(b, k as never));
  const text = numberText(at[0], prefix, copy);
  const w = prefix || copy > 1 ? v.text.label(text, frame[F.font]).width : v.text.numberWidth(frame[F.font], at[0]);
  const half = Math.round((w * dpr) / 2) | 0;
  // A label cut by the data area's edge (rings have none), or reaching the next block's name, is left out (its tick stays).
  if (x - half < frame[F.clear] || (!v.ring && x - half < Math.round(v.left * dpr)) || x + half > frame[F.limit])
    return;
  const t = nums.add(text);
  nums.x[t] = x;
  nums.y[t] = frame[F.y];
  frame[F.clear] = (x + half + Math.round(6 * dpr)) | 0;
}

/** Block b's shown run [ks, ks + i1 − i0) over block indices [i0, i1): its round numbers in shown order. */
function run(v: DrawView, s: RulerState, b: number, ks: number, i0: number, i1: number) {
  const step = frame[F.step];
  const asc = i1 - i0 === 1 || v.indexOf(b, ks as never) === i0;
  // Shown k of block index i: ks + i − i0 ascending, ks + i1 − 1 − i descending.
  const o = asc ? ks - i0 : ks + i1 - 1;
  const sg = asc ? 1 : -1;
  const n = s.sources[b] ? null : (v.numbering[b] ?? null);
  if (!s.sources[b] && !v.numberingPath[b] && (n === null || (n.method === "offset" && !n.skipZero))) {
    // Numbers that run by one: from round number to round number.
    const d = n ? n.first : 1;
    at[1] = -1;
    const m = frame[F.minor];
    if (m > 0) {
      const th = Math.round(2 * v.dpr) | 0;
      rects.begin(SkinToken.inkQuiet);
      for (let i = i0 + ((m - ((((i0 + d) % m) + m) % m)) % m); i < i1; i += m)
        if ((i + d) % step !== 0) {
          const x = Math.round(v.xd(b, o + sg * i) + (v.cell * v.span[b] * v.dpr) / 2) | 0;
          rects.push(x, frame[F.bottom] - th, 1, th);
        }
    }
    if (asc)
      for (let i = i0 + ((step - ((((i0 + d) % step) + step) % step)) % step); i < i1; i += step) {
        at[0] = i + d;
        tick(v, s, b, o + sg * i);
      }
    else
      for (let i = i1 - 1 - ((((i1 - 1 + d) % step) + step) % step); i >= i0; i -= step) {
        at[0] = i + d;
        tick(v, s, b, o + sg * i);
      }
    return;
  }
  let p = 1;
  while (step % (p * 10) === 0) p *= 10;
  if (p === 1) {
    for (let x = 0; x < i1 - i0; x++) {
      const i = asc ? i0 + x : i1 - 1 - x;
      if (startAt(v, s, b, i) && at[0] % step === 0) tick(v, s, b, o + sg * i);
    }
    return;
  }
  const t = ticksOf(v, s, b);
  const L = listOf(v, s, b, t, Math.round(Math.log10(p)) - 1);
  const list = t.lists[L];
  let lo = 0;
  for (let hi = list.length; lo < hi; ) {
    const mid = (lo + hi) >> 1;
    if (list[mid] < i0) lo = mid + 1;
    else hi = mid;
  }
  let hi = lo;
  while (hi < list.length && list[hi] < i1) hi++;
  for (let x = 0; x < hi - lo; x++) {
    const e = asc ? lo + x : hi - 1 - x;
    if (t.nums[L][e] % step !== 0) continue;
    at[0] = t.nums[L][e];
    at[1] = t.srcs[L][e];
    tick(v, s, b, o + sg * list[e]);
  }
}

function draw(v: DrawView, s: RulerState) {
  const ax = v.axis;
  if (!ax) return;
  const tm = v.text;
  const dpr = v.dpr;
  const font = tm.font("ui", 11);
  const bold = tm.font("ui", 11, 600);
  const left = Math.round(v.left * dpr) | 0;
  const right = v.ring ? 0x3fffffff : Math.round(v.right * dpr) | 0;
  frame[F.bottom] = Math.round((v.top + v.height) * dpr) - 1;
  frame[F.y] = Math.round((v.top + 10) * dpr);
  frame[F.font] = font;
  rects.reset();
  nums.reset();
  names.reset();
  nums.style(font, TextAlign.center, TextBaseline.middle);
  nums.begin(SkinToken.inkQuiet);
  names.style(bold, TextAlign.left, TextBaseline.middle);
  names.begin(SkinToken.ink);
  // The right end of the last label drawn (side by side: across blocks).
  let clear = -0x40000000;
  for (let b = 0; b < s.sources.length; b++) {
    const lo = v.visLo[b];
    const hi = v.visHi[b];
    if (ax.stacked) clear = -0x40000000;
    if (lo >= hi || v.span[b] > 1 || (v.onlyBlock >= 0 && b !== v.onlyBlock)) continue;
    rects.begin(SkinToken.grid);
    rects.push(v.xd(b, lo), frame[F.bottom], v.xd(b, hi) - v.xd(b, lo), 1);
    // An alignment's insert columns: a short mark under each, on the axis line.
    const ins = ax.alignment[b] >= 0 ? (ax.set.alignments[ax.alignment[b]]?.insertColumn ?? null) : null;
    if (ins !== null && v.cell * dpr >= 2) {
      rects.begin(SkinToken.inkQuiet);
      const th = Math.max(2, Math.round(3 * dpr)) | 0;
      for (let q = lo; q < hi; q++)
        if (ins[ax.indexOf(b, q as never)] === 1)
          rects.push(v.xd(b, q) + 1, frame[F.bottom] - th + 1, Math.max(1, v.xd(b, q + 1) - v.xd(b, q) - 2), th);
    }
    // Labels of this block end before the next block's name.
    const limit = (!ax.stacked && b + 1 < ax.count ? Math.min(right, v.xd(b + 1, 0) - Math.round(4 * dpr)) : right) | 0;
    const name = s.names[b];
    if (ax.count > 1 && name) {
      const w = Math.ceil(tm.label(name, bold).width * dpr) | 0;
      let x = (Math.max(v.xd(b, 0), left) + Math.round(2 * dpr)) | 0;
      // Past the limit (the last block's name at the data area's edge), it moves back into the gap before the block.
      const back = (!ax.stacked && b > 0 ? v.xd(b - 1, ax.shown[b - 1]) : left) + (Math.round(2 * dpr) | 0);
      if (x + w > limit) x = Math.max(back, limit - w) | 0;
      if (x >= clear && x + w <= limit) {
        const t = names.add(name);
        names.x[t] = x;
        names.y[t] = frame[F.y];
        clear = (x + w + Math.round(6 * dpr)) | 0;
      }
    }
    // The step fits the wider of the piece's end labels (a reversed piece ends on its lowest number).
    const wide = labelAt(v, s, b, ax.indexOf(b, 0 as never)) ? tm.numberWidth(font, at[0]) : 0;
    const last = labelAt(v, s, b, ax.indexOf(b, (ax.shown[b] - 1) as never)) ? at[0] : ax.length[b];
    perPx[0] = v.cell * v.span[b];
    frame[F.step] = stepFor(Math.max(wide, tm.numberWidth(font, last)) + 14);
    // On rings (seqspec:viewer#the-render-backend) the numbers sit outside the ring, so it has room for unlabeled ticks
    // between them: a fifth of the step (a half of 2 × 10^k), at least 4 px apart (inline: a helper's result would be
    // boxed).
    const st = frame[F.step];
    const mi = st % 5 === 0 ? st / 5 : st % 2 === 0 ? st / 2 : 0;
    frame[F.minor] = v.ring && mi >= 1 && mi * perPx[0] >= 4 ? mi : 0;
    frame[F.broke] = -1;
    frame[F.clear] = clear;
    frame[F.limit] = limit;
    let n = v.shownToRuns(b, lo, hi, runs);
    if (2 * n > runs.length) {
      runs = new Int32Array(4 * n);
      n = v.shownToRuns(b, lo, hi, runs);
    }
    for (let r = 0, ks = lo; r < n; r++) {
      const i0 = runs[2 * r];
      const i1 = runs[2 * r + 1];
      // A break: marked, and its first position labeled whatever its number.
      if (r > 0 && v.breakAt(b, ks as never)) {
        rects.begin(SkinToken.inkFaint);
        rects.push(v.xd(b, ks), Math.round(v.top * dpr) | 0, 1, (frame[F.bottom] - Math.round(v.top * dpr)) | 0);
        const first = v.indexOf(b, ks as never);
        if (labelAt(v, s, b, first)) tick(v, s, b, ks);
        frame[F.broke] = ks;
      }
      run(v, s, b, ks, i0, i1);
      ks += i1 - i0;
    }
    clear = frame[F.clear];
  }
  if (rects.n > 0) v.backend.rects(rects);
  if (names.n > 0) v.backend.text(names);
  if (nums.n > 0) v.backend.text(nums);
}

/** No row band where it is hovered (seqspec:interaction#hover): the ruler is the axis, not a track's data. */
const rowHover = () => false;
export const RulerRenderer = kernelEntry<RulerState>("ruler", 26, { miniature: false, draw, tiny: draw, rowHover });
