/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type Location } from '@molstar/seqx-core/query';
import type { PairsData } from '@molstar/seqx-core/resolve';
import { NoColor, type PassName, type PassWriter } from './draw-lists.js';
import { aOf, bOf, clipTo, LitOf, type PairsLay, type PairsPrep } from './pairs-lay.js';
import { type Emphasis, type EmphasisHighlight, type ItemDescription } from './layer-type.js';

// The pairs layer's emphasis (_spec/README.md "Render", The pairs layer): hover, selection and the highlight channels
// over the arcs and dots laid out, and the pairs drawn as items.

// ---- emphasis ----

/** Whether a pair has an end holding position `p` on `system` (a hovered base: its partner lights). */
function endHolds(prep: PairsPrep, k: number, system: string, p: number): boolean {
  const d = prep.data!;
  return (
    (system === aOf(prep, k) && p >= d.aStart[k] && p < d.aEnd[k]) ||
    (system === bOf(prep, k) && p >= d.bStart[k] && p < d.bEnd[k])
  );
}

/** Whether the selection's elements hold a pair (few keys: read in order). */
function selected(prep: PairsPrep, e: Emphasis, k: number): boolean {
  const keys = e.elements?.keys;
  if (!keys) return false;
  const key = prep.data!.keys[k];
  for (let i = 0; i < keys.length; i++) if (keys[i].owner === prep.owner && keys[i].key === key) return true;
  return false;
}

/** Whether a pair's `a` end holds `x` and its `b` end `y`. */
function endsAt(d: PairsData, k: number, x: number, y: number): boolean {
  return x >= d.aStart[k] && x < d.aEnd[k] && y >= d.bStart[k] && y < d.bEnd[k];
}

/** Whether a pair has an end in a Location's ranges (a projected hover: a cell's two positions on linked sections). */
function endIn(prep: PairsPrep, k: number, loc: Location): boolean {
  const d = prep.data!;
  for (let i = 0; i < loc.parts.length; i++) {
    const part = loc.parts[i];
    const a = part.system === aOf(prep, k);
    const b = part.system === bOf(prep, k);
    if (!a && !b) continue;
    const r = part.ranges;
    for (let j = 0; j + 1 < r.length; j += 2)
      if ((a && d.aStart[k] < r[j + 1] && r[j] < d.aEnd[k]) || (b && d.bStart[k] < r[j + 1] && r[j] < d.bEnd[k]))
        return true;
  }
  return false;
}

/**
 * Whether a pair is lit in a pass. The hover: the hovered pair; as arcs, every pair with an end at a hovered position
 * (a base, or either position of a cell: its partners light) or in the hover projected here; as dots, the pair at the
 * hovered cell. The selection: a selected pair.
 */
function lit(lay: PairsLay, pass: PassName, e: Emphasis, k: number, hovered: number): boolean {
  const prep = lay.prep!;
  if (pass !== 'hover') return selected(prep, e, k);
  if (k === hovered) return true;
  const h = e.hover;
  if (h?.element) return false;
  const x = h?.position === undefined ? -1 : Math.floor(h.position);
  const y = h?.yPosition === undefined ? -1 : Math.floor(h.yPosition);
  if (lay.mode === 'dots') {
    if (!h || x < 0 || y < 0 || h.ySystem === undefined) return false;
    return (
      (h.system === prep.aSystem && h.ySystem === prep.bSystem && endsAt(prep.data!, k, x, y)) ||
      (prep.mirror && h.system === prep.bSystem && h.ySystem === prep.aSystem && endsAt(prep.data!, k, y, x))
    );
  }
  if (h && x >= 0 && endHolds(prep, k, h.system, x)) return true;
  if (h && y >= 0 && h.ySystem !== undefined && endHolds(prep, k, h.ySystem, y)) return true;
  const p = e.projected?.hover;
  return p !== undefined && endIn(prep, k, p);
}

/**
 * The hover: the hovered pair (its element), or the pairs at the hovered positions (`lit`), drawn over in the hover's
 * ink; the selection: each selected pair in the accent.
 */
export function drawEmphasis(lay: PairsLay, pass: PassName, out: PassWriter, e: Emphasis) {
  const prep = lay.prep!;
  const tokens = lay.tokens!;
  const h = e.hover;
  let hovered = -1;
  if (pass === 'hover') {
    if (!h && !e.projected?.hover) return;
    if (h?.element) {
      if (h.element.owner !== prep.owner) return;
      hovered = prep.index.get(h.element.key) ?? -1;
      if (hovered < 0) return;
    }
  } else if (!e.elements?.keys.length) return;
  const ink = pass === 'hover' ? tokens.colors['seq-hover'] : tokens.colors['msui-accent'];
  const width = pass === 'hover' ? lay.config!.outline.hover + 1 : lay.config!.outline.selection + 1;
  const g = out.args;
  if (lay.mode === 'arcs') {
    if (lay.clip) clipTo(out, lay.clip);
    for (let i = 0; i < lay.n; i++) {
      if (!lit(lay, pass, e, lay.pair[i], hovered)) continue;
      g[4] = NoColor;
      g[5] = ink;
      g[6] = width;
      g[7] = 0;
      out.curveArgs(lay.curve, 8 * i);
    }
    if (lay.clip) out.unclip();
    return;
  }
  for (let i = 0; i < lay.nd; i++) {
    if (!lit(lay, pass, e, lay.dotPair[i], hovered)) continue;
    g[0] = lay.dots[4 * i];
    g[1] = lay.dots[4 * i + 1];
    g[2] = lay.dots[4 * i + 2];
    g[3] = lay.dots[4 * i + 3];
    g[4] = NoColor;
    g[5] = ink;
    g[6] = width - 1;
    g[7] = 0;
    g[8] = 0;
    out.rectArgs();
  }
  const p = lay.seg;
  for (let i = 0; i < lay.ns; i++) {
    if (!lit(lay, pass, e, lay.segPair[i], hovered)) continue;
    p[0] = lay.segs[5 * i];
    p[1] = lay.segs[5 * i + 1];
    p[2] = lay.segs[5 * i + 2];
    p[3] = lay.segs[5 * i + 3];
    g[4] = NoColor;
    g[5] = ink;
    g[6] = lay.segs[5 * i + 4] + 2;
    g[7] = 0;
    out.pathArgs(p, 2, false);
  }
}

/**
 * Each pair's highlight channel, made again when the channels change (by identity: the viewer keeps the list): the
 * last channel whose element keys name the pair (a key names its first pair, as `index` does).
 */
function litOf(lay: PairsLay, highlights: readonly EmphasisHighlight[]): Int32Array {
  const prep = lay.prep!;
  let l = LitOf.get(lay);
  if (!l) LitOf.set(lay, (l = { for: undefined, prep: undefined, lit: new Int32Array(0) }));
  if (l.for === highlights && l.prep === prep) return l.lit;
  const n = prep.data?.n ?? 0;
  if (l.lit.length < n) l.lit = new Int32Array(n);
  const lit = l.lit;
  lit.fill(0);
  for (let c = 0; c < highlights.length; c++) {
    const keys = highlights[c].elements?.keys;
    if (!keys) continue;
    for (let i = 0; i < keys.length; i++) {
      if (keys[i].owner !== prep.owner) continue;
      const k = prep.index.get(keys[i].key);
      if (k !== undefined) lit[k] = c + 1;
    }
  }
  l.for = highlights;
  l.prep = prep;
  return lit;
}

/**
 * The highlight channels' pairs (_spec/README.md "Render", The pairs layer, Highlight): each pair a channel lists by its
 * key drawn over in the channel's color, as the selection's are in the accent (arcs, chords, dots and stretches).
 */
export function drawHighlights(lay: PairsLay, out: PassWriter, highlights: readonly EmphasisHighlight[]) {
  if (!lay.prep?.data) return;
  let any = false;
  for (let c = 0; c < highlights.length && !any; c++) any = !!highlights[c].elements?.keys.length;
  if (!any) return;
  const lit = litOf(lay, highlights);
  const width = lay.config!.outline.selection + 1;
  const g = out.args;
  if (lay.mode === 'arcs') {
    if (lay.clip) clipTo(out, lay.clip);
    for (let i = 0; i < lay.n; i++) {
      const c = lit[lay.pair[i]];
      if (!c) continue;
      g[4] = NoColor;
      g[5] = highlights[c - 1].color;
      g[6] = width;
      g[7] = 0;
      out.curveArgs(lay.curve, 8 * i);
    }
    if (lay.clip) out.unclip();
    return;
  }
  for (let i = 0; i < lay.nd; i++) {
    const c = lit[lay.dotPair[i]];
    if (!c) continue;
    g[0] = lay.dots[4 * i];
    g[1] = lay.dots[4 * i + 1];
    g[2] = lay.dots[4 * i + 2];
    g[3] = lay.dots[4 * i + 3];
    g[4] = NoColor;
    g[5] = highlights[c - 1].color;
    g[6] = width - 1;
    g[7] = 0;
    g[8] = 0;
    out.rectArgs();
  }
  const p = lay.seg;
  for (let i = 0; i < lay.ns; i++) {
    const c = lit[lay.segPair[i]];
    if (!c) continue;
    p[0] = lay.segs[5 * i];
    p[1] = lay.segs[5 * i + 1];
    p[2] = lay.segs[5 * i + 2];
    p[3] = lay.segs[5 * i + 3];
    g[4] = NoColor;
    g[5] = highlights[c - 1].color;
    g[6] = lay.segs[5 * i + 4] + 2;
    g[7] = 0;
    out.pathArgs(p, 2, false);
  }
}

/** The pairs drawn (the arcs, the dots) as items: a pair's label (else its key) and its first end's extent. */
export function pairItems(lay: PairsLay): ItemDescription[] {
  const d = lay.prep?.data;
  if (!d) return [];
  const out: ItemDescription[] = [];
  const add = (k: number, x0: number, x1: number) => {
    const key = d.keys[k];
    out.push({
      key,
      label: d.labels?.[k] || String(key),
      start: d.aStart[k],
      end: d.aEnd[k],
      x0: Math.min(x0, x1),
      x1: Math.max(x0, x1),
    });
  };
  if (lay.mode === 'arcs') for (let i = 0; i < lay.n; i++) add(lay.pair[i], lay.curve[8 * i], lay.curve[8 * i + 6]);
  else if (lay.mode === 'dots')
    for (let i = 0; i < lay.nd; i++) add(lay.dotPair[i], lay.dots[4 * i], lay.dots[4 * i] + lay.dots[4 * i + 2]);
  return out;
}
