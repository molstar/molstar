/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Helpers of the small layers' tests: a layer's layout out of a frame, a pass's paths and rects read back, a hit at a
// point of a track.

import type { Node } from '@molstar/core/repr/tree';
import { DataStore, type FormatRuntime } from '@molstar/seqx-core/data';
import { resolveSpec } from '@molstar/seqx-core/resolve';
import { resolveViews } from '@molstar/seqx-core/view';
import {
  attachSources,
  bundledTextMetrics,
  type DrawLists,
  type Frame,
  type FrameInput,
  hitView,
  LightTokens,
  type PassName,
  PathCommandPoints,
  Renderer,
  renderConfig,
  walkPass,
} from '../index.js';
import { BuiltInLayerTypes, BuiltInLayouts } from '../built-ins.js';

/** The layout of the layer `index` of track `track` of section `section` of a frame. */
export function layerLay<T = any>(f: Frame, track: number, index = 0, section = 0): T {
  return f.layout.sections[section].tracks[track].layers[index].lay as T;
}

export interface ReadPath {
  /** The points, `[x, y, …]`. */
  readonly xy: number[];
  readonly closed: boolean;
  readonly fill: number;
  readonly stroke: number;
  readonly strokeWidth: number;
  readonly dash: number;
}

/** The paths of a pass, in order, with their points. */
export function readPaths(lists: DrawLists, pass: PassName = 'data'): ReadPath[] {
  const l = lists.passes[pass].paths;
  const out: ReadPath[] = [];
  for (let i = 0; i < l.n; i++) {
    const xy: number[] = [];
    let at = l.pointStart[i];
    let closed = false;
    for (let c = l.cmdStart[i]; c < l.cmdStart[i] + l.cmdCount[i]; c++) {
      const cmd = l.commands[c];
      if (cmd === 3) closed = true;
      for (let k = 0; k < PathCommandPoints[cmd]; k++, at++) xy.push(l.points[2 * at], l.points[2 * at + 1]);
    }
    out.push({ xy, closed, fill: l.fill[i], stroke: l.stroke[i], strokeWidth: l.strokeWidth[i], dash: l.dash[i] });
  }
  return out;
}

export interface ReadRect {
  readonly x: number;
  readonly y: number;
  readonly w: number;
  readonly h: number;
  readonly fill: number;
  readonly stroke: number;
}

/** The rects of a pass, in order. */
export function readRects(lists: DrawLists, pass: PassName = 'data'): ReadRect[] {
  const l = lists.passes[pass].rects;
  const out: ReadRect[] = [];
  for (let i = 0; i < l.n; i++)
    out.push({ x: l.x[i], y: l.y[i], w: l.w[i], h: l.h[i], fill: l.fill[i], stroke: l.stroke[i] });
  return out;
}

/** The hit of the layout at a point (the layer's track and x, y within it). */
export function hitAt(f: Frame, x: number, y: number) {
  return hitView(f.layout, x, y);
}

/** A seeded generator (a 32-bit LCG). */
export function rng(seed: number): () => number {
  let s = seed >>> 0;
  return () => {
    s = (Math.imul(s, 1664525) + 1013904223) >>> 0;
    return s / 0x100000000;
  };
}

export type { Node };

/** A settled frame of a spec drawn by a target without the layer types named (their fallbacks are drawn instead). */
export async function frameWithout(
  spec: Node,
  drop: readonly string[],
  input: FrameInput,
  runtimes: readonly FormatRuntime[] = [],
) {
  const model = resolveSpec(spec, { runtimes });
  const view = resolveViews(model)[0];
  const store = new DataStore();
  await attachSources(model, store);
  const types = Object.fromEntries(Object.entries(BuiltInLayerTypes).filter(([name]) => !drop.includes(name)));
  const renderer = new Renderer({
    layouts: BuiltInLayouts,
    model,
    view,
    data: store,
    tokens: LightTokens,
    text: bundledTextMetrics(),
    config: renderConfig(),
    types,
  });
  let f = renderer.frame(input);
  for (let i = 0; i < 4; i++) {
    store.want(f.wants);
    await store.whenSettled();
    const next = renderer.frame(input);
    if (!next.rebuilt) break;
    f = next;
  }
  return renderer.frame(input);
}

/** The labels-pass texts and rects whose visible part (under the clips around it) leaves `[top, bottom]`. */
export function outside(lists: DrawLists, pass: PassName, top: number, bottom: number): { drawn: number; bad: number } {
  const stack: [number, number][] = [];
  let bad = 0;
  let drawn = 0;
  const check = (y0: number, y1: number) => {
    drawn++;
    let lo = Number.NEGATIVE_INFINITY;
    let hi = Number.POSITIVE_INFINITY;
    for (const [a, b] of stack) {
      lo = Math.max(lo, a);
      hi = Math.min(hi, b);
    }
    if (Math.max(y0, lo) < top - 0.01 || Math.min(y1, hi) > bottom + 0.01) bad++;
  };
  walkPass(lists.passes[pass], {
    rect: (l, i) => check(l.y[i], l.y[i] + l.h[i]),
    path: () => undefined,
    glyphs: () => undefined,
    text: (l, i) => check(l.y[i] - 10, l.y[i] + 3),
    image: () => undefined,
    clip: (l, i) => stack.push([l.y[i], l.y[i] + l.h[i]]),
    unclip: () => stack.pop(),
  });
  return { drawn, bad };
}
