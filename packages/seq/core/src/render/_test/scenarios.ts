/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The frame scenarios the performance targets are measured on: a view, a size (1440 CSS px at dpr 2) and a pan, a zoom
// and a hover walk, each a list of inputs made up front (so a step allocates what the renderer does, not its inputs).
// Shared by the heavy time bounds (`frame-time.test.ts`) and the allocation case (`frame.quick.ts`).

import type { Node } from '@molstar/core/repr/tree';
import type { FormatRuntime } from '@molstar/seqx-core/data';
import { type Emphasis, type FrameInput, type Hit, hitView, type RenderConfigOverrides } from '../index.js';
import { BIG, BigRuntime, bigSpec, p53GeneSpec, p53TracksSpec, PIECES, piecesSpec } from './fixtures.js';
import { type Settled, settle } from './helpers.js';

export const SIZE = { width: 1440, height: 400 } as const;
export const DPR = 2;
/** Steps a walk takes: its inputs repeat after this many (a walk goes out and back). */
export const STEPS = 256;

export interface Scenario {
  readonly name: string;
  readonly spec: () => Node;
  readonly runtimes: readonly FormatRuntime[];
  /** The widest window, in display units, and its centre; a pan moves it, a zoom ranges between half of it and all. */
  readonly width: number;
  readonly center: number;
  /** A pan step, in display units. */
  readonly panStep: number;
  /** The render config's overrides (the defaults otherwise). */
  readonly config?: RenderConfigOverrides;
}

// the data box is 1304 CSS px: letters need 10 px a position (130 positions), cells 1.5 px, bars down to 1/128 px
export const Scenarios: readonly Scenario[] = [
  // a reversed axis of 150,000 single-position pieces (EC-33): a pan over 1,000 of them, their breaks in view. First:
  // what it leaves compiled must not make the others' warm draws allocate
  { name: 'reversed, 150,000 pieces', spec: piecesSpec, runtimes: [], width: 1000, center: PIECES / 2, panStep: 2 },
  { name: 'p53 tracks, letters', spec: p53TracksSpec, runtimes: [], width: 100, center: 196, panStep: 0.3 },
  { name: 'p53 tracks, whole protein', spec: p53TracksSpec, runtimes: [], width: 393, center: 196, panStep: 0.5 },
  { name: 'p53 on its gene, letters', spec: p53GeneSpec, runtimes: [], width: 100, center: 196, panStep: 0.3 },
  { name: 'p53 on its gene, whole protein', spec: p53GeneSpec, runtimes: [], width: 393, center: 196, panStep: 0.5 },
  // the gene reversed and collapsed (its coding pieces), the protein's letters a codon wide on it
  {
    name: 'p53 on its gene, coding',
    spec: () => p53GeneSpec({ coding: true }),
    runtimes: [],
    width: 60,
    center: 196,
    panStep: 0.3,
  },
  {
    name: '3M letters (120 positions)',
    spec: bigSpec,
    runtimes: [BigRuntime],
    width: 120,
    center: BIG / 2,
    panStep: 0.25,
  },
  {
    name: '3M cells (800 positions)',
    spec: bigSpec,
    runtimes: [BigRuntime],
    width: 800,
    center: BIG / 2,
    panStep: 0.5,
  },
  {
    // the default level (`values.coarser` 1): the values layer places the 100,000 positions in view one by one,
    // directly (the data is on the section's own system)
    name: '3M bars (100,000 positions)',
    spec: bigSpec,
    runtimes: [BigRuntime],
    width: 100_000,
    center: BIG / 2,
    panStep: 100,
  },
  // a level up to two device pixels' positions (`values.coarser` 2): the values layer reads the 64 level's bins
  {
    name: '3M bars (100,000 positions), level of two pixels',
    spec: bigSpec,
    runtimes: [BigRuntime],
    width: 100_000,
    center: BIG / 2,
    panStep: 100,
    config: { values: { coarser: 2 } },
  },
  // above the source's second level (64 positions) a device pixel spans more than a level: the frame reads bins
  {
    name: '3M bars (200,000 positions)',
    spec: bigSpec,
    runtimes: [BigRuntime],
    width: 200_000,
    center: BIG / 2,
    panStep: 200,
  },
];

/** A triangle wave of `STEPS` samples: 0 up to 1 and back. */
const tri = (i: number) => {
  const t = (i % STEPS) / (STEPS / 2);
  return t <= 1 ? t : 2 - t;
};

export type WalkKind = 'pan' | 'zoom';

/** The inputs of a walk: a pan slides the window a step at a time, a zoom narrows it to half and back around its centre. */
export function walk(s: Scenario, kind: WalkKind): FrameInput[] {
  const out: FrameInput[] = [];
  for (let i = 0; i < STEPS; i++) {
    let start: number;
    let end: number;
    if (kind === 'pan') {
      start = s.center - s.width / 2 + tri(i) * s.panStep * (STEPS / 2);
      end = start + s.width;
    } else {
      const w = s.width * (1 - 0.5 * tri(i));
      start = s.center - w / 2;
      end = s.center + w / 2;
    }
    out.push({ size: SIZE, dpr: DPR, viewport: { start, end } });
  }
  return out;
}

/** A scenario's renderer with the data of every window of its walks in the store. */
export async function prepare(s: Scenario): Promise<Settled> {
  const settled = await settle(s.spec(), s.runtimes, s.config);
  for (const kind of ['pan', 'zoom'] as const) {
    const inputs = walk(s, kind);
    for (let i = 0; i < STEPS; i += 16) await settled.frame(inputs[i]);
    await settled.frame(inputs[STEPS / 2]);
  }
  return settled;
}

/** The frame input a hover walk draws at: the scenario's widest window, as the viewer's frame is. */
export function hoverInput(s: Scenario): FrameInput {
  return { size: SIZE, dpr: DPR, viewport: { start: s.center - s.width / 2, end: s.center + s.width / 2 } };
}

/** Hits along every track, one per 5 CSS px of the data box, for a hover walk: the emphasis of each, `STEPS` of them. */
export function hovers(
  settled: Settled,
  input: FrameInput,
): { hits: Hit[]; points: { x: number; y: number }[]; emphasis: Emphasis[] } {
  const layout = settled.renderer.frame(input).layout;
  const hits: Hit[] = [];
  const points: { x: number; y: number }[] = [];
  for (let x = 130; x < 1430; x += 5)
    for (const sec of layout.sections)
      for (const t of sec.tracks) {
        const y = t.box.y + t.box.height / 2;
        const h = hitView(layout, x, y);
        if (h) {
          hits.push(h);
          points.push({ x, y });
        }
      }
  const emphasis: Emphasis[] = [];
  for (let i = 0; i < STEPS; i++) emphasis.push({ hover: hits[i % hits.length] });
  return { hits, points, emphasis };
}

/**
 * The hover as the viewer draws it (viewer.ts, `_updateEmphasis`): one emphasis record holding one hit record, the
 * hit's fields set to each hit in turn (`set(i)`), so a hover step is the viewer's: a field write and a frame.
 */
export function reusedHover(hits: readonly Hit[]): { emphasis: Emphasis; set(i: number): void } {
  const hit = { layer: '', section: '', system: '' } as { -readonly [K in keyof Hit]: Hit[K] };
  hit.position = undefined;
  hit.element = undefined;
  hit.extent = undefined;
  return {
    emphasis: { hover: hit },
    set(i) {
      const h = hits[i % hits.length];
      hit.layer = h.layer;
      hit.section = h.section;
      hit.system = h.system;
      hit.position = h.position;
      hit.element = h.element;
      hit.extent = h.extent;
    },
  };
}

/** Median, 95th percentile and maximum, ms. */
export function stats(times: number[]): { median: number; p95: number; max: number } {
  const t = [...times].sort((a, b) => a - b);
  return { median: t[t.length >> 1], p95: t[Math.floor(t.length * 0.95)], max: t[t.length - 1] };
}

/** Times `step` over `count` calls after `warm` unmeasured ones. */
export function time(step: (i: number) => void, count = 200, warm = 1000): number[] {
  for (let i = 0; i < warm; i++) step(i);
  const out: number[] = [];
  for (let i = 0; i < count; i++) {
    const t0 = performance.now();
    step(i);
    out.push(performance.now() - t0);
  }
  return out;
}
