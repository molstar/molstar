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
import {
  BIG,
  BigRuntime,
  bigSpec,
  holeySpec,
  p53GeneSpec,
  p53TextSpec,
  p53TracksSpec,
  PIECES,
  piecesSpec,
  smallLayersSpec,
} from './fixtures.js';
import { GffRuntime } from '../../runtimes/gff.js';
import { geneModelsSpec } from './gene-fixtures.js';
import { type Settled, settle } from './helpers.js';
import { moleculesSpec } from './molecules-fixtures.js';

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
  /** The bands' own scroll by track address (a track whose lanes scroll inside it). */
  readonly bands?: Readonly<Record<string, { readonly bandTop: number }>>;
}

// the data box is 1304 CSS px: letters need 10 px a position (130 positions), cells 1.5 px, bars down to 1/128 px
export const Scenarios: readonly Scenario[] = [
  // a reversed axis of 150,000 single-position pieces (EC-33): a pan over 1,000 of them, their breaks in view. First:
  // what it leaves compiled must not make the others' warm draws allocate
  { name: 'reversed, 150,000 pieces', spec: piecesSpec, runtimes: [], width: 1000, center: PIECES / 2, panStep: 2 },
  // a light chain's rulers over a numbering that stops at 107 of 207: the unnumbered tail marked, and labelled by ordinals
  { name: 'holey rulers', spec: holeySpec, runtimes: [], width: 207, center: 103.5, panStep: 0.5 },
  { name: 'p53 tracks, letters', spec: p53TracksSpec, runtimes: [], width: 100, center: 196, panStep: 0.3 },
  { name: 'p53 tracks, whole protein', spec: p53TracksSpec, runtimes: [], width: 393, center: 196, panStep: 0.5 },
  // the small layers: a ruler, secondary structure, variants, a line, an area, ranges dashed with arrows and highlighted
  {
    name: 'p53 small layers, whole protein',
    spec: smallLayersSpec,
    runtimes: [],
    width: 393,
    center: 196,
    panStep: 0.5,
  },
  { name: 'p53 small layers, letters', spec: smallLayersSpec, runtimes: [], width: 60, center: 175, panStep: 0.2 },
  { name: 'p53 on its gene, letters', spec: p53GeneSpec, runtimes: [], width: 100, center: 196, panStep: 0.3 },
  // the p53 text in a box of 80 px: a walk's start is the anchor, so a pan scrolls its lines (strips through their
  // translation warp, the lines laid out lazily)
  {
    name: 'p53 text',
    spec: () => p53TextSpec(),
    runtimes: [],
    width: 40,
    center: 100,
    panStep: 0.5,
    config: { text: { maxHeight: 80 } },
  },
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
  // a complex: 200 copies of a chain merged into one block (EC-67), a glycan, ATP and an ion as chips, the modified
  // residues marked, the links as arcs; the window over the block and the chips
  {
    name: 'complex, 200 copies',
    spec: () => moleculesSpec({ chains: 200 }),
    runtimes: [],
    width: 91,
    center: 45.5,
    panStep: 0.25,
  },
  // gene models (the gene line): a generated GFF3's transcripts over 2 Mb, a window of 60 kb (a few genes, their
  // boxes, hairlines and chevrons; below a position a pixel, never merged)
  {
    name: 'gene models (60 kb)',
    spec: geneModelsSpec,
    runtimes: [GffRuntime],
    width: 60_000,
    center: 400_000,
    panStep: 60,
  },
  // the same models in a track whose lanes scroll inside it (the track's `height` 120 holds eight lanes): a band's top
  // lane in the middle of the lanes, the lanes out of the band no marks
  {
    name: 'gene models, scrolled (60 kb)',
    spec: geneModelsSpec,
    runtimes: [GffRuntime],
    width: 60_000,
    center: 400_000,
    panStep: 60,
    bands: { '/views/main/chr1/tx': { bandTop: 2.5 } },
  },
  // and scrolled far (another band top: other lanes in view, the cut edge at its top)
  {
    name: 'gene models, scrolled far (60 kb)',
    spec: geneModelsSpec,
    runtimes: [GffRuntime],
    width: 60_000,
    center: 400_000,
    panStep: 60,
    bands: { '/views/main/chr1/tx': { bandTop: 20 } },
  },
  // the lanes rail dragged: the band's top moves a fraction of a lane a frame (the stored `bandTop` a drag writes), the rail
  // drawn with its window where it is
  {
    name: 'gene models, rail dragged (60 kb)',
    spec: geneModelsSpec,
    runtimes: [GffRuntime],
    width: 60_000,
    center: 400_000,
    panStep: 60,
    bands: { '/views/main/chr1/tx': { bandTop: 11.25 } },
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
  // dense items with arrowheads, at the step before bars: 11 to 23 px wide, each label measured too wide for its mark
  // before the early skip (range-labels.ts, the labels inside marks)
  {
    name: '3M items with arrowheads (8,000 positions)',
    spec: () => bigSpec({ arrow: true }),
    runtimes: [BigRuntime],
    width: 8000,
    center: BIG / 2,
    panStep: 8,
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
  // the values layer's looks (VL step 1, internal): the whole protein as a profile, and 40 residues as bars
  ...(['profile', 'bars', 'quiet'] as const).flatMap((look) => [
    {
      name: `p53 tracks, whole protein, values look ${look}`,
      spec: p53TracksSpec,
      runtimes: [],
      width: 393,
      center: 196,
      panStep: 0.5,
      config: { values: { look } },
    },
    {
      name: `p53 tracks, 40 residues, values look ${look}`,
      spec: p53TracksSpec,
      runtimes: [],
      width: 40,
      center: 196,
      panStep: 0.3,
      config: { values: { look } },
    },
  ]),
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
    out.push({ size: SIZE, dpr: DPR, viewport: { start, end, ...(s.bands ? { bands: s.bands } : {}) } });
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
  return {
    size: SIZE,
    dpr: DPR,
    viewport: { start: s.center - s.width / 2, end: s.center + s.width / 2, ...(s.bands ? { bands: s.bands } : {}) },
  };
}

/** Hits along every track, one per 5 CSS px of the data box, for a hover walk: the emphasis of each, `STEPS` of them. */
export function hovers(
  settled: Settled,
  input: FrameInput,
): { hits: Hit[]; points: { x: number; y: number }[]; emphasis: Emphasis[] } {
  const layout = settled.renderer.frame(input).layout;
  const hits: Hit[] = [];
  const points: { x: number; y: number }[] = [];
  // a section of several strips (a text's lines): along each band of each line
  for (const sec of layout.sections) {
    const st = sec.placement.strips;
    if (st.n < 2) continue;
    for (let i = 0; i < st.n; i++)
      for (let j = 0; j < st.bands; j++) {
        const k = i * st.bands + j;
        if (st.bandH[k] <= 0) continue;
        const a = st.axes[i];
        const y = st.dy[i] + st.bandY[k] + st.bandH[k] / 2;
        for (let x = a.x0 + 2; x < a.x0 + a.width; x += 5) {
          const h = hitView(layout, x + st.dx[i], y);
          if (h) {
            hits.push(h);
            points.push({ x: x + st.dx[i], y });
          }
        }
      }
  }
  for (let x = 130; x < 1430; x += 5)
    for (const sec of layout.sections)
      for (const t of sec.placement.strips.n < 2 ? sec.tracks : []) {
        // near the track's top, else its middle: a tall band of lanes starts them at its top (core _spec, Ranges)
        for (const y of [t.box.y + Math.min(t.box.height / 2, 8), t.box.y + t.box.height / 2]) {
          const h = hitView(layout, x, y);
          if (h) {
            hits.push(h);
            points.push({ x, y });
            break;
          }
        }
      }
  const emphasis: Emphasis[] = [];
  for (let i = 0; i < STEPS; i++) emphasis.push({ hover: hits[i % hits.length] });
  return { hits, points, emphasis };
}

/** A hit as one shape: every field, in one order (what the viewer's HoverModel record is). */
function oneShape(h: Hit): Hit {
  return {
    layer: h.layer,
    section: h.section,
    system: h.system,
    position: h.position,
    element: h.element,
    extent: h.extent,
    block: h.block,
    ySystem: h.ySystem,
    yPosition: h.yPosition,
    target: h.target,
    state: h.state,
  };
}

/**
 * The hover as the viewer draws it (viewer-emphasis.ts, `hoverHit`): one emphasis record holding one hit record, the
 * hit's fields set to each hit in turn (`set(i)`), so a hover step is the viewer's: a field write and a frame. The
 * fields are read from records of one shape, as the viewer reads them from its HoverModel's one record. Why: read from
 * the layers' hits (a shape per layer type) its loads went megamorphic once the complex's shapes had joined (the
 * feedback is the one closure site's, every scenario's; `--log-ic`), and `set` allocated 16 B on two hits in three (the
 * heap sampler): every scenario after the complex's measured 6 to 11 B a hover. A finding of the test, not the renderer.
 */
export function reusedHover(source: readonly Hit[]): { emphasis: Emphasis; set(i: number): void } {
  const hits = source.map(oneShape);
  const hit = { layer: '', section: '', system: '' } as { -readonly [K in keyof Hit]: Hit[K] };
  hit.position = undefined;
  hit.element = undefined;
  hit.extent = undefined;
  hit.block = undefined;
  hit.ySystem = undefined;
  hit.yPosition = undefined;
  hit.target = undefined;
  hit.state = undefined;
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
      hit.block = h.block;
      hit.ySystem = h.ySystem;
      hit.yPosition = h.yPosition;
      hit.target = h.target;
      // a rows column cell's name (the tip reads the cell it names)
      hit.state = h.state;
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
