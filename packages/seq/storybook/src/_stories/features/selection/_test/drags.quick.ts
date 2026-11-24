/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A drag-select's work a move and a frame, a drag kind each (the viewer's _spec, "SelectionModel and HighlightModel",
// Drags): a rect on the 4,000 × 4,000 dot plot, a range on its ruler, a box on the PAE grid, a sweep across 12 chromosomes on the 3,000-contig ring
// and on its row. Each reports the work a pointer move does (ms, bytes, selection writes), a move and its frame, the
// layouts, and the SVG twin's export (what a story beside the canvas draws again). By hand, never in CI:
// `pnpm -s -C packages/seq/storybook test:perf`.

import type { GridLay, RingPlacement } from '@molstar/seqx-core/render';
import { loadGenome } from '@molstar/seqx-core/resolve';
import { ComplexViews, TrnaScaleView } from '@molstar/seqx-examples';
import type { SeqViewer } from '@molstar/seqx-viewer/viewer';
import { measureBytesPerStep } from '@molstar/test-harness/quick';
import { afterEach, beforeAll, describe, it, vi } from 'vitest';
import { openExample } from '../../../examples/_test/draw-checks.js';
import { GenomeViews } from '../../genome/views.js';

type Point = readonly [number, number];

afterEach(() => {
  vi.useRealTimers();
});

beforeAll(async () => {
  await loadGenome();
});

/** A grid section's lay (the first layer with two axes). */
function gridOf(v: SeqViewer): GridLay {
  for (const s of v.frame!.layout.sections)
    for (const t of s.lays) for (const l of t) if ((l.lay as GridLay | undefined)?.axisY) return l.lay as GridLay;
  throw new Error('no grid');
}

/** The page centre of a grid's cell. */
function cellCentre(g: GridLay, i: number, j: number): Point {
  g.axis!.cell(i);
  const x = (g.axis!.sx0 + g.axis!.sx1) / 2;
  g.axisY!.cell(j);
  return [x, (g.axisY!.sx0 + g.axisY!.sx1) / 2];
}

/** Opens a view for the profile, wider than the stories' default (the page the author dragged on). */
async function open(view: Parameters<typeof openExample>[0]): Promise<SeqViewer> {
  const v = await openExample(view, { realClock: true, delayMs: 0, wait: 3000, fixtures: true });
  v.resize(1440);
  await vi.advanceTimersByTimeAsync(500);
  v.draw();
  return v;
}

/** The mean ms of `n` calls of `f`. */
function msOf(f: (i: number) => void, n: number): number {
  const t0 = performance.now();
  for (let i = 0; i < n; i++) f(i);
  return (performance.now() - t0) / n;
}

/**
 * A drag along `points` (its first the press): what a move does alone (several a frame, as a pointer sends them), and
 * a move with its frame; frames slower than 20 ms are timed only (their bytes would take minutes to read).
 */
async function profile(name: string, v: SeqViewer, points: readonly Point[]): Promise<void> {
  let writes = 0;
  const off = v.selection.value.listen(() => writes++);
  const n = points.length;
  v.input.down(...points[0]);
  v.input.move(...points[1]);
  v.draw();
  const layouts = v.counts.layouts;
  const w0 = writes;
  const moveMs = msOf((i) => v.input.move(...points[i]), n);
  const perMove = (writes - w0) / n;
  const frame = (i: number) => {
    v.input.move(...points[i % n]);
    v.draw();
  };
  const frameMs = msOf(frame, n);
  const moveBytes = measureBytesPerStep((i) => v.input.move(...points[i % n])).bytesPerStep;
  const frameBytes = frameMs < 20 ? `${measureBytesPerStep(frame).bytesPerStep.toFixed(0)} B` : '(not read)';
  const laid = v.counts.layouts - layouts;
  v.input.up();
  off();
  const t = performance.now();
  const svg = v.renderSvg({ idPrefix: 'd-' });
  await vi.advanceTimersByTimeAsync(400);
  await svg;
  const svgMs = performance.now() - t;
  console.info(
    `${name}: a move ${moveMs.toFixed(3)} ms, ${moveBytes.toFixed(0)} B, ${perMove} selection writes; a move and its frame ${frameMs.toFixed(2)} ms, ${frameBytes}; ${laid} layouts; the SVG twin's export ${svgMs.toFixed(0)} ms`,
  );
}

describe('a drag-select’s work a move and a frame', () => {
  it('a rect on the 4,000 × 4,000 dot plot, and a range on its ruler', async () => {
    const v = await open(TrnaScaleView);
    const g = gridOf(v);
    const rect: Point[] = [];
    for (let k = 0; k < 64; k++) rect.push(cellCentre(g, 100 + k * 40, 100 + k * 40));
    await profile('dot plot rect', v, rect);
    const ruler = v.frame!.layout.sections[0];
    const range: Point[] = [];
    for (let k = 0; k < 64; k++) range.push([v.dataX + 20 + k * 10, ruler.y + 6]);
    await profile('dot plot ruler range', v, range);
    v.retire();
  });

  it('a box on the PAE grid (CDK2, 298 × 298)', async () => {
    const v = await open(ComplexViews.find((x) => x.id === 'pae-monomer')!);
    const g = gridOf(v);
    const box: Point[] = [];
    for (let k = 0; k < 64; k++) box.push(cellCentre(g, 10 + k * 4, 10 + k * 4));
    await profile('PAE box', v, box);
    v.retire();
  });

  it('a sweep across 12 chromosomes of the 3,000 contigs, on the ring and on the row', async () => {
    const v = await open(GenomeViews.find((x) => x.id === 'contigs')!);
    const blocks = v.view.resolved.value!.sections.find((x) => x.address === '/views/main/genome')!.blocks!;
    // chr11 into chr22, a megabase into each
    const from = (blocks[10].offset ?? 0) + 1e6;
    const to = (blocks[21].offset ?? 0) + 1e6;
    const ring = v.frame!.layout.sections.find((s) => s.section.address === '/views/main/ring')!;
    const p = ring.placement as RingPlacement;
    const y = p.strips.bandY[0] + p.strips.bandH[0] / 2;
    const round: Point[] = [];
    for (let k = 0; k < 64; k++) {
      p.warp.map(p.frameX(from + ((to - from) * k) / 63), y);
      round.push([p.warp.px, p.warp.py]);
    }
    await profile('ring sweep (12 chromosomes)', v, round);
    v.selection.clear();
    const row = v.frame!.layout.sections.find((s) => s.section.address === '/views/main/genome')!;
    const w = v.view.get('/views/main/genome')!;
    const xOf = (u: number) => v.dataX + ((u - w.start) / (w.end - w.start)) * w.width;
    const along: Point[] = [];
    for (let k = 0; k < 64; k++) along.push([xOf(from + ((to - from) * k) / 63), row.y + row.height / 2]);
    await profile('row sweep (12 chromosomes)', v, along);
    v.retire();
  });
});
