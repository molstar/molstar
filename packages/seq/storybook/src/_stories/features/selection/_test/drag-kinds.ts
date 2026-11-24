/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The drag kinds the drag rule is held to (the viewer's _spec, "SelectionModel and HighlightModel", Drags), each a view
// as its story composes it and the pointer's path: a rect on the 4,000 × 4,000 dot plot, a range on its ruler, a box on
// the PAE grid (CDK2, 298 × 298), a sweep across 12 chromosomes of the 3,000 contigs on the ring and on the row.

import type { GridLay, RingPlacement } from '@molstar/seqx-core/render';
import { ComplexViews, type ExampleView, TrnaScaleView } from '@molstar/seqx-examples';
import type { SeqViewer } from '@molstar/seqx-viewer/viewer';
import { vi } from 'vitest';
import { openExample } from '../../../examples/_test/draw-checks.js';
import { GenomeViews } from '../../genome/views.js';

export type Point = readonly [number, number];

export interface DragKind {
  readonly name: string;
  readonly view: ExampleView;
  /** The pointer's path on the opened view, its first point the press; `n` points. */
  path(v: SeqViewer, n: number): Point[];
}

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

const pae = ComplexViews.find((x) => x.id === 'pae-monomer')!;
const contigs = GenomeViews.find((x) => x.id === 'contigs')!;

/** The genome section's units of chr11 and of chr22, a megabase into each. */
function chromosomes(v: SeqViewer): [number, number] {
  const blocks = v.view.resolved.value!.sections.find((x) => x.address === '/views/main/genome')!.blocks!;
  return [(blocks[10].offset ?? 0) + 1e6, (blocks[21].offset ?? 0) + 1e6];
}

export const DragKinds: readonly DragKind[] = [
  {
    name: 'a rect on the 4,000 × 4,000 dot plot',
    view: TrnaScaleView,
    path: (v, n) => {
      const g = gridOf(v);
      return Array.from({ length: n }, (_, k) => cellCentre(g, (100 + (k * 2560) / n) | 0, (100 + (k * 2560) / n) | 0));
    },
  },
  {
    name: 'a range on the dot plot’s ruler',
    view: TrnaScaleView,
    path: (v, n) => {
      const ruler = v.frame!.layout.sections[0];
      return Array.from({ length: n }, (_, k): Point => [v.dataX + 20 + (k * 640) / n, ruler.y + 6]);
    },
  },
  {
    name: 'a box on the PAE grid',
    view: pae,
    path: (v, n) => {
      const g = gridOf(v);
      return Array.from({ length: n }, (_, k) => cellCentre(g, (10 + (k * 256) / n) | 0, (10 + (k * 256) / n) | 0));
    },
  },
  {
    name: 'a sweep across 12 chromosomes on the 3,000-contig ring',
    view: contigs,
    path: (v, n) => {
      const [from, to] = chromosomes(v);
      const ring = v.frame!.layout.sections.find((s) => s.section.address === '/views/main/ring')!;
      const p = ring.placement as RingPlacement;
      const y = p.strips.bandY[0] + p.strips.bandH[0] / 2;
      return Array.from({ length: n }, (_, k): Point => {
        p.warp.map(p.frameX(from + ((to - from) * k) / (n - 1)), y);
        return [p.warp.px, p.warp.py];
      });
    },
  },
  {
    name: 'a sweep across 12 chromosomes on the 3,000-contig row',
    view: contigs,
    path: (v, n) => {
      const [from, to] = chromosomes(v);
      const row = v.frame!.layout.sections.find((s) => s.section.address === '/views/main/genome')!;
      const w = v.view.get('/views/main/genome')!;
      const xOf = (u: number) => v.dataX + ((u - w.start) / (w.end - w.start)) * w.width;
      return Array.from(
        { length: n },
        (_, k): Point => [xOf(from + ((to - from) * k) / (n - 1)), row.y + row.height / 2],
      );
    },
  },
];

/** Opens a drag kind's view on the fake clock (`realClock`: the timings read `performance`), 1,440 px wide. */
export async function openDragView(view: ExampleView, realClock = false): Promise<SeqViewer> {
  const v = await openExample(view, { realClock, delayMs: 0, wait: 3000, fixtures: true });
  v.resize(1440);
  await vi.advanceTimersByTimeAsync(500);
  v.draw();
  return v;
}
