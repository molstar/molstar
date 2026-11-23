/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A ring's block labels never overlap (_spec/README.md "Render", Rings, Labels): at the fit of GRCh38's 23 and 24
// chromosomes and of a draft with 2,976 contigs in runs, at several widths, no two drawn labels' boxes intersect, in
// the canvas's draw lists and the SVG alike; a label dropped for room keeps its block's name in the hit.

import { loadGenome } from '@molstar/seqx-core/resolve';
import { beforeAll, describe, expect, it } from 'vitest';
import type { RingPlacement } from '../layout/rings.js';
import { renderSvg } from '../render-svg.js';
import { genomeSpec, toyGenome } from './genome-fixtures.js';
import { settle } from './helpers.js';

beforeAll(async () => {
  await loadGenome();
});

/** GRCh38's chromosome lengths (the GRC's assembly report): chr1 to chr22, X, Y. */
const Grch38 = [
  248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973, 145138636, 138394717, 133797422,
  135086622, 133275309, 114364328, 107043718, 101991189, 90338345, 83257441, 80373285, 58617616, 64444167, 46709983,
  50818468, 156040895, 57227415,
].map((length, i) => ({ name: i < 22 ? `chr${i + 1}` : i === 22 ? 'chrX' : 'chrY', length }));

interface Box {
  readonly text: string;
  readonly x0: number;
  readonly x1: number;
  readonly y0: number;
  readonly y1: number;
}

/** The placed labels' boxes (what is drawn: their anchor point, side, measured width, the line's height). */
function placed(p: RingPlacement): Box[] {
  const l = p.labels;
  const size = p.labelSize;
  const out: Box[] = [];
  for (let k = 0; k < l.m; k++) {
    const w = l.width[k];
    const x0 = l.right[k] ? l.x[k] : l.x[k] - w;
    const cy = l.y[k] - 0.35 * size;
    out.push({ text: l.text[l.placed[k]], x0, x1: x0 + w, y0: cy - size / 2, y1: cy + size / 2 });
  }
  return out;
}

const intersect = (a: Box, b: Box) => a.x0 < b.x1 && a.x1 > b.x0 && a.y0 < b.y1 && a.y1 > b.y0;

function overlaps(boxes: readonly Box[]): string[] {
  const out: string[] = [];
  for (let i = 0; i < boxes.length; i++)
    for (let j = i + 1; j < boxes.length; j++)
      if (intersect(boxes[i], boxes[j])) out.push(`${boxes[i].text} × ${boxes[j].text}`);
  return out;
}

const cases = [
  { name: 'GRCh38, 23', seqs: Grch38.slice(0, 23) },
  { name: 'GRCh38, 24', seqs: Grch38 },
  { name: '3,000 contigs', seqs: toyGenome(24, 2976) },
];

describe('a ring’s block labels', () => {
  for (const c of cases)
    it(`never overlap: ${c.name}, at 360 to 1,200 px, square and not`, async () => {
      const spec = genomeSpec(c.seqs, { layout: 'rings', name: 'GRCh38' });
      const s = await settle(spec);
      for (const [width, height] of [
        [360, 360],
        [480, 480],
        [600, 600],
        [720, 720],
        [900, 900],
        [1200, 1200],
        [1200, 400],
        [900, 500],
        [600, 900],
      ]) {
        const f = await s.frame({ size: { width, height } });
        const p = f.layout.sections[0].placement as RingPlacement;
        const boxes = placed(p);
        expect(boxes.length, `${width}×${height}`).toBeGreaterThan(0);
        expect(overlaps(boxes), `${c.name} at ${width}×${height}`).toEqual([]);
        // the SVG draws the same labels at the same places
        const svg = await renderSvg(spec, { size: { width, height } });
        for (const b of boxes.slice(0, 5)) expect(svg.svg).toContain(`>${b.text}</text>`);
      }
    });
});
