/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The ESPript-like figure (figure-docs.ts) over the core's globin fixture read from disk (Node only), and the render
// props of a line per 60 columns.

import { readFileSync } from 'node:fs';
import type { SeqViewSpecDoc } from '@molstar/seq-core/format';
import { type AlnRow, type FigureOptions, figureFromRows, parseClustal } from './figure-docs.js';

export { type AlnRow, FIGURE_FRAMES, type FrameRow, GLOBIN_NWK, sixtyColumns } from './figure-docs.js';

/** The fixture's rows: name (`ID/first-last`), gapped residues, and each block's printed end numbers. */
export function globinRows(): AlnRow[] {
  return parseClustal(readFileSync(new URL('../../core/_test/fixtures/pf00042-seed-8.aln', import.meta.url), 'utf8'));
}

/** The figure over the fixture's rows (`FigureOptions`). */
export const figureDoc = (o: FigureOptions = {}): SeqViewSpecDoc => figureFromRows(globinRows(), o);

/**
 * Whether the SVG's items ("x,y,w,h" rects or "x0,y0,x1,y1" segments) are the recording's, all moved by one (dx, dy)
 * within a few px (an export places its section in its document; nothing else may differ).
 */
export function translated(svg: readonly string[], drawn: readonly string[], segments = false): boolean {
  if (svg.length !== drawn.length || svg.length === 0) return false;
  const key = (s: string) => s.split(',').map(Number);
  const order = (a: number[], b: number[]) => a[1] - b[1] || a[0] - b[0] || a[3] - b[3] || a[2] - b[2];
  const A = svg.map(key).sort(order);
  const B = drawn.map(key).sort(order);
  const dx = A[0][0] - B[0][0];
  const dy = A[0][1] - B[0][1];
  // A rect's size stays; a segment's other end moves with it.
  const shift = segments ? [dx, dy, dx, dy] : [dx, dy, 0, 0];
  return Math.abs(dx) <= 8 && Math.abs(dy) <= 8 && A.every((a, i) => a.every((v, j) => v === B[i][j] + shift[j]));
}

/** An SVG's rects of a fill (`<path fill="…" d="Mx yhWvHh-Wz…">`, the SVG backend's batches) as "x,y,w,h". */
export function svgRects(xml: string, fill: string): string[] {
  const out: string[] = [];
  for (const m of xml.matchAll(new RegExp(`<path fill="${fill}"[^>]* d="([^"]*)"`, 'gi')))
    for (const r of m[1].matchAll(/M(-?[\d.]+) (-?[\d.]+)h(-?[\d.]+)v(-?[\d.]+)h-?[\d.]+z/g))
      out.push(`${r[1]},${r[2]},${r[3]},${r[4]}`);
  return out;
}
/** An SVG's line segments of a stroke (`<path stroke="…" … d="Mx yLx y…">`, the SVG backend's) as "x0,y0,x1,y1". */
export function svgSegments(xml: string, stroke: string): string[] {
  const out: string[] = [];
  for (const m of xml.matchAll(new RegExp(`<path stroke="${stroke}"[^>]* d="([^"]*)"`, 'gi')))
    for (const r of m[1].matchAll(/M(-?[\d.]+) (-?[\d.]+)L(-?[\d.]+) (-?[\d.]+)/g))
      out.push(`${r[1]},${r[2]},${r[3]},${r[4]}`);
  return out;
}
