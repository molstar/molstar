/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { HeadersBehavior } from '@molstar/seqx-viewer/behaviors/headers';
import type { NavigatorModel } from '@molstar/seqx-viewer/navigator';
import type { SeqViewer } from '@molstar/seqx-viewer/viewer';
import { ribbonsOf, strongestStem } from './views.js';

// The genome ring stories' gestures (Features/Genome Ring, the square and the stem menu): each an opening a story runs
// once the view is drawn, and its headless twin replays (`_test/genome-ring-stories.test.ts`). They act through the
// public models only: the square's navigator model (a point found by its hit), the headers' menu, the input's keys.

export const RingSection = '/views/main/ring';

/** The square's model, sized when nothing mounted it (headless), drawn. */
export function squareOf(viewer: SeqViewer): NavigatorModel {
  const n = viewer.navigatorModel('square')!;
  if (!(n.width > 0)) n.resize(220, 220);
  viewer.draw();
  n.draw();
  return n;
}

/** The lattice index of a chromosome on the square (its block's), by its name. */
function blockIndex(viewer: SeqViewer, n: NavigatorModel, name: string): number {
  const systems = viewer.spec.model.value.systems;
  for (let x = 0; x < n.width; x += 2) {
    const h = n.hit(x, n.height - 4);
    if (h?.system && systems.get(h.system)?.name === name) return h.block!;
  }
  return -1;
}

/** A point on the square whose hit is the cell of two chromosomes (x then y), and the part there. */
export function cellPoint(viewer: SeqViewer, n: NavigatorModel, xName: string, yName: string): [number, number] {
  const i = blockIndex(viewer, n, xName);
  const j = blockIndex(viewer, n, yName);
  let sx = 0;
  let sy = 0;
  let k = 0;
  for (let y = 0; y < n.height; y += 1)
    for (let x = 0; x < n.width; x += 1) {
      const h = n.hit(x, y);
      if (h?.block !== i || h.yBlock !== j) continue;
      sx += x;
      sy += y;
      k++;
    }
  return [sx / k, sy / k];
}

/** A point of the run's box's part (`body` its middle; `end` its bottom-right corner's grip), on the diagonal. */
export function gripPoint(n: NavigatorModel, part: 'body' | 'end'): [number, number] | undefined {
  let first: number | undefined;
  let last: number | undefined;
  for (let d = 0; d < Math.min(n.width, n.height); d += 0.5) {
    const h = n.hit(d, d);
    if (h?.part !== part) continue;
    first ??= d;
    last = d;
  }
  if (first === undefined || last === undefined) return undefined;
  const at = part === 'end' ? last : (first + last) / 2;
  return [at, at];
}

/** A drag on the square from a point to another, in steps, released. */
export function dragSquare(n: NavigatorModel, from: readonly [number, number], to: readonly [number, number]): void {
  n.down(from[0], from[1]);
  for (let s = 1; s <= 6; s++) n.move(from[0] + ((to[0] - from[0]) * s) / 6, from[1] + ((to[1] - from[1]) * s) / 6);
  n.up();
  n.click(to[0], to[1]);
}

/** A click on the square (with Cmd or Ctrl: `toggle`). */
export function clickSquare(n: NavigatorModel, at: readonly [number, number], toggle = false): boolean {
  n.down(at[0], at[1]);
  n.up();
  return n.click(at[0], at[1], 1, toggle ? { metaKey: true, ctrlKey: true } : {});
}

/** The ring shows a run (chr4 to chr6), then its box dragged three chromosomes on: chr7 to chr9. */
export async function squareDragged(viewer: SeqViewer): Promise<void> {
  await viewer.execute('section.include', { address: RingSection, include: ['chr4', 'chr5', 'chr6'] });
  viewer.draw();
  const n = squareOf(viewer);
  const body = gripPoint(n, 'body');
  const to = cellPoint(viewer, n, 'chr8', 'chr8');
  if (body) dragSquare(n, body, to);
  viewer.draw();
}

/** The ring shows chr4 and chr5, then the box's corner dragged out to chr9: chr4 to chr9. */
export async function squareResized(viewer: SeqViewer): Promise<void> {
  await viewer.execute('section.include', { address: RingSection, include: ['chr4', 'chr5'] });
  viewer.draw();
  const n = squareOf(viewer);
  const end = gripPoint(n, 'end');
  const c = cellPoint(viewer, n, 'chr9', 'chr9');
  if (end) dragSquare(n, end, [c[0] + 3, c[1] + 3]);
  viewer.draw();
}

/** One click on the cell chr4 × chr7: the pair at once. */
export async function squareClicked(viewer: SeqViewer): Promise<void> {
  const n = squareOf(viewer);
  clickSquare(n, cellPoint(viewer, n, 'chr4', 'chr7'));
  viewer.draw();
}

/** chr4 alone (its cell), then Cmd+click on chr9's cell: chr4 and chr9, a scattered set outlined cell by cell. */
export async function squareScattered(viewer: SeqViewer): Promise<void> {
  const n = squareOf(viewer);
  clickSquare(n, cellPoint(viewer, n, 'chr4', 'chr4'));
  viewer.draw();
  const m = squareOf(viewer);
  clickSquare(m, cellPoint(viewer, m, 'chr9', 'chr9'), true);
  viewer.draw();
}

/** The strongest stem's menu (its partners by count, select them, show it alone), as a secondary click opens it. */
export async function stemMenuOpen(viewer: SeqViewer): Promise<boolean> {
  const headers = viewer.behavior<HeadersBehavior>('headers');
  const l = ribbonsOf(viewer);
  if (!headers || !l || l.ns === 0) return false;
  await headers.loaded();
  const { x, y } = strongestStem(l);
  return headers.context(x, y);
}

/** The keyboard cursor on the strongest stem, then `.` twice: on to its second partner by count, said. */
export function partnerKeys(viewer: SeqViewer): void {
  const l = ribbonsOf(viewer);
  if (!l || l.ns === 0) return;
  const { stem } = strongestStem(l);
  viewer.selection.setCursor({
    section: l.section,
    system: l.systems[l.stemSystem[stem]],
    position: l.stemStart[stem],
  });
  viewer.input.key({ key: '.' });
  viewer.input.key({ key: '.' });
  viewer.draw();
}
