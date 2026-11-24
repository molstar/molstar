/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { LightTokens } from '@molstar/seqx-core/render';
import type { SquareNavigator } from '@molstar/seqx-core/render/navigator-square';
import { loadGenome } from '@molstar/seqx-core/resolve';
import type { ExampleView } from '@molstar/seqx-examples';
import { CellsType } from '@molstar/seqx-ext-rondo';
import { rondoViewer } from '@molstar/seqx-ext-rondo/viewer';
import type { NavigatorModel } from '@molstar/seqx-viewer/navigator';
import type { SeqViewer } from '@molstar/seqx-viewer/viewer';
import { afterEach, beforeAll, describe, expect, it, vi } from 'vitest';
import { openExample } from '../../../examples/_test/draw-checks.js';
import { ringStoryView } from '../views.js';

// The square beside the ring draws its cells (`rondo:square`, `rondo:cells`): a block pair by its count on a square-root
// ramp, mirrored for one assembly; hover names the pair and its links; kept when the ring shows others (the square made
// again over the same data); without the extension the lattice, said.

let viewer: SeqViewer | undefined;

beforeAll(async () => {
  await loadGenome();
});

afterEach(() => {
  viewer?.retire();
  viewer = undefined;
  vi.useRealTimers();
});

async function open(view: ExampleView, rondo = true) {
  viewer?.retire();
  viewer = await openExample(view, {
    delayMs: 0,
    wait: 3000,
    fixtures: true,
    ...(rondo ? { extensions: [rondoViewer] } : {}),
  });
  viewer.draw();
  const n = viewer.navigatorModel('square')!;
  n.resize(220, 220);
  // the cells' bins land: the square laid out until it waits on nothing
  for (let i = 0; i < 20; i++) {
    n.draw();
    await vi.advanceTimersByTimeAsync(100);
    viewer.draw();
  }
  n.draw();
  return { v: viewer, n };
}

const square = (n: NavigatorModel) => (n as unknown as { renderer: SquareNavigator }).renderer;
/** The data pass's fills in a ramp colour (the cells). */
const cellFills = (n: NavigatorModel) => {
  const ramp = new Set([0, 1, 2, 3, 4].map((k) => LightTokens.colors[`seq-ramp-${k}` as 'seq-ramp-0']));
  const r = n.frame!.lists.passes.data.rects;
  let k = 0;
  for (let i = 0; i < r.n; i++) if (ramp.has(r.fill[i])) k++;
  return k;
};
const centre = (s: SquareNavigator, i: number, j: number) =>
  [(s.edgePx(i, true) + s.edgePx(i + 1, true)) / 2, (s.edgePx(j, false) + s.edgePx(j + 1, false)) / 2] as const;

describe('The square’s cells beside the ring', () => {
  it('a cell a block pair by count, the diagonal included, mirrored; hover “chr4 × chr7 · n links”', async () => {
    const { n } = await open(ringStoryView('ring'));
    const s = square(n);
    expect(n.issues).toEqual([]);
    expect(cellFills(n)).toBeGreaterThan(20);
    const counts: number[][] = [];
    for (let i = 0; i < s.blockCount; i++) {
      counts.push([]);
      for (let j = 0; j < s.blockCount; j++)
        counts[i].push(Number(/^([\d,]+)/.exec(s.cellText(i, j) ?? '0')![1].replace(/,/g, '')));
    }
    // mirrored (one assembly: unsided), some links within a chromosome
    for (let i = 0; i < s.blockCount; i++)
      for (let j = 0; j < s.blockCount; j++) expect(counts[i][j]).toBe(counts[j][i]);
    expect(counts.some((row, i) => row[i] > 0)).toBe(true);
    // hover names the pair and its links
    const i4 = [...Array(s.blockCount).keys()].find((k) => s.blockName(k) === 'chr4')!;
    const i7 = [...Array(s.blockCount).keys()].find((k) => s.blockName(k) === 'chr7')!;
    n.move(...centre(s, i4, i7));
    expect(n.readout.value).toBe(`chr4 × chr7 · ${counts[i4][i7].toLocaleString('en-US')} links`);
  });

  it('the cells stay drawn whatever the ring shows: after “chr4 and chr7 only”, the same counts, not read again', async () => {
    const { v, n } = await open(ringStoryView('ring'));
    const s0 = square(n);
    const before = s0.cellText(3, 6);
    // the square is made again over the new view (the ring's `include`, as a release of its box writes it); its cells
    // are made of the same data, through the same store: kept, never read and counted again (12k pairs a release)
    const laid = vi.spyOn(CellsType, 'layout');
    await v.execute('section.include', { address: '/views/main/ring', include: ['chr4', 'chr7'] });
    for (let i = 0; i < 20; i++) {
      v.draw();
      n.draw();
      await vi.advanceTimersByTimeAsync(100);
    }
    const s = square(n);
    expect(s).not.toBe(s0);
    expect(laid).not.toHaveBeenCalled();
    expect(s.cellText(3, 6)).toBe(before);
    expect(cellFills(n)).toBeGreaterThan(20);
  });

  it('EC-235: without the extension the square is the lattice, its picture said as the axis', async () => {
    const { v, n } = await open(ringStoryView('ring'), false);
    expect(cellFills(n)).toBe(0);
    expect(square(n).cellText(0, 0)).toBeUndefined();
    const said = v.view.resolved.value!.issues.map((i) => i.message);
    expect(said).toContain('No picture `rondo:square` here: drawn as the axis');
  });
});
