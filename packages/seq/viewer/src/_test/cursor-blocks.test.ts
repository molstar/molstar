/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The keyboard cursor across the gaps between blocks on a grid over two sections of blocks (_spec/README.md
// "Behaviors", Cursor; EC-154, amended): ← → and ↑ ↓ cross from a block's last position to the next block's first, and
// stay at the ends, as the cursor of a row of blocks does.

import { createSeqSpec } from '@molstar/seqx-core/lang';
import { afterEach, describe, expect, it } from 'vitest';
import { DefaultSeqViewerSpec } from '../default-spec.js';
import { createSeqViewer, type SeqViewer } from '../viewer.js';
import { type Frames, installFrames } from './frames.js';

let frames: Frames | undefined;
const viewers: SeqViewer[] = [];

afterEach(() => {
  for (const v of viewers.splice(0)) v.retire();
  frames?.restore();
  frames = undefined;
});

const Plot = '/views/v/plot';

/** Blocks `s0` (5) and `s1` (4) across, `t0` (3) and `t1` (4) down, one grid over each pair. */
function spec() {
  const b = createSeqSpec({ title: 'cursor over blocks' });
  const xs = [5, 4].map((n, k) => b.data.sequence(`s${k}`, { alphabet: 'dna', letters: 'ACGT'.repeat(3).slice(0, n) }));
  const ys = [3, 4].map((n, k) => b.data.sequence(`t${k}`, { alphabet: 'dna', letters: 'TGCA'.repeat(3).slice(0, n) }));
  const grids: string[] = [];
  [5, 4].forEach((nx, x) =>
    [3, 4].forEach((ny, y) =>
      grids.push(b.data.grid(`g${x}${y}`, { x: xs[x], y: ys[y], values: new Float32Array(nx * ny) })),
    ),
  );
  b.views.view('v', (v) => {
    const a = v.section('a', {});
    xs.forEach((q, k) => a.block(`b${k}`, { of: q }));
    a.track('letters').layer('letters', { type: 'letters' });
    const bs = v.section('b', {});
    ys.forEach((q, k) => bs.block(`b${k}`, { of: q }));
    bs.track('letters').layer('letters', { type: 'letters' });
    const g = v.section('plot', { layout: 'grid', follow: a.address, follow_y: bs.address });
    g.track('plot', {}, (t) => grids.forEach((d, k) => t.layer(`g${k}`, { type: 'grid', data: d })));
  });
  return b.build();
}

async function make() {
  frames = installFrames();
  const viewer = createSeqViewer({ ...DefaultSeqViewerSpec(), ui: undefined });
  viewers.push(viewer);
  await viewer.load(spec());
  viewer.resize(800, 800);
  await viewer.settle();
  frames.flush();
  viewer.draw();
  return viewer;
}

describe('the cursor on a grid over blocks', () => {
  it('EC-154: ← → cross from a block’s last position to the next block’s first, and stay at the ends', async () => {
    const v = await make();
    const c = { section: Plot, system: '/data/s0', position: 4, ySystem: '/data/t0', yPosition: 1 };
    v.selection.setCursor(c);
    v.runCommand('cursor.next', {});
    expect(v.selection.cursor.value).toMatchObject({
      system: '/data/s1',
      position: 0,
      ySystem: '/data/t0',
      yPosition: 1,
    });
    v.runCommand('cursor.prev', {});
    expect(v.selection.cursor.value).toMatchObject({ system: '/data/s0', position: 4 });
    v.selection.setCursor({ ...c, system: '/data/s1', position: 3 });
    v.runCommand('cursor.next', {});
    expect(v.selection.cursor.value).toMatchObject({ system: '/data/s1', position: 3 });
  });

  it('EC-154: ↑ ↓ cross from a block’s last row to the next block’s first, and stay at the ends', async () => {
    const v = await make();
    const c = { section: Plot, system: '/data/s0', position: 2, ySystem: '/data/t0', yPosition: 2 };
    v.selection.setCursor(c);
    v.runCommand('cursor.down', {});
    expect(v.selection.cursor.value).toMatchObject({
      system: '/data/s0',
      position: 2,
      ySystem: '/data/t1',
      yPosition: 0,
    });
    v.runCommand('cursor.up', {});
    expect(v.selection.cursor.value).toMatchObject({ ySystem: '/data/t0', yPosition: 2 });
    v.selection.setCursor({ ...c, ySystem: '/data/t1', yPosition: 3 });
    v.runCommand('cursor.down', {});
    expect(v.selection.cursor.value).toMatchObject({ ySystem: '/data/t1', yPosition: 3 });
  });

  it('EC-154: ↓ from the section above enters the first block’s first row', async () => {
    const v = await make();
    v.selection.setCursor({ section: '/views/v/b', system: '/data/t0', position: 1 });
    v.runCommand('cursor.up', {});
    v.runCommand('cursor.down', {});
    v.runCommand('cursor.down', {});
    expect(v.selection.cursor.value).toMatchObject({ section: Plot });
    expect(v.selection.cursor.value?.ySystem).toBe('/data/t0');
    expect(v.selection.cursor.value?.yPosition).toBe(0);
  });
});
