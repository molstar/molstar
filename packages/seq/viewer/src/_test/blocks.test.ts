/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A section of blocks on several systems in the viewer (_spec/README.md "ViewModel", Blocks): its axis is the
// section's own units on no system, so a drag-select across blocks selects each block's shown positions on its own
// system (never the section's units as a system's), a double-click focuses the block under it (its offset in the
// section), and the cursor crosses into the next block's first position (EC-154, amended).

import { createSeqSpec } from '@molstar/seqx-core/lang';
import { rangesLocation } from '@molstar/seqx-core/query';
import type { BlockPlacement } from '@molstar/seqx-core/render';
import { afterEach, describe, expect, it } from 'vitest';
import { DefaultSeqViewerSpec } from '../default-spec.js';
import { createSeqViewer, type SeqViewer } from '../viewer.js';
import { installFrames, type Frames } from './frames.js';

let frames: Frames | undefined;
const viewers: SeqViewer[] = [];

afterEach(() => {
  for (const v of viewers.splice(0)) v.retire();
  frames?.restore();
  frames = undefined;
});

const Section = '/views/v/complex';

/** Two chains as blocks of one section: CDK2 (120 residues, offsets 0–120), cyclin (80, offsets 120–200). */
function blocksSpec() {
  const b = createSeqSpec({ title: 'blocks' });
  const a = b.data.sequence('cdk2', { alphabet: 'protein', letters: 'ACDEFGHIKLMNPQRSTVWY'.repeat(6) });
  const c = b.data.sequence('cyclin', { alphabet: 'protein', letters: 'YWVTSRQPNMLKIHGFEDCA'.repeat(4) });
  b.views.view('v', (v) => {
    const s = v.section('complex', {});
    s.block('a', { of: a, label: 'CDK2' });
    s.block('b', { of: c, label: 'cyclin A' });
    s.track('letters').layer('letters', { type: 'letters' });
  });
  return b.build();
}

async function make() {
  frames = installFrames();
  const viewer = createSeqViewer({ ...DefaultSeqViewerSpec(), ui: undefined });
  viewers.push(viewer);
  await viewer.load(blocksSpec());
  viewer.resize(800);
  await viewer.settle();
  frames.flush();
  const p = viewer.frame!.layout.sections[0].placement as BlockPlacement;
  const s = p.strips;
  /** The page point of a block's position, on its letters. */
  const at = (block: number, position: number) => ({
    x: s.axes[block].xOfUnit(position + 0.5),
    y: s.bandY[block * s.bands] + s.bandH[block * s.bands] / 2,
  });
  return { viewer, at };
}

describe('a section of blocks on several systems in the viewer', () => {
  it('its axis is its own units, on no system', async () => {
    const { viewer } = await make();
    const axis = viewer.view.resolved.value!.sections[0].axis;
    expect(axis.system).toBe(Section);
    expect(axis.ranges).toEqual([0, 200]);
  });

  it('a drag-select across blocks selects each block’s shown positions on its own system, never the section’s units', async () => {
    const { viewer, at } = await make();
    viewer.input.mac = false;
    const from = at(0, 100);
    const to = at(1, 30);
    expect(viewer.hit(from.x, from.y)).toMatchObject({ system: '/data/cdk2', position: 100 });
    expect(viewer.input.down(from.x, from.y)).toBe(true);
    viewer.input.move(to.x, to.y);
    viewer.input.up();
    const parts = viewer.selection.value.value.shown.parts;
    // CDK2 from the anchor to its end, cyclin from its start to the pointer: each on its own system
    expect(parts.map((p) => [p.system, Array.from(p.ranges)])).toEqual([
      ['/data/cdk2', [100, 120]],
      ['/data/cyclin', [0, 31]],
    ]);
  });

  it('a double-click focuses the block under it, at its offset in the section', async () => {
    const { viewer, at } = await make();
    const p = at(1, 30);
    expect(viewer.hit(p.x, p.y)).toMatchObject({ system: '/data/cyclin', position: 30, block: 1 });
    expect(viewer.input.doubleClick(p.x, p.y)).toBe(true);
    const v = viewer.view.get(Section)!;
    expect(v.start).toBeLessThanOrEqual(150);
    expect(v.end).toBeGreaterThanOrEqual(151);
    expect(v.end - v.start).toBeLessThan(200);
  });

  it('a focus on a block’s positions moves to them at the block’s offset; a bookmark names no one system', async () => {
    const { viewer } = await make();
    viewer.addresses.focusLocation(rangesLocation('/data/cyclin', [10, 20]));
    const v = viewer.view.get(Section)!;
    expect(v.start).toBeLessThanOrEqual(130);
    expect(v.end).toBeGreaterThanOrEqual(140);
    expect(v.end - v.start).toBeLessThan(200);
    expect(viewer.addresses.shownAddress()).toBeUndefined();
  });

  it('EC-154 (amended): the cursor crosses into the next block’s first position, and stays at the section’s end', async () => {
    const { viewer, at } = await make();
    const q = at(0, 119);
    viewer.input.click(q.x, q.y);
    expect(viewer.selection.cursor.value).toMatchObject({ system: '/data/cdk2', position: 119 });
    viewer.input.key({ key: 'ArrowRight' });
    expect(viewer.selection.cursor.value).toMatchObject({ system: '/data/cyclin', position: 0 });
    const p = at(1, 79);
    viewer.input.click(p.x, p.y);
    expect(viewer.selection.cursor.value).toMatchObject({ system: '/data/cyclin', position: 79 });
    viewer.input.key({ key: 'ArrowRight' });
    expect(viewer.selection.cursor.value).toMatchObject({ system: '/data/cyclin', position: 79 });
  });
});
