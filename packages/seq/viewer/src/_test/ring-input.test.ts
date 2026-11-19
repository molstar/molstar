/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A ring in the viewer (_spec/README.md "Input", Rings), headless: input read through the ring's placement (its `at`
// and `offsetAt` by angle): a zoom about the pointer's angle over the origin, Space+drag moving the window with the
// pointer round the ring, a drag selecting an arc (over the origin: two runs), a double-click focusing an item over the
// origin by the shorter arc, the window's grips dragged by the `window` behavior, a history step restored by position
// across a layout patch (tracks to rings and back).

import type { Node } from '@molstar/core/repr/tree';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import type { RingPlacement } from '@molstar/seqx-core/render';
import { afterEach, describe, expect, it, vi } from 'vitest';
import { viewerItem } from '../config.js';
import { DefaultSeqViewerSpec } from '../default-spec.js';
import { createSeqViewer, type SeqViewer } from '../viewer.js';
import { type Frames, installFrames } from './frames.js';

const Ring = '/views/v/ring';
const Detail = '/views/v/detail';
const Length = 600;

function spec(): Node {
  const b = createSeqSpec({ title: 'ring' });
  const c = b.data.sequence('c', { alphabet: 'dna', letters: 'ACGT'.repeat(Length / 4), cyclic: true });
  const f = b.data.ranges('f', {
    target: c,
    rows: { start: [560, 200], end: [40, 260], id: ['ori', 'gene'], label: ['ori', 'gene'] },
  });
  b.views.view('v', (v) => {
    const ring = v.section('ring', { of: c, layout: 'rings', labels: 'none' });
    ring.track('f').layer('f', { type: 'ranges', data: f });
    ring.track('l').layer('l', { type: 'letters' });
    v.section('detail', { of: c, follow: ring }).track('l').layer('l', { type: 'letters' });
  });
  return b.build();
}

let frames: Frames | undefined;
let viewer: SeqViewer | undefined;

afterEach(() => {
  viewer?.retire();
  viewer = undefined;
  frames?.restore();
  frames = undefined;
  vi.useRealTimers();
});

async function make() {
  vi.useFakeTimers({ toFake: ['Date', 'performance'] });
  frames = installFrames();
  const v = createSeqViewer({ ...DefaultSeqViewerSpec(), ui: undefined, config: [[viewerItem('view.jumpMs'), 0]] });
  viewer = v;
  await v.load(spec());
  v.resize(640);
  await v.settle();
  frames.flush();
  v.input.mac = false;
  return v;
}

const ringOf = (v: SeqViewer) =>
  v.frame!.layout.sections.find((s) => s.section.address === Ring)!.placement as RingPlacement;

/** The page point of a unit on the ring, `y` px into its rings (negative: on the window's rim). */
function pointAt(v: SeqViewer, unit: number, y: number) {
  const p = ringOf(v);
  p.warp.map(p.frameX(unit), y);
  return { x: p.warp.px, y: p.warp.py };
}

/** The data's middle depth of the ring's first band. */
const bandY = (v: SeqViewer, band = 0) => {
  const s = ringOf(v).strips;
  return s.bandY[band] + s.bandH[band] / 2;
};

describe('a ring in the viewer, its input through its placement', () => {
  it('a zoom keeps the unit under the pointer’s angle, over the origin too', async () => {
    const v = await make();
    v.view.setRange(Ring, 560, 640);
    frames!.flush();
    const q = pointAt(v, 10.5, bandY(v));
    expect(v.input.wheel(q.x, q.y, { deltaX: 0, deltaY: -100, ctrlKey: true, deltaMode: 0 })).toBe(true);
    frames!.flush();
    const r = v.view.get(Ring)!;
    expect(r.end - r.start).toBeLessThan(80);
    // the unit stays in the window, a turn past the total
    expect(r.start).toBeLessThan(610.5);
    expect(r.end).toBeGreaterThan(610.5);
  });

  it('Space+drag moves the window with the pointer round the ring', async () => {
    const v = await make();
    v.view.setRange(Ring, 100, 160);
    frames!.flush();
    v.input.setSpace(true);
    const a = pointAt(v, 130, bandY(v));
    const b = pointAt(v, 170, bandY(v));
    expect(v.input.down(a.x, a.y)).toBe(true);
    v.input.move(b.x, b.y);
    v.input.up();
    v.input.setSpace(false);
    const r = v.view.get(Ring)!;
    expect(r.start).toBeCloseTo(140, 0);
    expect(r.end).toBeCloseTo(200, 0);
  });

  it('a drag along the ring selects an arc: over the origin, two runs', async () => {
    const v = await make();
    const a = pointAt(v, 590.5, bandY(v, 1));
    const b = pointAt(v, 9.5, bandY(v, 1));
    expect(v.input.down(a.x, a.y)).toBe(true);
    v.input.move(b.x, b.y);
    v.input.up();
    const parts = v.selection.value.value.shown.parts;
    expect(parts.length).toBe(1);
    expect(Array.from(parts[0].ranges)).toEqual([0, 10, 590, 600]);
  });

  it('a double-click on an item over the origin focuses it by the shorter arc', async () => {
    const v = await make();
    const q = pointAt(v, 20.5, ringOf(v).strips.bandY[0] + 3);
    expect(v.hit(q.x, q.y)?.element?.key).toBe('ori');
    expect(v.input.doubleClick(q.x, q.y)).toBe(true);
    frames!.flush();
    const r = v.view.get(Ring)!;
    expect(r.start).toBeCloseTo(560, 6);
    expect(r.end).toBeCloseTo(640, 6);
  });

  it('the window’s grip on the rim is dragged by the window behavior, across the origin', async () => {
    const v = await make();
    v.view.setRange(Ring, 20, 80);
    frames!.flush();
    const a = pointAt(v, 80, -5);
    expect(v.hit(a.x, a.y)).toMatchObject({ target: 'window', part: 'end' });
    const b = pointAt(v, 100, -5);
    expect(v.input.down(a.x, a.y)).toBe(true);
    v.input.move(b.x, b.y);
    v.input.up();
    expect(v.view.get(Ring)!.end).toBeCloseTo(100, 0);
    expect(v.view.get(Ring)!.start).toBeCloseTo(20, 6);
    // the start grip past the origin, backwards
    frames!.flush();
    const c = pointAt(v, 20, -5);
    const d = pointAt(v, 590, -5);
    expect(v.input.down(c.x, c.y)).toBe(true);
    v.input.move(d.x, d.y);
    v.input.up();
    const r = v.view.get(Ring)!;
    expect(r.end - r.start).toBeCloseTo(110, 0);
  });

  it('a history step is restored by position across a layout patch (rings to tracks and back)', async () => {
    const v = await make();
    v.view.setRange(Ring, 100, 160);
    frames!.flush();
    expect(v.view.jump(() => v.view.setRange(Ring, 300, 340))).toBe('moved');
    v.spec.applyPatch({ ops: [{ op: 'set', at: ['views', 'v', 'ring'], key: 'layout', value: 'tracks' }] });
    frames!.flush();
    expect(v.view.back()).toBe(true);
    frames!.flush();
    expect(v.view.get(Ring)).toMatchObject({ start: 100, end: 160 });
    expect(v.view.forward()).toBe(true);
    v.spec.applyPatch({ ops: [{ op: 'set', at: ['views', 'v', 'ring'], key: 'layout', value: 'rings' }] });
    frames!.flush();
    expect(v.view.back()).toBe(true);
    frames!.flush();
    expect(v.view.get(Ring)).toMatchObject({ start: 100, end: 160 });
    expect(v.view.get(Detail)).toMatchObject({ start: 100, end: 160 });
  });
});
