/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A ring nobody follows has no window (_spec/README.md "ViewModel", Layouts): Cmd/Ctrl+wheel and a plain wheel over it
// leave its section's viewport at the whole axis, whether the wheel is taken staying the binding table's; a ring a
// section follows still moves.

import type { Node } from '@molstar/core/repr/tree';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import type { RingPlacement } from '@molstar/seqx-core/render';
import { afterEach, describe, expect, it, vi } from 'vitest';
import { viewerItem } from '../config.js';
import { DefaultSeqViewerSpec } from '../default-spec.js';
import { createSeqViewer, type SeqViewer } from '../viewer.js';
import { type Frames, installFrames } from './frames.js';

const Ring = '/views/v/ring';
const Length = 600;

function spec(followed: boolean): Node {
  const b = createSeqSpec({ title: 'ring' });
  const c = b.data.sequence('c', { alphabet: 'dna', letters: 'ACGT'.repeat(Length / 4), cyclic: true });
  b.views.view('v', (v) => {
    const ring = v.section('ring', { of: c, layout: 'rings', labels: 'none' });
    ring.track('l').layer('l', { type: 'letters' });
    if (followed) v.section('detail', { of: c, follow: ring }).track('l').layer('l', { type: 'letters' });
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

async function make(followed: boolean) {
  vi.useFakeTimers({ toFake: ['Date', 'performance'] });
  frames = installFrames();
  const v = createSeqViewer({ ...DefaultSeqViewerSpec(), ui: undefined, config: [[viewerItem('view.jumpMs'), 0]] });
  viewer = v;
  await v.load(spec(followed));
  v.resize(640);
  await v.settle();
  frames.flush();
  v.input.mac = false;
  return v;
}

/** A page point on the ring's first band, a unit's angle. */
function pointAt(v: SeqViewer, unit: number) {
  const p = v.frame!.layout.sections.find((s) => s.section.address === Ring)!.placement as RingPlacement;
  const s = p.strips;
  p.warp.map(p.frameX(unit), s.bandY[0] + s.bandH[0] / 2);
  return { x: p.warp.px, y: p.warp.py };
}

describe('a ring nobody follows', () => {
  it('Cmd/Ctrl+wheel in and out and a plain wheel leave its viewport at the whole axis', async () => {
    const v = await make(false);
    const whole = { ...v.view.get(Ring)! };
    expect(whole.start).toBe(0);
    expect(whole.end).toBe(Length);
    const q = pointAt(v, 130);
    for (const e of [
      { deltaX: 0, deltaY: -100, ctrlKey: true, deltaMode: 0 },
      { deltaX: 0, deltaY: 100, ctrlKey: true, deltaMode: 0 },
      { deltaX: 0, deltaY: 100, deltaMode: 0 },
    ]) {
      v.input.wheel(q.x, q.y, e);
      frames!.flush();
      expect(v.view.get(Ring)).toEqual(whole);
    }
    // the model's own zoom and pan are idle too
    expect(v.view.zoom(Ring, 4, 130)).toBe(false);
    expect(v.view.pan(Ring, 20)).toBe(false);
    expect(v.view.panPx(Ring, 40)).toBe(false);
    expect(v.view.get(Ring)).toEqual(whole);
  });

  it('a followed ring still moves', async () => {
    const v = await make(true);
    const q = pointAt(v, 130);
    expect(v.input.wheel(q.x, q.y, { deltaX: 0, deltaY: -100, ctrlKey: true, deltaMode: 0 })).toBe(true);
    frames!.flush();
    const r = v.view.get(Ring)!;
    expect(r.end - r.start).toBeLessThan(Length);
    expect(v.view.pan(Ring, 20)).toBe(true);
  });
});
