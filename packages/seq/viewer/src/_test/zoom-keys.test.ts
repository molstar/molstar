/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { createSeqSpec } from '@molstar/seqx-core/lang';
import { afterEach, describe, expect, it } from 'vitest';
import { DefaultSeqViewerSpec } from '../default-spec.js';
import { createSeqViewer, type SeqViewer } from '../viewer.js';
import { type Frames, installFrames } from './frames.js';

// The zoom keys on the canvas (_spec/README.md "Bindings"; `view.keys`: `view.zoom-in` is = and +, `view.zoom-out` is
// -, `view.fit` is 0): the novice passes found + and = doing nothing, a tool's keys or the viewer's; they zoom a step.

let frames: Frames | undefined;
let viewer: SeqViewer | undefined;

afterEach(() => {
  viewer?.retire();
  viewer = undefined;
  frames?.restore();
  frames = undefined;
});

async function make() {
  frames = installFrames();
  const b = createSeqSpec();
  const p = b.data.sequence('p', { alphabet: 'protein', letters: 'MEEPQSDPSV'.repeat(30) });
  b.views.view('v', (v) => {
    const s = v.section('protein', { of: p });
    s.track('letters').layer('letters', { type: 'letters' });
  });
  viewer = createSeqViewer({ ...DefaultSeqViewerSpec(), ui: undefined });
  await viewer.load(b.build());
  viewer.resize(600, 200);
  await viewer.settle();
  frames.flush();
  viewer.draw();
  const section = viewer.view.resolved.value!.sections[0].address;
  const width = () => {
    const w = viewer!.view.get(section)!;
    return w.end - w.start;
  };
  return { v: viewer, width };
}

describe('the zoom keys', () => {
  it('= and + zoom in a step, - zooms out, 0 fits', async () => {
    const { v, width } = await make();
    const whole = width();
    for (const key of ['=', '+']) {
      v.view.fit(v.view.resolved.value!.sections[0].address);
      const before = width();
      expect(v.input.key({ key })).toBe(true);
      await v.settle();
      expect(width()).toBeLessThan(before);
    }
    const zoomed = width();
    expect(v.input.key({ key: '-' })).toBe(true);
    await v.settle();
    expect(width()).toBeGreaterThan(zoomed);
    expect(v.input.key({ key: '0' })).toBe(true);
    await v.settle();
    expect(width()).toBe(whole);
  });
});
