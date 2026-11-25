/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// `height: 'content'` in a viewer (core's _spec/README.md "Render", Track heights and Section scroll): the canvas still
// follows `canvas.maxHeight` (a section past it scrolls the canvas), and a patch switches a track's `height` to
// `content` and back, the band laid again each time.

import type { Node } from '@molstar/core/repr/tree';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { afterEach, describe, expect, it, vi } from 'vitest';
import { viewerItem } from '../config.js';
import { DefaultSeqViewerSpec } from '../default-spec.js';
import { createSeqViewer, type SeqViewer } from '../viewer.js';
import { type Frames, installFrames } from './frames.js';

function spec(track: number | 'content' = 240): Node {
  const b = createSeqSpec({ title: 'content height' });
  const s = b.data.sequence('s', { alphabet: 'dna', letters: 'ACGT'.repeat(100) });
  const tall = b.data.ranges('tall', {
    target: s,
    rows: {
      start: Array.from({ length: 40 }, (_, i) => i),
      end: Array.from({ length: 40 }, () => 150),
      id: Array.from({ length: 40 }, (_, i) => `r${i}`),
    },
  });
  b.views.view('v', (v) => {
    const sec = v.section('s', { of: s, height: 'content' } as never);
    sec.track('tall', { height: track } as never).layer('l', { type: 'ranges', data: tall });
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

async function make(config: [string, unknown][] = [], track: number | 'content' = 240) {
  vi.useFakeTimers({ toFake: ['Date', 'performance'] });
  frames = installFrames();
  const v = createSeqViewer({
    ...DefaultSeqViewerSpec(),
    ui: undefined,
    config: config.map(([k, x]) => [viewerItem(k as never), x]) as never,
  });
  viewer = v;
  await v.load(spec(track));
  v.resize(600);
  await v.settle();
  frames.flush();
  v.draw();
  return v;
}

const trackBox = (v: SeqViewer) => v.frame!.layout.sections[0].tracks[0].box.height;

const redraw = async (v: SeqViewer) => {
  await v.settle();
  frames!.flush();
  v.draw();
  frames!.flush();
  v.draw();
};

describe('height: content in a viewer', () => {
  it('canvas.maxHeight still bounds the canvas: a section past it is its content, the canvas scrolls', async () => {
    const v = await make([['canvas.maxHeight', 300]], 'content');
    const s = v.frame!.layout.sections[0];
    expect(s.placement.contentHeight).toBeGreaterThan(300);
    expect(s.placement.box.height).toBe(s.placement.contentHeight);
    expect(v.height).toBe(300);
  });

  it('a patch switches a track’s height to content and back', async () => {
    const v = await make();
    expect(trackBox(v)).toBe(240);
    const at = ['views', 'v', 's', 'tall'];
    const on = await v.spec.patch('Content', { ops: [{ op: 'set', at, key: 'height', value: 'content' }] });
    expect(on.ok && on.changed).toBe(true);
    await redraw(v);
    expect(trackBox(v)).toBeGreaterThan(240);
    const off = await v.spec.patch('Px', { ops: [{ op: 'set', at, key: 'height', value: 240 }] });
    expect(off.ok && off.changed).toBe(true);
    await redraw(v);
    expect(trackBox(v)).toBe(240);
  });
});
