/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A track's lanes scrolling inside it in the viewer (core's _spec/README.md "Render", Track heights): the wheel's
// vertical delta over a band whose lanes overflow it scrolls them (session state, no undo step), the header and the
// tracks around stay where they are, a band at its end gives the wheel back, and a band that shows all its lanes takes
// no vertical wheel.

import type { Node } from '@molstar/core/repr/tree';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import type { RangesLay } from '@molstar/seqx-core/render';
import { afterEach, describe, expect, it, vi } from 'vitest';
import { viewerItem } from '../config.js';
import { DefaultSeqViewerSpec } from '../default-spec.js';
import { createSeqViewer, type SeqViewer } from '../viewer.js';
import { type Frames, installFrames } from './frames.js';

const Tall = '/views/v/s/tall';
const Short = '/views/v/s/short';

function stacked(n: number) {
  return {
    start: Array.from({ length: n }, (_, i) => i),
    end: Array.from({ length: n }, () => 150),
    id: Array.from({ length: n }, (_, i) => `r${i}`),
  };
}

function spec(): Node {
  const b = createSeqSpec({ title: 'lanes' });
  const s = b.data.sequence('s', { alphabet: 'dna', letters: 'ACGT'.repeat(100) });
  const tall = b.data.ranges('tall', { target: s, rows: stacked(40) });
  const short = b.data.ranges('short', { target: s, rows: stacked(2) });
  b.views.view('v', (v) => {
    const sec = v.section('s', { of: s, heights: 'fixed' } as never);
    sec.track('tall', { height: 240 } as never).layer('l', { type: 'ranges', data: tall });
    sec.track('short').layer('l', { type: 'ranges', data: short });
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
  v.resize(600);
  await v.settle();
  frames.flush();
  v.input.mac = false;
  return v;
}

const lay = (v: SeqViewer, track: string) => {
  for (const s of v.frame!.layout.sections)
    for (const t of s.tracks) for (const l of t.layers) if (l.key.startsWith(track)) return l.lay as RangesLay;
  throw new Error(`no lanes in ${track}`);
};

describe('a track’s lanes scrolling in the viewer', () => {
  it('the wheel over an overflowing band scrolls its lanes (session state), the header and the tracks around stay', async () => {
    const v = await make();
    const l = lay(v, Tall);
    expect(l.scrolls).toBe(true);
    expect(l.scrollMax).toBeGreaterThan(0);
    const tracks = () => v.frame!.layout.sections[0].tracks.map((t) => [t.box.y, t.box.height]);
    const before = tracks();
    const revision = v.spec.revision.value;
    const x = l.box.x + 20;
    const y = l.box.y + 30;
    v.input.move(x, y);
    // 42 px of wheel at 14 px a lane: three lanes
    expect(v.input.wheel(x, y, { deltaX: 0, deltaY: 42, deltaMode: 0 })).toBe(true);
    frames!.flush();
    v.draw();
    expect(v.view.band(Tall)?.bandTop).toBe(3);
    expect(lay(v, Tall).top).toBe(3);
    // no undo step; nothing around it moved
    expect(v.spec.revision.value).toBe(revision);
    expect(tracks()).toEqual(before);
  });

  it('at the bottom the wheel is given back (scroll chaining): the band keeps its end', async () => {
    const v = await make();
    const l = lay(v, Tall);
    const x = l.box.x + 20;
    const y = l.box.y + 30;
    v.input.move(x, y);
    v.view.setBand(Tall, { bandTop: l.scrollMax });
    frames!.flush();
    v.draw();
    expect(lay(v, Tall).top).toBe(l.scrollMax);
    v.input.wheel(x, y, { deltaX: 0, deltaY: 42, deltaMode: 0 });
    frames!.flush();
    expect(v.view.band(Tall)?.bandTop).toBe(l.scrollMax);
    // and upward it scrolls again
    expect(v.input.wheel(x, y, { deltaX: 0, deltaY: -28, deltaMode: 0 })).toBe(true);
    frames!.flush();
    expect(v.view.band(Tall)?.bandTop).toBe(l.scrollMax - 2);
  });

  it('a band that shows all its lanes takes no vertical wheel', async () => {
    const v = await make();
    const l = lay(v, Short);
    expect(l.scrollMax).toBe(0);
    const x = l.box.x + 20;
    const y = l.box.y + 10;
    v.input.move(x, y);
    v.input.wheel(x, y, { deltaX: 0, deltaY: 42, deltaMode: 0 });
    frames!.flush();
    expect(v.view.band(Short)).toBeUndefined();
  });
});
