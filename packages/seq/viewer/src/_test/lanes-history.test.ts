/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A lanes band in the view's history (_spec/README.md "ViewModel", Bands and History): a step keeps the band by the item
// at its top lane, so a back after the lanes were laid again (rows removed above it) returns to the same item, not to
// the same lane number.

import type { Node } from '@molstar/core/repr/tree';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import type { RangesLay } from '@molstar/seqx-core/render';
import { afterEach, describe, expect, it, vi } from 'vitest';
import { viewerItem } from '../config.js';
import { DefaultSeqViewerSpec } from '../default-spec.js';
import { createSeqViewer, type SeqViewer } from '../viewer.js';
import { type Frames, installFrames } from './frames.js';

const Tall = '/views/v/s/tall';
const Section = '/views/v/s';

/** `n` items overlapping everywhere, one lane each, named from `from`. */
function stacked(from: number, n: number) {
  return {
    start: Array.from({ length: n }, (_, i) => i),
    end: Array.from({ length: n }, () => 150),
    id: Array.from({ length: n }, (_, i) => `r${from + i}`),
  };
}

function spec(): Node {
  const b = createSeqSpec({ title: 'lanes history' });
  const s = b.data.sequence('s', { alphabet: 'dna', letters: 'ACGT'.repeat(100) });
  const tall = b.data.ranges('tall', { target: s, rows: stacked(0, 40) });
  b.views.view('v', (v) => {
    const sec = v.section('s', { of: s, heights: 'fixed' } as never);
    sec.track('tall', { height: 240 } as never).layer('l', { type: 'ranges', data: tall });
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
  return v;
}

const lay = (v: SeqViewer) => {
  for (const s of v.frame!.layout.sections)
    for (const t of s.tracks) for (const l of t.layers) if (l.key.startsWith(Tall)) return l.lay as RangesLay;
  throw new Error('no lanes');
};

const redraw = (v: SeqViewer) => {
  vi.advanceTimersByTime(v.config.view.jumpMs + 1);
  frames!.flush();
  v.draw();
  frames!.flush();
};

describe('a lanes band in the view’s history', () => {
  it('back returns the band to the same item when the lanes were laid again between the steps', async () => {
    const v = await make();
    v.view.setBand(Tall, { bandTop: 10 });
    redraw(v);
    const top = lay(v);
    expect(top.top).toBe(10);
    const item = top.ids[top.lanes.findIndex((l, i) => i < top.items && l === 10)];
    expect(item).toMatch(/r10/);

    // a jump pushes the step holding the band at r10
    v.runCommand('view.focus', { section: Section, focus: { start: 20, end: 60 } });
    redraw(v);

    // the five first rows leave the data: r10 is in lane 5 once the lanes are packed again
    const o = await v.spec.patch('Fewer rows', {
      ops: [{ op: 'set', at: ['data', 'tall'], key: 'rows', value: stacked(5, 35) }],
    });
    expect(o.ok && o.changed).toBe(true);
    await v.settle();
    redraw(v);
    expect(lay(v).sticky.laneOf(item)).toBe(5);

    expect(v.view.back()).toBe(true);
    redraw(v);
    // the lane number the step was taken at (10) shows r15; the item's own shows r10
    expect(v.view.band(Tall)?.bandTop).toBe(5);
    expect(lay(v).top).toBe(5);
  });

  it('a band kept at an item the data no longer holds keeps its lane number (clamped)', async () => {
    const v = await make();
    v.view.setBand(Tall, { bandTop: 10 });
    redraw(v);
    v.runCommand('view.focus', { section: Section, focus: { start: 20, end: 60 } });
    redraw(v);
    // r10 leaves the data
    const rows = stacked(0, 40);
    const keep = rows.id.map((id, i) => (id === 'r10' ? -1 : i)).filter((i) => i >= 0);
    const o = await v.spec.patch('Without r10', {
      ops: [
        {
          op: 'set',
          at: ['data', 'tall'],
          key: 'rows',
          value: {
            start: keep.map((i) => rows.start[i]),
            end: keep.map((i) => rows.end[i]),
            id: keep.map((i) => rows.id[i]),
          },
        },
      ],
    });
    expect(o.ok && o.changed).toBe(true);
    await v.settle();
    redraw(v);
    expect(v.view.back()).toBe(true);
    redraw(v);
    expect(v.view.band(Tall)?.bandTop).toBe(10);
  });
});
