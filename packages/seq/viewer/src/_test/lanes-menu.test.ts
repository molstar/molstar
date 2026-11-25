/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A track's menu offers every lane (_spec/README.md "Behaviors", Headers): on a track whose lanes scroll inside (its
// content past its cap) "Show every lane" sets its `height` to `content`; on a track at `height: 'content'` "Cap the
// lanes" removes it; each one patch, one undo step; neither where it does not apply.

import { type Node, resolveAddress } from '@molstar/core/repr/tree';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { afterEach, describe, expect, it } from 'vitest';
import { DefaultSeqViewerSpec } from '../default-spec.js';
import { headerMenu } from '../header-entries.js';
import { type UiHistory, uiHistory } from '../ui-history.js';
import { createSeqViewer, type SeqViewer } from '../viewer.js';
import { type Frames, installFrames } from './frames.js';

const Section = '/views/v/s';
const Tall = `${Section}/tall`;
const Short = `${Section}/short`;

const stacked = (n: number) => ({
  start: Array.from({ length: n }, (_, i) => i),
  end: Array.from({ length: n }, () => 150),
  id: Array.from({ length: n }, (_, i) => `r${i}`),
});

function ringSpec(): Node {
  const b = createSeqSpec({ title: 'ring lanes' });
  const c = b.data.sequence('c', { alphabet: 'dna', letters: 'ACGT'.repeat(150), cyclic: true });
  const r = b.data.ranges('r', { target: c, rows: stacked(80) });
  b.views.view('v', (v) => {
    v.section('ring', { of: c, layout: 'rings', labels: 'none' }).track('t').layer('l', { type: 'ranges', data: r });
  });
  return b.build();
}

function spec(tallHeight?: 'content'): Node {
  const b = createSeqSpec({ title: 'every lane' });
  const s = b.data.sequence('s', { alphabet: 'dna', letters: 'ACGT'.repeat(100) });
  const tall = b.data.ranges('tall', { target: s, rows: stacked(80) });
  const short = b.data.ranges('short', { target: s, rows: stacked(2) });
  b.views.view('v', (v) => {
    const sec = v.section('s', { of: s });
    sec.track('tall', tallHeight ? ({ height: tallHeight } as never) : {}).layer('l', { type: 'ranges', data: tall });
    sec.track('short').layer('l', { type: 'ranges', data: short });
  });
  return b.build();
}

let frames: Frames | undefined;
let viewer: SeqViewer | undefined;
afterEach(() => {
  viewer?.retire();
  frames?.restore();
  viewer = undefined;
  frames = undefined;
});

async function make(tallHeight?: 'content', root: Node = spec(tallHeight)) {
  frames = installFrames();
  const v = createSeqViewer({ ...DefaultSeqViewerSpec(), ui: undefined, history: uiHistory() });
  viewer = v;
  await v.load(root);
  v.resize(600);
  await settled(v);
  return v;
}

async function settled(v: SeqViewer) {
  await v.settle();
  frames!.flush();
  v.draw();
  frames!.flush();
}

const items = (v: SeqViewer, track: string, section = Section) =>
  headerMenu(v, { section, track })!
    .groups.flat()
    .filter((i) => i.name === 'seq:every-lane' || i.name === 'seq:cap-lanes');
const heightOf = (v: SeqViewer, id: string) =>
  resolveAddress(v.spec.root.value, ['views', 'v', 's', id])!.params?.height;
const history = (v: SeqViewer) => (v.spec.history as UiHistory).stack;

describe('every lane from the track’s menu', () => {
  it('a track whose lanes scroll offers "Show every lane": one patch sets `height: content`, one undo step', async () => {
    const v = await make();
    const [first, ...rest] = items(v, Tall);
    expect(rest).toEqual([]);
    expect(first.label).toBe('Show every lane');
    expect(first.enabled).toBe(true);
    const before = history(v).steps.value.length;
    await v.execute(first.command!, first.args);
    await settled(v);
    expect(heightOf(v, 'tall')).toBe('content');
    expect(history(v).steps.value.length).toBe(before + 1);
    // at `content` now: the other entry, and "Show every lane" gone
    expect(items(v, Tall).map((i) => i.label)).toEqual(['Cap the lanes']);
  });

  it('"Cap the lanes" removes `height`, one undo step; undo brings back what was', async () => {
    const v = await make('content');
    const [cap, ...rest] = items(v, Tall);
    expect(rest).toEqual([]);
    expect(cap.label).toBe('Cap the lanes');
    const before = history(v).steps.value.length;
    await v.execute(cap.command!, cap.args);
    await settled(v);
    expect(heightOf(v, 'tall')).toBeUndefined();
    expect(history(v).steps.value.length).toBe(before + 1);
    await history(v).undo();
    await settled(v);
    expect(heightOf(v, 'tall')).toBe('content');
  });

  it('a track whose lanes all fit offers neither', async () => {
    const v = await make();
    expect(items(v, Short)).toEqual([]);
  });

  it('a ring’s track offers neither', async () => {
    const v = await make(undefined, ringSpec());
    expect(headerMenu(v, { section: '/views/v/ring', track: '/views/v/ring/t' })).toBeDefined();
    expect(items(v, '/views/v/ring/t', '/views/v/ring')).toEqual([]);
  });
});
