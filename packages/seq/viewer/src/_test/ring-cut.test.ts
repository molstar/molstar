/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The zoom of a ring nobody follows is a cut (_spec/README.md "Input", The ring's zoom is a cut): Cmd/Ctrl+wheel in over
// a chromosome opens it (`section.include`, one patch, one undo step), out over an opened one shows every block again,
// one cut a gesture, the hole says once what the ring shows and changes nothing; a followed ring keeps its window's zoom.

import type { Node } from '@molstar/core/repr/tree';
import { assembly, createSeqSpec, genome } from '@molstar/seqx-core/lang';
import type { RingPlacement } from '@molstar/seqx-core/render';
import { afterEach, describe, expect, it, vi } from 'vitest';
import { DefaultSeqViewerSpec } from '../default-spec.js';
import { type UiHistory, uiHistory } from '../ui-history.js';
import { createSeqViewer, type SeqViewer } from '../viewer.js';
import { type Frames, installFrames } from './frames.js';

const Ring = '/views/main/ring';
const Names = Array.from({ length: 6 }, (_, i) => `chr${i + 1}`);
const Hole = 'The ring shows the whole genome: open a chromosome from its label’s menu or the square';

function spec(followed = false, include?: readonly string[]): Node {
  const b = createSeqSpec({ title: 'the ring cut' });
  const g = b.data.assembly(
    'g',
    assembly(
      'Six',
      Names.map((name) => ({ name, length: 1_000_000 })),
    ),
  );
  b.views.view('main', (v) => {
    genome(v, 'ring', g, { layout: 'rings', ...(include ? { include } : {}) } as never);
    if (followed) v.section('detail', { of: g, follow: Ring } as never);
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

async function make(followed = false, include?: readonly string[]) {
  vi.useFakeTimers({ toFake: ['Date', 'performance'] });
  frames = installFrames();
  const v = createSeqViewer({ ...DefaultSeqViewerSpec(), ui: undefined, history: uiHistory() });
  viewer = v;
  await v.load(spec(followed, include));
  v.resize(600, 600);
  await v.settle();
  frames.flush();
  v.draw();
  v.input.mac = false;
  return v;
}

const placement = (v: SeqViewer) => v.placementOf(Ring) as RingPlacement;
const blocks = (v: SeqViewer) => v.view.resolved.value!.sections.find((s) => s.address === Ring)!.blocks!;

/** A page point on the ring's first band at the middle of block `k` (of what is shown). */
function over(v: SeqViewer, k: number) {
  const p = placement(v);
  const b = blocks(v)[k];
  const s = p.strips;
  p.warp.map(p.frameX((b.offset ?? 0) + 1000), s.bandY[0] + s.bandH[0] / 2);
  return { x: p.warp.px, y: p.warp.py };
}

/** The ring's hole: its box's centre. */
function hole(v: SeqViewer) {
  const { box } = placement(v);
  return { x: box.x + box.width / 2, y: box.y + box.height / 2 };
}

const include = (v: SeqViewer) => {
  const views = v.spec.root.value.children!.find((c) => c.kind === 'views')!;
  return views.children![0].children!.find((c) => c.id === 'ring')!.params?.include as string[] | undefined;
};
const undo = (v: SeqViewer) => (v.spec.history as UiHistory).stack.undo();
const zoom = (v: SeqViewer, at: { x: number; y: number }, deltaY: number) =>
  v.input.wheel(at.x, at.y, { deltaX: 0, deltaY, deltaMode: 0, ctrlKey: true, metaKey: false });
/** Time passes: the wheel rests. */
const rest = async (v: SeqViewer) => {
  vi.advanceTimersByTime(1000);
  await v.settle();
  frames!.flush();
  v.draw();
};

describe('the ring’s zoom is a cut', () => {
  it('in over a chromosome opens it: one patch, one undo step, the event taken', async () => {
    const v = await make();
    expect(include(v)).toBeUndefined();
    expect(zoom(v, over(v, 2), -100)).toBe(true);
    await rest(v);
    expect(include(v)).toEqual(['chr3']);
    expect(blocks(v).length).toBe(1);
    expect((await undo(v)).ok).toBe(true);
    expect(include(v)).toBeUndefined();
  });

  it('out over an opened chromosome shows every block again: one patch, one undo step', async () => {
    const v = await make(false, ['chr2']);
    expect(blocks(v).length).toBe(1);
    expect(zoom(v, over(v, 0), 100)).toBe(true);
    await rest(v);
    expect(include(v)).toBeUndefined();
    expect(blocks(v).length).toBe(Names.length);
    await undo(v);
    expect(include(v)).toEqual(['chr2']);
  });

  it('out on the whole genome and in on what is shown alone change nothing, and are taken', async () => {
    const v = await make();
    expect(zoom(v, over(v, 1), 100)).toBe(true);
    await rest(v);
    expect(include(v)).toBeUndefined();
  });

  it('in on a chromosome shown alone changes nothing, and is taken', async () => {
    const v = await make(false, ['chr4']);
    expect(zoom(v, over(v, 0), -100)).toBe(true);
    await rest(v);
    expect(include(v)).toEqual(['chr4']);
  });

  it('one cut a gesture: more deltas wait for the wheel to rest', async () => {
    const v = await make();
    const at = over(v, 0);
    // under a notch: nothing yet
    expect(zoom(v, at, -20)).toBe(true);
    await v.settle();
    expect(include(v)).toBeUndefined();
    vi.advanceTimersByTime(50);
    expect(zoom(v, at, -40)).toBe(true);
    await v.settle();
    expect(include(v)).toEqual(['chr1']);
    // the same gesture goes on: out is no second cut
    vi.advanceTimersByTime(50);
    expect(zoom(v, over(v, 0), 200)).toBe(true);
    await v.settle();
    expect(include(v)).toEqual(['chr1']);
    // the wheel rests: out shows all, a new gesture
    await rest(v);
    expect(zoom(v, over(v, 0), 100)).toBe(true);
    await rest(v);
    expect(include(v)).toBeUndefined();
  });

  it('over the hole nothing changes and the first such wheel of a gesture says so once', async () => {
    const v = await make();
    const said: string[] = [];
    v.notice.listen(({ text }) => said.push(text));
    const at = hole(v);
    expect(zoom(v, at, -100)).toBe(true);
    expect(zoom(v, at, -100)).toBe(true);
    await v.settle();
    expect(include(v)).toBeUndefined();
    expect(said).toEqual([Hole]);
    // a new gesture says it again
    await rest(v);
    expect(zoom(v, at, -100)).toBe(true);
    expect(said).toEqual([Hole, Hole]);
  });

  it('off the exact centre the hole is the hole too, and a pinch cuts as a wheel does', async () => {
    const v = await make();
    const c = hole(v);
    const said: string[] = [];
    v.notice.listen(({ text }) => said.push(text));
    expect(zoom(v, { x: c.x + 15, y: c.y + 9 }, -100)).toBe(true);
    expect(include(v)).toBeUndefined();
    expect(said).toEqual([Hole]);
    await rest(v);
    // a trackpad's pinch (Ctrl, small pixel deltas): in
    const at = over(v, 4);
    expect(v.input.wheel(at.x, at.y, { deltaX: 0, deltaY: -60, deltaMode: 0, ctrlKey: true })).toBe(true);
    await rest(v);
    expect(include(v)).toEqual(['chr5']);
  });

  it('a ring a section follows keeps the zoom of its window', async () => {
    const v = await make(true);
    const before = { ...v.view.get(Ring)! };
    expect(zoom(v, over(v, 2), -100)).toBe(true);
    frames!.flush();
    await v.settle();
    expect(include(v)).toBeUndefined();
    const after = v.view.get(Ring)!;
    expect(after.end - after.start).toBeLessThan(before.end - before.start);
  });
});
