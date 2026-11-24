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
const Whole = 'The ring shows the whole genome: open a chromosome from its label’s menu or the square';
const Alone = (name: string) => `The ring shows ${name}: wheel out over it to show every chromosome`;
const Filtered = 'The ring shows the chromosomes its view selects: open one from its label’s menu or the square';

function spec(followed = false, include?: unknown): Node {
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

async function make(followed = false, include?: unknown) {
  vi.useFakeTimers({ toFake: ['Date', 'performance'] });
  frames = installFrames();
  const v = createSeqViewer({
    ...DefaultSeqViewerSpec(),
    ui: undefined,
    history: uiHistory(),
  });
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

  /** What a zoom says, once the gesture is over. */
  const saying = (v: SeqViewer) => {
    const said: string[] = [];
    v.notice.listen(({ text }) => said.push(text));
    return said;
  };

  it('out on the whole genome changes nothing and says what the ring shows, once a gesture', async () => {
    const v = await make();
    const said = saying(v);
    expect(zoom(v, over(v, 1), 100)).toBe(true);
    expect(zoom(v, over(v, 1), 100)).toBe(true);
    await rest(v);
    expect(include(v)).toBeUndefined();
    expect(said).toEqual([Whole]);
  });

  it('in on a chromosome shown alone changes nothing and says which, wheel out to show every one', async () => {
    const v = await make(false, ['chr4']);
    const said = saying(v);
    expect(zoom(v, over(v, 0), -100)).toBe(true);
    await rest(v);
    expect(include(v)).toEqual(['chr4']);
    expect(said).toEqual([Alone('chr4')]);
  });

  it('out over the hole of an opened ring changes nothing and says which chromosome it shows', async () => {
    const v = await make(false, ['chr4']);
    const said = saying(v);
    const c = hole(v);
    expect(zoom(v, { x: c.x + 12, y: c.y + 7 }, 100)).toBe(true);
    await rest(v);
    expect(include(v)).toEqual(['chr4']);
    expect(said).toEqual([Alone('chr4')]);
  });

  it('outside the ring’s box the zoom is not taken and nothing is said', async () => {
    const v = await make();
    const said = saying(v);
    const { box } = placement(v);
    // the margin beside the ring, inside the section's height
    expect(box.x).toBeGreaterThan(2);
    // the bindings' own zoom of the window (taken, as any zoom is), never the cut
    zoom(v, { x: box.x - 2, y: box.y + box.height / 2 }, -100);
    zoom(v, { x: box.x - 2, y: box.y + box.height / 2 }, 100);
    await rest(v);
    expect(include(v)).toBeUndefined();
    expect(said).toEqual([]);
  });

  it('an author’s filter is never cleared: a zoom over a ring whose include is a filter says so and keeps it', async () => {
    const v = await make(false, { min_length: 1 });
    const said = saying(v);
    expect(zoom(v, over(v, 0), 100)).toBe(true);
    await rest(v);
    expect(zoom(v, over(v, 1), -100)).toBe(true);
    await rest(v);
    expect(include(v)).toEqual({ min_length: 1 });
    expect(said).toEqual([Filtered, Filtered]);
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
    expect(said).toEqual([Whole]);
    // a new gesture says it again
    await rest(v);
    expect(zoom(v, at, -100)).toBe(true);
    expect(said).toEqual([Whole, Whole]);
  });

  it('off the exact centre the hole is the hole too, and a pinch cuts as a wheel does', async () => {
    const v = await make();
    const c = hole(v);
    const said: string[] = [];
    v.notice.listen(({ text }) => said.push(text));
    expect(zoom(v, { x: c.x + 15, y: c.y + 9 }, -100)).toBe(true);
    expect(include(v)).toBeUndefined();
    expect(said).toEqual([Whole]);
    await rest(v);
    // a trackpad's pinch (Ctrl, small pixel deltas): in
    const at = over(v, 4);
    expect(v.input.wheel(at.x, at.y, { deltaX: 0, deltaY: -60, deltaMode: 0, ctrlKey: true })).toBe(true);
    await rest(v);
    expect(include(v)).toEqual(['chr5']);
  });

  it('a detached follower ring keeps the zoom of its window (the resolved section follows, the view model does not)', async () => {
    // a section of blocks resolves no follow in this spec (a genome's blocks have no one system a path joins), so the
    // detached state is made as the view model has it: `follow` on the resolved section, the view's own group detached
    const v = await make();
    const s = v.view.resolved.value!.sections.find((x) => x.address === Ring)! as { follow?: unknown };
    s.follow = { leader: '/views/main/other' };
    expect(v.view.isFollowing(Ring)).toBe(false);
    const said: string[] = [];
    v.notice.listen(({ text }) => said.push(text));
    zoom(v, over(v, 1), -100);
    await rest(v);
    expect(include(v)).toBeUndefined();
    expect(said).toEqual([]);
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

describe('a double-click on a chromosome opens it', () => {
  const dbl = (v: SeqViewer, at: { x: number; y: number }, shiftKey = false) =>
    v.input.doubleClick(at.x, at.y, { shiftKey });

  it('over a chromosome’s arc: the cut, one patch, one undo step, the event taken', async () => {
    const v = await make();
    const steps = (v.spec.history as UiHistory).stack.steps.value.length;
    expect(dbl(v, over(v, 3))).toBe(true);
    await rest(v);
    expect(include(v)).toEqual(['chr4']);
    expect(blocks(v).length).toBe(1);
    expect((v.spec.history as UiHistory).stack.steps.value.length).toBe(steps + 1);
    expect((await undo(v)).ok).toBe(true);
    expect(include(v)).toBeUndefined();
  });

  it('over the hole it is not the cut’s; on a ring that shows one chromosome it is not either', async () => {
    const v = await make();
    dbl(v, hole(v));
    await rest(v);
    expect(include(v)).toBeUndefined();
    const one = await make(false, ['chr2']);
    dbl(one, over(one, 0));
    await rest(one);
    expect(include(one)).toEqual(['chr2']);
  });

  it('a Shift+double-click is the bindings’ own (zoom out), not the cut', async () => {
    const v = await make();
    dbl(v, over(v, 1), true);
    await rest(v);
    expect(include(v)).toBeUndefined();
  });

  it('an author’s filter is kept and said', async () => {
    const v = await make(false, { min_length: 1 });
    const said: string[] = [];
    v.notice.listen(({ text }) => said.push(text));
    expect(dbl(v, over(v, 1))).toBe(true);
    await rest(v);
    expect(include(v)).toEqual({ min_length: 1 });
    expect(said).toEqual([Filtered]);
  });

  it('a followed ring keeps its window’s double-click', async () => {
    const v = await make(true);
    dbl(v, over(v, 2));
    await rest(v);
    expect(include(v)).toBeUndefined();
  });
});
