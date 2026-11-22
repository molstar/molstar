/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Node } from '@molstar/core/repr/tree';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { rangesLocation } from '@molstar/seqx-core/query';
import { afterEach, describe, expect, it } from 'vitest';
import { parseChord } from '../bindings.js';
import { AlignKey, BridgeChangeKinds, bridgeAria, carriedSections } from '../bridge-view.js';
import { changeKindOf } from '../commands.js';
import { DefaultViewerConfig } from '../config.js';
import { DefaultSeqViewerSpec } from '../default-spec.js';
import { describerFor } from '../describe.js';
import { createSeqViewer, type SeqViewer } from '../viewer.js';
import { type Frames, installFrames } from './frames.js';

// The viewer across a bridge (_spec/README.md "Bridges"; across-sections, Moving two sections joined by a bridge): a
// drag and a wheel over the band move both sides, `lock: px` carries a side's gesture to the other (each section moved
// once), `bridge.align`, `bridge.selectAcross`, the lock toggle, the cursor across the band, a ribbon's card, the
// bands' ARIA words, and Shift+A against the bindings.

let frames: Frames | undefined;
const viewers: SeqViewer[] = [];

afterEach(() => {
  for (const v of viewers.splice(0)) v.retire();
  frames?.restore();
  frames = undefined;
});

const Top = '/views/v/peach';
const Bottom = '/views/v/grape';
const Bridge = '/views/v/br';
const X = '/data/x';
const Y = '/data/y';

const letters = (n: number) => 'ACGT'.repeat(Math.ceil(n / 4)).slice(0, n);

/**
 * Two sections (x, 400; y, 300) joined by a bridge of three ribbons: x [0, 100) to y [0, 100), x [100, 200) to y
 * [200, 300) reversed, x [300, 400) to y [100, 200). With `self`, both sections show x (one genome against itself).
 */
function spec(o: { lock?: 'px'; self?: boolean; middle?: boolean } = {}): Node {
  const b = createSeqSpec({ title: 'bridge' });
  const x = b.data.sequence('x', { alphabet: 'dna', letters: letters(400), name: 'chrX' });
  const y = b.data.sequence('y', { alphabet: 'dna', letters: letters(300), name: 'chrY' });
  const yb = o.self ? x : y;
  const p = b.data.pairs('p', {
    a: x,
    b: yb,
    rows: {
      a_start: [0, 100, 300],
      a_end: [100, 200, 400],
      b_start: [0, 200, 100],
      b_end: [100, 300, 200],
      b_strand: ['+', '-', '+'],
      id: ['p0', 'p1', 'p2'],
    },
  });
  b.views.view('v', (v) => {
    v.section('peach', { of: x }, (s) => void s.track('t').layer('l', { type: 'letters' }));
    v.bridge('br', { a: Top, b: Bottom, ...(o.lock ? { lock: o.lock } : {}) } as never, (br) =>
      br.layer('ribbons', { type: 'pairs', data: p } as never),
    );
    v.section('grape', { of: yb }, (s) => void s.track('t').layer('l', { type: 'letters' }));
  });
  return b.build();
}

async function make(node: Node): Promise<SeqViewer> {
  frames = installFrames();
  const viewer = createSeqViewer({ ...DefaultSeqViewerSpec(), ui: undefined });
  viewers.push(viewer);
  await viewer.load(node);
  viewer.resize(800);
  await viewer.settle();
  frames.flush();
  viewer.draw();
  return viewer;
}

const range = (v: SeqViewer, s: string) => {
  const w = v.view.get(s)!;
  return [w.start, w.end];
};
const band = (v: SeqViewer) => v.frame!.layout.bridges!.list[0].box;
/** Re-draws after a move (the frame the input reads next). */
const redraw = (v: SeqViewer) => {
  frames!.flush();
  v.draw();
};

describe('a bridge in the viewer: its pointer, its lock, its commands (EC-206, EC-209)', () => {
  it('a drag over the band pans both sides by the same px, whatever the lock; selects nothing', async () => {
    const v = await make(spec());
    v.view.zoom(Top, 4);
    v.view.zoom(Bottom, 4);
    const [t0] = range(v, Top);
    const [b0] = range(v, Bottom);
    const box = band(v);
    const y = box.y + box.height / 2;
    // on a ribbon too: the band's, never a drag-select
    expect(v.input.down(box.x + 200, y)).toBe(true);
    v.input.move(box.x + 160, y);
    v.input.up();
    const pxTop = (range(v, Top)[0] - t0) / ((range(v, Top)[1] - range(v, Top)[0]) / v.view.get(Top)!.width);
    const pxBottom =
      (range(v, Bottom)[0] - b0) / ((range(v, Bottom)[1] - range(v, Bottom)[0]) / v.view.get(Bottom)!.width);
    expect(pxTop).toBeCloseTo(40, 6);
    expect(pxBottom).toBeCloseTo(40, 6);
    expect(v.selection.value.value.location.parts.length).toBe(0);
  });

  it('a wheel over the band zooms both sides about the pointer’s x', async () => {
    const v = await make(spec());
    const box = band(v);
    const x = box.x + box.width / 4;
    const before = { top: range(v, Top), bottom: range(v, Bottom) };
    v.input.wheel(x, box.y + 10, { deltaX: 0, deltaY: -100, deltaMode: 0, metaKey: true, ctrlKey: false });
    redraw(v);
    for (const [s, total] of [
      [Top, 400],
      [Bottom, 300],
    ] as const) {
      const [a, b] = range(v, s);
      expect(b - a).toBeLessThan(total);
      // the offset under the pointer stays: a quarter across, on both
      const w = before[s === Top ? 'top' : 'bottom'];
      const at = w[0] + ((w[1] - w[0]) * (x - box.x)) / v.view.get(s)!.width;
      expect(a + ((b - a) * (x - box.x)) / v.view.get(s)!.width).toBeCloseTo(at, 6);
    }
  });

  it('without a lock a side pans alone; with `lock: px` the other pans by the same px', async () => {
    const free = await make(spec());
    free.view.zoom(Top, 4);
    free.view.zoom(Bottom, 4);
    const b0 = range(free, Bottom);
    free.view.panPx(Top, 50);
    expect(range(free, Bottom)).toEqual(b0);

    const locked = await make(spec({ lock: 'px' }));
    locked.view.zoom(Top, 4);
    const bl = range(locked, Bottom);
    expect(bl[1] - bl[0]).toBeCloseTo(75, 6);
    locked.view.panPx(Top, 50);
    const w = locked.view.get(Bottom)!;
    expect(((w.start - bl[0]) * w.width) / (w.end - w.start)).toBeCloseTo(50, 6);
  });

  it('one move a section a gesture: a lock between sections a group already moves moves nothing twice', async () => {
    const v = await make(spec({ lock: 'px', self: true }));
    // one genome against itself: equal axes, one scroll group; the lock is redundant
    expect(carriedSections(v.view, [Top])).toEqual([Top]);
    v.view.zoom(Top, 4);
    const before = range(v, Top);
    v.view.panPx(Top, 40);
    const w = v.view.get(Top)!;
    expect(((w.start - before[0]) * w.width) / (w.end - w.start)).toBeCloseTo(40, 6);
    expect(range(v, Bottom)).toEqual(range(v, Top));
  });

  it('`bridge.align` (Shift+A at a side) brings the other side to the image of the window, one jump', async () => {
    const v = await make(spec());
    v.view.setRange(Top, 110, 190);
    redraw(v);
    expect(changeKindOf('bridge.align', BridgeChangeKinds)).toBe('navigate');
    const jumps = v.view.jumps.value;
    expect(v.runCommand('bridge.align', { section: Top })).toBe('Aligned');
    const [s, e] = range(v, Bottom);
    // the reversed ribbon's far end: y [200, 300)
    expect(s).toBeCloseTo(200, 6);
    expect(e).toBeCloseTo(300, 6);
    expect(v.view.jumps.value).toBe(jumps + 1);
    expect(v.view.back()).toBe(true);
    expect(range(v, Bottom)).toEqual([0, 300]);
  });

  it('`bridge.align` with nothing in view says so and moves nothing', async () => {
    const v = await make(spec());
    v.view.setRange(Top, 210, 290);
    redraw(v);
    expect(v.runCommand('bridge.align', { section: Top })).toBe('Nothing aligned in view');
    expect(range(v, Bottom)).toEqual([0, 300]);
  });

  it('`bridge.selectAcross` adds the far ends of the ribbons with an end in the selection, one change', async () => {
    const v = await make(spec());
    v.selection.add(rangesLocation(X, [20, 40]));
    expect(changeKindOf('bridge.selectAcross', BridgeChangeKinds)).toBe('session');
    expect(v.runCommand('bridge.selectAcross', { section: Top })).toBe('Selected across');
    const parts = v.selection.value.value.location.parts.map((p) => [p.system, [...p.ranges]]);
    expect(parts).toContainEqual([Y, [0, 100]]);
    expect(parts).toContainEqual([X, [20, 40]]);
  });

  it('the lock toggle is a patch: one undo step', async () => {
    const v = await make(spec());
    expect(changeKindOf('bridge.lock', BridgeChangeKinds)).toBe('patch');
    const revision = v.spec.head.value.revision;
    expect(await v.execute('bridge.lock', { bridge: Bridge })).toBe('Locked');
    await v.settle();
    expect(v.view.resolved.value!.bridges![0].lock).toBe('px');
    expect(v.spec.head.value.revision).not.toBe(revision);
    expect(await v.execute('bridge.lock', { bridge: Bridge })).toBe('Unlocked');
    await v.settle();
    expect(v.view.resolved.value!.bridges![0].lock).toBe('none');
  });

  it('the cursor crosses the band through the topmost ribbon covering it (a reversed one backwards)', async () => {
    const v = await make(spec());
    v.selection.setCursor({ section: Top, system: X, position: 120, keyboard: true });
    v.runCommand('cursor.down', {});
    const c = v.selection.cursor.value!;
    expect(c.section).toBe(Bottom);
    expect(c.system).toBe(Y);
    // 20 into [100, 200) is 20 from the end of [200, 300)
    expect(c.position).toBe(279);
    v.runCommand('cursor.up', {});
    expect(v.selection.cursor.value!.position).toBe(120);
  });

  it('a ribbon’s card names both ends, its strand and its aligned length', async () => {
    const v = await make(spec());
    const box = band(v);
    const ribbons = v.frame!.layout.bridges!.list[0].ribbons();
    const r = ribbons.find((x) => x.key === 'p1')!;
    // near its top end, where no other ribbon is
    const hit = v.hit((r.ax0 + r.ax1) / 2, box.y + 2);
    expect(hit?.bridge).toBe(0);
    v.hover.set(hit, 'test');
    const card = describerFor({ ...v, layout: () => v.frame?.layout }).describe(v.hover.current);
    expect(card?.title).toBe('peach chrX:101–200 ↔ grape chrY:201–300');
    expect(card?.text).toContain('reversed');
    expect(card?.text).toContain('100 bp');
  });

  it('the band’s ARIA words count the ribbons in view', async () => {
    const v = await make(spec());
    expect(bridgeAria(v.view, v.frame!.layout)).toEqual([
      { label: 'peach to grape', description: '3 ribbons between peach and grape in view' },
    ]);
    // the layer hidden: an empty band
    await v.execute('layer.hide', { address: `${Bridge}/ribbons` });
    await v.settle();
    redraw(v);
    expect(bridgeAria(v.view, v.frame!.layout)[0].description).toBe('no ribbons in view');
  });

  it('Shift+A is no other key’s: not a built-in command’s, not a binding’s', async () => {
    const v = await make(spec());
    const chord = parseChord('Shift+A');
    expect(chord).toEqual(AlignKey);
    for (const [id, keys] of Object.entries(DefaultViewerConfig.view.keys))
      for (const k of keys) expect(parseChord(k), id).not.toEqual(chord);
    for (const c of v.commands)
      if (c.id !== 'bridge.align' && c.shortcut) expect({ ...c.shortcut }, c.id).not.toEqual({ ...AlignKey });
    // a key typed into the view runs it, at the pointer's side or the keyboard cursor's
    v.view.setRange(Top, 110, 190);
    redraw(v);
    v.selection.setCursor({ section: Top, system: X, position: 150, keyboard: true });
    expect(v.input.key({ key: 'A', shiftKey: true })).toBe(true);
    expect(range(v, Bottom)[0]).toBeCloseTo(200, 6);
  });
});
