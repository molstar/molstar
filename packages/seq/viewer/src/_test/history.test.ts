/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Zooming out (_spec/README.md "ViewModel", History; "Input"; "Behaviors"), headless on a fake clock: the view's
// history of jumps (a focus, a fit, a bookmark's focus) with back and forward, a move by hand pushing nothing; a
// double-click on what is focused going back; Shift+double-click zooming out about the pointer; the zoom keys about the
// keyboard cursor or the centre; the side buttons; a jump's motion, cut under reduced motion.

import type { Node } from '@molstar/core/repr/tree';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import type { RangesLay } from '@molstar/seqx-core/render';
import { afterEach, describe, expect, it, vi } from 'vitest';
import { parseChord } from '../bindings.js';
import { viewerItem } from '../config.js';
import { DefaultSeqViewerSpec } from '../default-spec.js';
import { createSeqViewer, type SeqViewer } from '../viewer.js';
import { type Frames, installFrames } from './frames.js';

const Protein = '/views/main/protein';
const Other = '/views/main/other';

function spec(): Node {
  const b = createSeqSpec({ title: 'history' });
  const p = b.data.sequence('p', { alphabet: 'protein', letters: 'MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDD' });
  const q = b.data.sequence('q', { alphabet: 'protein', letters: 'MKV'.repeat(20) });
  const domains = b.data.ranges('domains', {
    target: p,
    rows: { start: new Int32Array([0, 20]), end: new Int32Array([10, 30]), id: ['A', 'B'] },
  });
  b.refs.ref('b', { to: `${domains}/(item :id B)`, label: 'Domain B' });
  b.views.view('main', (v) => {
    const s = v.section('protein', { of: p });
    s.track('letters').layer('l', { type: 'letters' });
    s.track('domains').layer('d', { type: 'ranges', data: domains });
    // a second leading section, on another system: a back restores both together
    v.section('other', { of: q }).track('letters').layer('l', { type: 'letters' });
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

async function make(config: [unknown, unknown][] = []) {
  vi.useFakeTimers({ toFake: ['Date', 'performance'] });
  frames = installFrames();
  const v = createSeqViewer({
    ...DefaultSeqViewerSpec(),
    ui: undefined,
    config: [[viewerItem('viewport.maxPxPerPosition'), 1000], ...(config as never[])],
  });
  viewer = v;
  await v.load(spec());
  v.resize(600);
  await v.settle();
  frames.flush();
  return v;
}

const range = (v: SeqViewer, s = Protein) => {
  const r = v.view.get(s)!;
  return [r.start, r.end];
};
const run = (v: SeqViewer, id: string, ctx: unknown = { section: Protein }) => v.runCommand(id, ctx);

/** Lets a jump's motion end and draws. */
function settleMotion(v: SeqViewer) {
  vi.advanceTimersByTime(v.config.view.jumpMs + 1);
  frames!.flush();
  v.draw();
}

describe('the view’s history of jumps', () => {
  it('back and forward over a focus, a fit and a bookmark’s focus; a move by hand pushes nothing', async () => {
    const v = await make();
    const fitted = range(v);
    expect(run(v, 'view.back')).toBe('No view to go back to');
    run(v, 'view.focus', { section: Protein, focus: { start: 2, end: 8 } });
    expect(range(v)).toEqual([2, 8]);
    run(v, 'view.fit');
    expect(range(v)).toEqual(fitted);
    expect(run(v, 'address.focus', { address: '/refs/b' })).toBe('Focused');
    expect(range(v)).toEqual([20, 30]);
    expect(v.view.history.value).toEqual({ back: 3, forward: 0 });

    expect(run(v, 'view.back')).toBe('Back');
    expect(range(v)).toEqual(fitted);
    run(v, 'view.back');
    expect(range(v)).toEqual([2, 8]);
    run(v, 'view.back');
    expect(range(v)).toEqual(fitted);
    expect(v.view.history.value).toEqual({ back: 0, forward: 3 });
    expect(run(v, 'view.back')).toBe('No view to go back to');
    run(v, 'view.forward');
    run(v, 'view.forward');
    expect(range(v)).toEqual(fitted);
    run(v, 'view.forward');
    expect(range(v)).toEqual([20, 30]);
    expect(run(v, 'view.forward')).toBe('No view to go forward to');

    // a pan by hand: no step; the jump after it pushes where the hand left the view, and drops the steps forward
    run(v, 'view.back');
    v.view.pan(Protein, 3);
    expect(v.view.history.value).toEqual({ back: 2, forward: 1 });
    const panned = range(v);
    run(v, 'view.focus', { section: Protein, focus: { start: 4, end: 9 } });
    expect(v.view.history.value).toEqual({ back: 3, forward: 0 });
    run(v, 'view.back');
    expect(range(v)).toEqual(panned);
  });

  it('a back restores every leading section’s range together; the history is bounded', async () => {
    const v = await make([[viewerItem('view.historyLimit'), 3]]);
    const before = [range(v), range(v, Other)];
    v.view.jump(() => v.view.setRange(Protein, 5, 15) && v.view.setRange(Other, 10, 20));
    expect([range(v), range(v, Other)]).toEqual([
      [5, 15],
      [10, 20],
    ]);
    expect(v.view.back()).toBe(true);
    expect([range(v), range(v, Other)]).toEqual(before);
    for (let i = 0; i < 6; i++) run(v, 'view.focus', { section: Protein, focus: { start: i, end: i + 5 } });
    expect(v.view.history.value.back).toBe(3);
    while (v.view.back());
    expect(range(v)).toEqual([2, 7]);
  });

  it('a step zoom animates and is no step; the same jump twice is one', async () => {
    const v = await make();
    run(v, 'view.zoom-in');
    expect(v.view.history.value.back).toBe(0);
    run(v, 'view.focus', { section: Protein, focus: { start: 2, end: 8 } });
    run(v, 'view.focus', { section: Protein, focus: { start: 2, end: 8 } });
    // the second found the view there: back (`view.backOnFocused`)
    expect(v.view.history.value).toEqual({ back: 0, forward: 1 });
  });
});

describe('double-clicks and keys', () => {
  /** The centre of domain B's mark, in canvas px. */
  function domainB(v: SeqViewer) {
    v.draw();
    const lay = v.frame!.layout.sections[0].tracks[1].layers[0].lay as RangesLay;
    const j = Array.from({ length: lay.m }, (_, k) => k).find((k) => lay.keys[lay.mitem[k]] === 'B')!;
    return [(lay.mx0[j] + lay.mx1[j]) / 2, lay.box.y + (lay.mlane[j] + 0.5) * lay.laneHeight] as const;
  }

  it('a double-click focuses an element; on the element already focused, it goes back', async () => {
    const v = await make();
    const fitted = range(v);
    const [x, y] = domainB(v);
    expect(v.input.doubleClick(x, y)).toBe(true);
    expect(range(v)).toEqual([20, 30]);
    settleMotion(v);
    const [x2, y2] = domainB(v);
    v.input.doubleClick(x2, y2);
    expect(range(v)).toEqual(fitted);
    // configured off: it stays
    viewer!.retire();
    const w = await make([[viewerItem('view.backOnFocused'), false]]);
    w.input.doubleClick(...domainB(w));
    settleMotion(w);
    w.input.doubleClick(...domainB(w));
    expect(range(w)).toEqual([20, 30]);
  });

  it('Shift+double-click zooms out about the pointer by the double-click’s factor', async () => {
    const v = await make();
    v.view.setRange(Protein, 10, 20);
    v.draw();
    const vp = v.view.get(Protein)!;
    // the pointer over offset 12
    const x = v.dataX + ((12 - vp.start) / (vp.end - vp.start)) * vp.width;
    const y = v.frame!.layout.sections[0].tracks[0].box.y + 2;
    expect(v.input.doubleClick(x, y, { shiftKey: true })).toBe(true);
    const [s, e] = range(v);
    expect(e - s).toBeCloseTo(20);
    // offset 12 stays under the pointer
    expect((12 - s) / (e - s)).toBeCloseTo(0.2);
    expect(v.view.history.value.back).toBe(0);
  });

  it('the keys: - and = (and +) zoom about the centre, or the keyboard cursor; 0 fits; back and forward', async () => {
    const v = await make();
    v.view.setRange(Protein, 10, 30);
    expect(v.input.key({ key: '-' })).toBe(true);
    expect(range(v)).toEqual([7.5, 32.5]);
    expect(v.input.key({ key: '=' })).toBe(true);
    expect(range(v)).toEqual([10, 30]);
    // `+` is typed with Shift on many layouts
    expect(v.input.key({ key: '+', shiftKey: true })).toBe(true);
    expect(range(v)).toEqual([12, 28]);
    // about the cursor at position 14 (its cell's middle, 14.5)
    v.selection.setCursor({ section: Protein, system: '/data/p', position: 14, keyboard: true });
    v.input.key({ key: '-' });
    const [s, e] = range(v);
    expect((14.5 - s) / (e - s)).toBeCloseTo((14.5 - 12) / 16);
    expect(v.input.key({ key: '0' })).toBe(true);
    expect(range(v)).toEqual([0, 49]);
    expect(v.input.key({ key: 'Backspace' })).toBe(true);
    const back = range(v);
    expect(back[1] - back[0]).toBeCloseTo(20);
    expect(v.input.key({ key: 'ArrowRight', altKey: true })).toBe(true);
    expect(range(v)).toEqual([0, 49]);
    expect(v.input.key({ key: 'ArrowLeft', altKey: true })).toBe(true);
    expect(range(v)).toEqual(back);
    // the pan keys moved to Shift+Alt
    expect(v.input.key({ key: 'ArrowRight', altKey: true, shiftKey: true })).toBe(true);
    expect(range(v)[0]).toBeGreaterThan(back[0]);
  });

  it('the mouse’s side buttons go back and forward, through the bindings', async () => {
    const v = await make();
    run(v, 'view.focus', { section: Protein, focus: { start: 2, end: 8 } });
    expect(v.input.sideButton('back', 300, 20)).toBe(true);
    expect(v.input.sideButtonTaken('back')).toBe(true);
    expect(range(v)).toEqual([0, 49]);
    expect(v.input.sideButton('forward', 300, 20)).toBe(true);
    expect(range(v)).toEqual([2, 8]);
  });

  it('chords: modifiers before the key, `+` itself last', () => {
    expect(parseChord('Alt+ArrowLeft')).toEqual({ key: 'ArrowLeft', alt: true });
    expect(parseChord('Shift+Alt+ArrowRight')).toEqual({ key: 'ArrowRight', shift: true, alt: true });
    expect(parseChord('+')).toEqual({ key: '+' });
    expect(parseChord('Mod++')).toEqual({ key: '+', mod: true });
    expect(parseChord('Hyper+x')).toBeUndefined();
    expect(parseChord('')).toBeUndefined();
  });
});

describe('a jump’s motion', () => {
  it('draws the window between where it was and the model’s, eased, then the model’s; a hand move ends it', async () => {
    const v = await make();
    run(v, 'view.focus', { section: Protein, focus: { start: 20, end: 30 } });
    expect(v.moving).toBe(true);
    // the model is there at once
    expect(range(v)).toEqual([20, 30]);
    vi.advanceTimersByTime(v.config.view.jumpMs / 2);
    v.draw();
    const mid = v.drawnViewport!;
    expect(mid.start).toBeGreaterThan(0);
    expect(mid.start).toBeLessThan(20);
    expect(mid.end - mid.start).toBeLessThan(50);
    expect(mid.end - mid.start).toBeGreaterThan(10);
    settleMotion(v);
    expect(v.moving).toBe(false);
    expect(v.drawnViewport).toEqual({ start: 20, end: 30 });
    // a pan by hand during a motion cuts it
    run(v, 'view.back');
    expect(v.moving).toBe(true);
    v.view.zoom(Protein, 2);
    expect(v.moving).toBe(false);
  });

  it('is a cut under reduced motion, or with `view.motion: cut`', async () => {
    const v = await make();
    v.reducedMotion = true;
    run(v, 'view.focus', { section: Protein, focus: { start: 20, end: 30 } });
    expect(v.moving).toBe(false);
    v.draw();
    expect(v.drawnViewport).toEqual({ start: 20, end: 30 });
    viewer!.retire();
    const w = await make([[viewerItem('view.motion'), 'cut']]);
    run(w, 'view.focus', { section: Protein, focus: { start: 20, end: 30 } });
    expect(w.moving).toBe(false);
  });
});
