/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The `window` behavior (_spec/README.md "Behaviors", Window) in Node, over a section's viewport as a window: its body
// dragged pans it, an end resizes it about the other, a floor holds, Esc puts it back where the press found it, its
// slider's text is said when it settles (a drag's end, a step), never while it moves; its steps are commands.

import { createSeqSpec } from '@molstar/seqx-core/lang';
import { TracksLayout } from '@molstar/seqx-core/render';
import { describe, expect, it } from 'vitest';
import { windowBehavior } from '../behaviors/window.js';
import { changeKindOf } from '../commands.js';
import { viewerConfig } from '../config.js';
import { SpecModel } from '../models/spec.js';
import { ViewModel } from '../models/view.js';

const ProteinSection = '/views/main/protein';

function setup() {
  const b = createSeqSpec({ title: 'window' });
  const p = b.data.sequence('p', { alphabet: 'protein', letters: 'MEEPQSDPSV'.repeat(10) });
  b.views.view('main', (v) => v.section('protein', { of: p }).track('letters').layer('l', { type: 'letters' }));
  const config = viewerConfig({ viewport: { maxPxPerPosition: 4 } });
  const view = new ViewModel(new SpecModel(b.build()), config);
  view.setWidth(ProteinSection, 40);
  const w = windowBehavior(view, config);
  const range = () => {
    const v = view.get(ProteinSection)!;
    return [v.start, v.end];
  };
  return { view, w, range };
}

describe('the window behavior', () => {
  it('a drag on the body pans the window; Esc puts it back where the press found it', () => {
    const { view, w, range } = setup();
    view.setRange(ProteinSection, 20, 40);
    expect(w.begin(ProteinSection, 'body', 25)).toBe(true);
    expect(w.dragTo(30)).toBe(true);
    expect(range()).toEqual([25, 45]);
    // the slider says nothing while it moves
    expect(w.value.value).toBeUndefined();
    w.cancel();
    expect(range()).toEqual([20, 40]);
    expect(w.dragging).toBe(false);
  });

  it('a body drag is a pan: 2x the axis past either end, the window keeps its width and stops at the end', () => {
    const { view, w, range } = setup();
    view.setRange(ProteinSection, 20, 40);
    w.begin(ProteinSection, 'body', 30);
    w.dragTo(30 + 200);
    expect(range()).toEqual([80, 100]);
    w.dragTo(30 - 200);
    expect(range()).toEqual([0, 20]);
    w.end();
  });

  it('a drag on an end resizes about the other; the floor holds; the end of a drag says where it settled', () => {
    const { view, w, range } = setup();
    view.setRange(ProteinSection, 20, 40);
    w.begin(ProteinSection, 'end', 40);
    w.dragTo(60);
    expect(range()).toEqual([20, 60]);
    // past the other end: no narrower than the floor (40 px at 4 px a position: 10), the start held
    w.dragTo(5);
    expect(range()).toEqual([20, 30]);
    w.end();
    expect(w.value.value).toMatchObject({ section: ProteinSection, start: 21, end: 30, text: '21 to 30 of 100' });
    w.begin(ProteinSection, 'start', 20);
    w.dragTo(0);
    expect(range()).toEqual([0, 30]);
    w.end();
  });

  it('its steps are commands that navigate: pan a step, grow and shrink about its middle', async () => {
    const { view, w, range } = setup();
    view.setRange(ProteinSection, 20, 40);
    const byId = new Map(w.commands.map((c) => [c.id, c]));
    const signal = new AbortController().signal;
    for (const c of w.commands) expect(changeKindOf(c.id)).toBe('navigate');
    await byId.get('window.pan-on')!.run({}, { section: ProteinSection }, signal);
    expect(range()).toEqual([22, 42]);
    await byId.get('window.grow')!.run({}, { section: ProteinSection }, signal);
    expect(range()[1] - range()[0]).toBeCloseTo(25, 6);
    await byId.get('window.shrink')!.run({}, { section: ProteinSection }, signal);
    expect(range()[1] - range()[0]).toBeCloseTo(20, 6);
    expect(w.value.value?.text).toMatch(/ of 100$/);
    expect(byId.get('window.grow')!.when!({ section: undefined })).toBe('No window to move');
  });

  it('on a cyclic axis a body drag over the origin is a short step on where the layout wraps, else the pointer’s place', () => {
    const ring = '/views/main/ring';
    const setup = (wraps: boolean) => {
      const b = createSeqSpec({ title: 'ring' });
      const c = b.data.sequence('c', { alphabet: 'dna', letters: 'ACGTACGTAC'.repeat(10), cyclic: true });
      b.views.view('main', (v) => v.section('ring', { of: c }).track('letters').layer('l', { type: 'letters' }));
      const config = viewerConfig({ viewport: { maxPxPerPosition: 4 } });
      const view = new ViewModel(new SpecModel(b.build()), config);
      // the tracks layout wraps a cyclic window; a layout that does not (text's) said so here
      view.setLayoutSource({ placement: () => undefined, layouts: { tracks: { ...TracksLayout, wraps } } });
      view.setWidth(ring, 40);
      view.setRange(ring, 80, 100);
      return { view, w: windowBehavior(view, config), range: () => [view.get(ring)!.start, view.get(ring)!.end] };
    };
    // the pointer at 95, then 99, then over the origin at 5: ten units on, not ninety back
    const r = setup(true);
    expect(r.w.begin(ring, 'body', 95)).toBe(true);
    r.w.dragTo(99);
    expect(r.range()).toEqual([84, 104]);
    r.w.dragTo(5);
    expect(r.range()).toEqual([90, 110]);
    // and back over it
    r.w.dragTo(90);
    expect(r.range()).toEqual([75, 95]);
    r.w.end();
    // a layout that draws no window over the origin: the unit under the pointer as it is (one rule: `view.wraps`)
    const t = setup(false);
    t.w.begin(ring, 'body', 95);
    t.w.dragTo(5);
    expect(t.range()).toEqual([0, 20]);
    t.w.dragTo(90);
    expect(t.range()).toEqual([75, 95]);
  });

  it('a window over a cyclic axis’s end is said as the arc over it, from its start to its end past the total', () => {
    const ring = '/views/main/ring';
    const b = createSeqSpec({ title: 'ring' });
    const c = b.data.sequence('c', { alphabet: 'dna', letters: 'ACGTACGTAC'.repeat(10), cyclic: true });
    b.views.view('main', (v) => v.section('ring', { of: c }).track('letters').layer('l', { type: 'letters' }));
    const config = viewerConfig({ viewport: { maxPxPerPosition: 4 } });
    const view = new ViewModel(new SpecModel(b.build()), config);
    view.setLayoutSource({ placement: () => undefined, layouts: { tracks: TracksLayout } });
    view.setWidth(ring, 40);
    view.setRange(ring, 90, 110);
    const w = windowBehavior(view, config);
    w.begin(ring, 'body', 95);
    w.dragTo(96);
    w.end();
    expect(w.value.value).toMatchObject({ start: 92, end: 11, text: '92 to 11 of 100' });
  });
});
