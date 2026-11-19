/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The `window` behavior (_spec/README.md "Behaviors", Window) in Node, over a section's viewport as a window: its body
// dragged pans it, an end resizes it about the other, a floor holds, Esc puts it back where the press found it, its
// slider's text is said when it settles (a drag's end, a step), never while it moves; its steps are commands.

import { createSeqSpec } from '@molstar/seqx-core/lang';
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
});
