/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The rows rail's gesture on a band whose rows or lanes change under a press (the rail is the frames' live object): the
// press starts again where the pointer is, within what there is now, and a cancel puts back what the band held.

import { createSeqSpec } from '@molstar/seqx-core/lang';
import { describe, expect, it } from 'vitest';
import { RowsRailGesture, type RowsRailWindow } from '../behaviors/rows-rail.js';
import { viewerConfig } from '../config.js';
import { SpecModel } from '../models/spec.js';
import { ViewModel } from '../models/view.js';

const Track = '/views/main/s/t';

function setup() {
  const b = createSeqSpec({ title: 'rail' });
  const p = b.data.sequence('p', { alphabet: 'protein', letters: 'MEEPQSDPSV'.repeat(10) });
  b.views.view('main', (v) => v.section('s', { of: p }).track('t').layer('l', { type: 'letters' }));
  const view = new ViewModel(new SpecModel(b.build()), viewerConfig());
  return { view, g: new RowsRailGesture(view, viewerConfig()) };
}

describe('the rows rail’s gesture over a changing band', () => {
  it('lanes that grow between the press and a move: the window follows the pointer from there, never past the new end', () => {
    const { view, g } = setup();
    // 40 lanes, 10 in view, the top at 5
    view.setBand(Track, { bandTop: 5 });
    const rail: { -readonly [K in keyof RowsRailWindow]: RowsRailWindow[K] } = {
      start: 5,
      end: 15,
      total: 40,
      width: 200,
      lanes: true,
    };
    expect(g.press(Track, 'body', 6, rail)).toBe(true);
    // the lanes grow to 100 under the press: the window is as long, the unit twice as short on the rail
    rail.total = 100;
    g.dragTo(6);
    expect(view.band(Track)?.bandTop).toBe(5);
    g.dragTo(16);
    expect(view.band(Track)?.bandTop).toBe(15);
    // and shrink to 20: the top never past the lanes there are
    rail.total = 20;
    g.dragTo(500);
    expect(view.band(Track)?.bandTop).toBe(10);
  });

  it('a cancel puts a band that had none back to none, and one that had a top back to it', () => {
    const { view, g } = setup();
    const rail = { start: 0, end: 10, total: 40, width: 200, lanes: true };
    g.press(Track, 'body', 1, rail);
    g.dragTo(9);
    expect(view.band(Track)?.bandTop).toBe(8);
    g.cancel();
    expect(view.band(Track)).toBeUndefined();
    view.setBand(Track, { bandTop: 3 });
    g.press(Track, 'body', 4, { ...rail, start: 3, end: 13 });
    g.dragTo(20);
    g.cancel();
    expect(view.band(Track)?.bandTop).toBe(3);
  });
});
