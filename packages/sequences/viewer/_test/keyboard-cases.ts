/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// IN4's key repeat, for its heap guard and its quick twin (keyboard.heap.ts, keyboard.quick.ts): case 3 focused with
// labels off (`family`; a label is a string per hovered element by design), one `KeyInput` for every press (`press`),
// and the steps: presses that leave the cursor at the end of the bottom row (`staying`, the menu key among them),
// `cursorRect` on data, a row and a header (`rects`), the zoom keys (`zooming`), Left and Right along a row (`along`),
// Up and Down, Home and End (`upDown`), and Down and Up over a million-row run (`millionRows`).

import { createSeqViewSpec } from '@molstar/seq-core/builder';
import { case3 } from '../../core/_test/format/cases.js';
import { Inputs, type KeyCursor, Modifier } from '../src/index.js';
import { cursorRect } from '../src/keyboard.js';
import { type Driven, driven } from './input-fixtures.js';
import { lazyLoaders } from './layout-fixtures.js';

const S = 'family.main';

export async function family() {
  const t = await driven(case3(), { size: { [S]: [1000, 150] }, labels: false });
  t.focus(S);
  return t;
}

/** One `KeyInput` for all presses, as an engine's (the fixtures' `key` makes a literal per call). */
const k = Inputs.key();
export function press(t: Driven, key: string, shift = false, section = S) {
  k.section = t.ctx.layout.sectionIndex(section);
  k.key = key;
  k.modifiers = shift ? Modifier.shift : 0;
  return t.ctx.input.key(k);
}

/** Presses that move the cursor: Right along a row across the header, then Left back. */
export const along = (t: Driven) => (i: number) => void press(t, i % 66 < 33 ? 'ArrowRight' : 'ArrowLeft');

/** The cursor moved to the end of the bottom row, where `staying`'s presses leave it. */
export function toBottomEnd(t: Driven) {
  for (const key of ['End', 'ArrowRight', 'ArrowDown', 'ArrowDown', 'ArrowDown', 'ArrowDown', 'ArrowDown'])
    press(t, key);
}
/** Presses that leave the cursor at the end of the bottom row: Right, Down, the menu key, End. */
export const staying = (t: Driven) => (i: number) =>
  void press(t, i % 4 === 0 ? 'ArrowRight' : i % 4 === 1 ? 'ArrowDown' : i % 4 === 2 ? 'F10' : 'End', i % 4 === 2);

/** `cursorRect`, read per frame by the ring and the tooltip: on data, a row and a header. */
export function rects(t: Driven) {
  const out = new Float64Array(6);
  const at: KeyCursor[] = [];
  for (const key of ['ArrowRight', 'ArrowDown', 'Home']) {
    press(t, key);
    at.push(t.ctx.ui.keyCursor.value as KeyCursor);
  }
  return (i: number) => void cursorRect(t.ctx, at[i % 3], out);
}

/** The zoom keys, anchored at the cursor (after a Right). */
export function zooming(t: Driven) {
  press(t, 'ArrowRight');
  return (i: number) => void press(t, i & 1 ? '+' : '-');
}

/** Up and Down along the tracks and rows, and Home and End. */
export const upDown = (t: Driven) => (i: number) =>
  void press(t, ['ArrowDown', 'ArrowDown', 'End', 'ArrowUp', 'ArrowUp', 'Home'][i % 6]);

/** The rows of a million-row run, focused: Down and Up, 200 each. */
export async function millionRows() {
  const [uri, n] = ['https://example.org/huge', 1_000_000];
  const b = createSeqViewSpec();
  const aln = b.data.alignment({
    ref: 'huge',
    row_count: n,
    column_count: 200,
    source: { uri, protocol: 'seqstar-tiles/1' },
  });
  b.view({ ref: 'v', label: 'V' }).block({ alignment: aln }).track({ ref: 'v.rows' }).rows();
  const t = await driven(b.getDocument(), { loaders: lazyLoaders(uri, n, 200), labels: false });
  t.focus('v.main');
  return (i: number) => void press(t, i % 400 < 200 ? 'ArrowDown' : 'ArrowUp', false, 'v.main');
}
