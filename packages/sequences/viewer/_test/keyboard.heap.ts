/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Heap guards for IN4's key repeat (seqspec:viewer#behaviors; the cursor is a published value, `ui.keyCursor`): on
// `heapScaling` (a step is one key press; labels off, since a label is a string per hovered element by design, bounded
// by IN2's churn guard). A press that leaves the cursor where it is allocates nothing: the edges of a row and of the
// rows, the zoom keys, the menu key. A press that moves it makes the one small record the state publishes and nothing
// else: a bound per press, whatever the number of positions, rows and blocks (a million-row run included). The positive
// controls, one per criterion: the presses that leave the cursor allocating an object each (`leaking`), primed and
// measured alike, are caught as not flat, and the presses that move it, allocating likewise, by the per-press bound.

import { createSeqViewSpec } from '@molstar/seq-core/builder';
import { describe, expect, it } from 'vitest';
import { case3 } from '../../core/_test/format/cases.js';
import { Inputs, type KeyCursor, Modifier } from '../src/index.js';
import { cursorRect } from '../src/keyboard.js';
import { controlCaught, describeScaling, heapScaling, isFlat, leaking } from './heap.js';
import { type Driven, driven } from './input-fixtures.js';
import { lazyLoaders } from './layout-fixtures.js';

const N = 20_000;
/** One published cursor (88 bytes; 93 measured, 2026-10-11) is all a press may make. */
const PER_PRESS = 112;

/** At least 300,000 presses, then until three chunks in a row allocate under 2 KB: V8's lower tiers box doubles. */
function prime(step: (i: number) => void) {
  const gc = (globalThis as { gc?: () => void }).gc;
  let k = 0;
  for (let chunk = 0, quiet = 0; chunk < 100 && (quiet < 3 || chunk < 15); chunk++) {
    gc?.();
    const before = process.memoryUsage().heapUsed;
    for (let j = 0; j < 20_000; j++) step(k++);
    quiet = process.memoryUsage().heapUsed - before < 2048 ? quiet + 1 : 0;
  }
}
const guard = (step: (i: number) => void) => {
  prime(step);
  const r = heapScaling(N, step);
  console.info(`key repeat: ${(r.scaled / (3 * N)).toFixed(1)} bytes a press`);
  return r;
};

async function family() {
  const t = await driven(case3(), { size: { 'family.main': [1000, 150] }, labels: false });
  t.focus('family.main');
  return t;
}
const S = 'family.main';
/** One `KeyInput` for all presses, as an engine's (the fixtures' `key` makes a literal per call). */
const k = Inputs.key();
function press(t: Driven, key: string, shift = false, section = S) {
  k.section = t.ctx.layout.sectionIndex(section);
  k.key = key;
  k.modifiers = shift ? Modifier.shift : 0;
  return t.ctx.input.key(k);
}

/** Presses that move the cursor: Right along a row across the header, then Left back. */
const along = (t: Driven) => (i: number) => void press(t, i % 66 < 33 ? 'ArrowRight' : 'ArrowLeft');
/** Presses that leave the cursor at the end of the bottom row: Right, Down, the menu key, End. */
const staying = (t: Driven) => (i: number) =>
  void press(t, i % 4 === 0 ? 'ArrowRight' : i % 4 === 1 ? 'ArrowDown' : i % 4 === 2 ? 'F10' : 'End', i % 4 === 2);

describe("IN4's key repeat", { timeout: 300_000 }, () => {
  it('a press that leaves the cursor where it is allocates nothing: the end of a row, the bottom, the top, the menu key', async () => {
    const t = await family();
    const keys = ['End', 'ArrowRight', 'ArrowDown', 'ArrowDown', 'ArrowDown', 'ArrowDown', 'ArrowDown'];
    for (const k of keys) press(t, k);
    const at = t.ctx.ui.keyCursor.value;
    const r = guard(staying(t));
    expect(isFlat(r), describeScaling(r, N)).toBe(true);
    expect(t.ctx.ui.keyCursor.value).toBe(at);
  });

  it('cursorRect, read per frame by the ring and the tooltip, allocates nothing: on data, a row and a header', async () => {
    const t = await family();
    const out = new Float64Array(6);
    const at: KeyCursor[] = [];
    for (const k of ['ArrowRight', 'ArrowDown', 'Home']) {
      press(t, k);
      at.push(t.ctx.ui.keyCursor.value as KeyCursor);
    }
    const r = guard((i) => void cursorRect(t.ctx, at[i % 3], out));
    expect(isFlat(r), describeScaling(r, N)).toBe(true);
  });

  it('zoom keys, anchored at the cursor, allocate nothing; the cursor is not republished', async () => {
    const t = await family();
    press(t, 'ArrowRight');
    const at = t.ctx.ui.keyCursor.value;
    const r = guard((i) => void press(t, i & 1 ? '+' : '-'));
    expect(isFlat(r), describeScaling(r, N)).toBe(true);
    expect(t.ctx.ui.keyCursor.value).toBe(at);
  });

  it('Left and Right along a row, across the header and back: a bounded record per press, flat in the count', async () => {
    const t = await family();
    const r = guard(along(t));
    expect(r.scaled / (3 * N), describeScaling(r, N)).toBeLessThan(PER_PRESS);
  });

  it('Up and Down along the tracks and rows, and Home and End: the same bound', async () => {
    const t = await family();
    const r = guard((i) => void press(t, ['ArrowDown', 'ArrowDown', 'End', 'ArrowUp', 'ArrowUp', 'Home'][i % 6]));
    expect(r.scaled / (3 * N), describeScaling(r, N)).toBeLessThan(PER_PRESS);
  });

  it('the rows of a million-row run: Down and Up keep the bound, however far the cursor is from the first row', async () => {
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
    const down = (i: number) => void press(t, i % 400 < 200 ? 'ArrowDown' : 'ArrowUp', false, 'v.main');
    const r = guard(down);
    expect(r.scaled / (3 * N), describeScaling(r, N)).toBeLessThan(PER_PRESS);
  });

  // The controls last: cases in one file share V8's feedback.
  it('the same presses allocating an object each are caught (positive control)', async () => {
    const t = await family();
    for (const k of ['End', 'ArrowRight', 'ArrowDown', 'ArrowDown', 'ArrowDown', 'ArrowDown', 'ArrowDown']) press(t, k);
    expect(controlCaught(() => guard(leaking(staying(t))))).toBe(true);
  });

  it('presses that move the cursor, allocating an object each, are caught by the per-press bound (positive control)', async () => {
    const t = await family();
    expect(controlCaught(() => guard(leaking(along(t))).scaled / (3 * N), PER_PRESS)).toBe(true);
  });
});
