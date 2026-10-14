/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Feature labels' geometry (labels.ts), without a DOM: a track's view as its xScale maps it, a label's box at a zoom,
// and the ink on a fill.

import { describe, expect, it } from 'vitest';
import {
  axisView,
  cutText,
  INK,
  inkOn,
  type LabelBox,
  type NightingaleAxis,
  placeLabel,
  positionAt,
  watchDimensions,
  xOf,
} from '../src/labels.js';

const view = (...a: [number | undefined, number | undefined, number, number, number, number]) =>
  axisView({ start: 0, end: 0, left: 0, right: 0 }, ...a);
const box = (): LabelBox => ({ mode: 'none', left: 0, width: 0 });

describe('the view', () => {
  it("is Nightingale's xScale: [start, end + 1] on the width less the margins; -1 or unset ends at the length", () => {
    const v = view(1, -1, 100, 220, 10, 10);
    expect(v).toEqual({ start: 1, end: 100, left: 10, right: 210 });
    expect([xOf(v, 1), xOf(v, 51), xOf(v, 101)]).toEqual([10, 110, 210]);
    expect(view(undefined, undefined, 50, 100, 0, 0)).toMatchObject({ start: 1, end: 50 });
    // Zoomed to 41..60: twenty positions over 200 px.
    const z = view(41, 60, 100, 220, 10, 10);
    expect([xOf(z, 41), xOf(z, 61)]).toEqual([10, 210]);
  });
});

describe('a position under the pointer', () => {
  it('is the position whose span holds x, at any zoom; none in the gutters', () => {
    const v = view(1, -1, 100, 220, 10, 10); // 2 px a position from x = 10
    expect([positionAt(v, 10), positionAt(v, 11.9), positionAt(v, 12), positionAt(v, 209.9)]).toEqual([1, 1, 2, 100]);
    expect([positionAt(v, 9.9), positionAt(v, 210), positionAt(v, Number.NaN)]).toEqual([-1, -1, -1]);
    const z = view(41, 60, 100, 220, 10, 10); // 10 px a position
    expect([positionAt(z, 10), positionAt(z, 25), positionAt(z, 209)]).toEqual([41, 42, 60]);
    // Both ways agree: a position's left edge maps back to it.
    for (const p of [41, 50, 60]) expect(positionAt(z, xOf(z, p))).toBe(p);
  });
});

describe("a label's box", () => {
  // 100 positions on 200 px: 2 px a position.
  const v: NightingaleAxis = view(1, 100, 100, 200, 0, 0);
  // A label `text` px wide whose narrowest cut ("Tra…") is 20 px.
  /** A label's fit: its text's width, its cut's, its feature's height. */
  const label = (text: number) => [text, Math.min(text, 20), 15] as const;

  it("is the feature's visible span, the text whole when it fits with its padding", () => {
    expect(placeLabel(box(), v, [{ start: 1, end: 60 }], ...label(40))).toEqual({ mode: 'full', left: 0, width: 120 });
    expect(placeLabel(box(), v, [{ start: 61, end: 100 }], ...label(74))).toEqual({
      mode: 'full',
      left: 120,
      width: 80,
    });
  });

  it('is cut short where it does not fit, down to its first characters and an ellipsis; hidden below that', () => {
    expect(placeLabel(box(), v, [{ start: 1, end: 20 }], ...label(60)).mode).toBe('cut'); // 40 px
    expect(placeLabel(box(), v, [{ start: 1, end: 13 }], ...label(60)).mode).toBe('cut'); // 26 px: 20 and the padding
    expect(placeLabel(box(), v, [{ start: 1, end: 12 }], ...label(60)).mode).toBe('none'); // 24 px
    // A short label is never cut: whole, or hidden.
    expect(placeLabel(box(), v, [{ start: 1, end: 12 }], 30, 30, 15).mode).toBe('none');
    expect(placeLabel(box(), v, [{ start: 1, end: 60 }], 10, 10, 6).mode).toBe('none');
    expect([cutText('Transactivation 1'), cutText('VH'), cutText('CDRs')]).toEqual(['Tra…', null, 'CDR…']);
  });

  it('follows zoom and pan: clipped to the viewport, in the widest visible fragment', () => {
    const zoomed = view(41, 60, 100, 200, 0, 0); // 10 px a position
    const out = box();
    // A block over 1..50 shows 41..50 at the left edge.
    expect(placeLabel(out, zoomed, [{ start: 1, end: 50 }], ...label(30))).toEqual({
      mode: 'full',
      left: 0,
      width: 100,
    });
    // The same box object, written again (no allocation per frame).
    expect(placeLabel(out, zoomed, [{ start: 70, end: 80 }], ...label(30))).toBe(out);
    expect(out.mode).toBe('none');
    const two = [
      { start: 30, end: 42 },
      { start: 50, end: 57 },
    ];
    expect(placeLabel(box(), zoomed, two, ...label(50))).toEqual({ mode: 'full', left: 90, width: 80 });
  });
});

describe('the ink', () => {
  it('reads on the fill: dark on light fills, light on dark ones, gray for unknown colors', () => {
    expect(inkOn('#e8e8e8')).toBe(INK.dark);
    expect(inkOn('#fdae61')).toBe(INK.dark);
    expect(inkOn('#1f77b4')).toBe(INK.dark); // a mid blue: dark has the higher contrast
    expect(inkOn('#253494')).toBe(INK.light);
    expect(inkOn('#000', 1)).toBe(INK.light);
    expect(inkOn('#000', 0.2)).toBe(INK.dark); // mostly the white ground
    expect(inkOn('rebeccapurple')).toBe(inkOn('#808080'));
  });
});

describe("a track's dimensions", () => {
  // Nightingale sets `width` from the measured size a frame or more after the layout, with no attribute or event.
  it("are watched through the element's own onDimensionsChange, its work kept; stopping restores it", () => {
    const calls: string[] = [];
    const proto = {
      onDimensionsChange() {
        calls.push('element');
      },
    };
    const node: { onDimensionsChange?(): void } = Object.create(proto);
    const stop = watchDimensions(node, () => calls.push('watcher'));
    node.onDimensionsChange?.();
    expect(calls).toEqual(['element', 'watcher']);
    stop();
    expect(Object.hasOwn(node, 'onDimensionsChange')).toBe(false);
    node.onDimensionsChange?.();
    expect(calls).toEqual(['element', 'watcher', 'element']);
  });

  it('are watched on an element without the hook, and an own hook is put back', () => {
    let seen = 0;
    const bare: { onDimensionsChange?(): void } = {};
    const stopBare = watchDimensions(bare, () => seen++);
    bare.onDimensionsChange?.();
    expect(seen).toBe(1);
    stopBare();
    expect(bare.onDimensionsChange).toBeUndefined();
    const own = { onDimensionsChange: () => seen++ };
    const stopOwn = watchDimensions(own, () => seen++);
    own.onDimensionsChange();
    expect(seen).toBe(3);
    const hook = own.onDimensionsChange;
    stopOwn();
    expect(own.onDimensionsChange).not.toBe(hook);
    own.onDimensionsChange();
    expect(seen).toBe(4);
  });
});
