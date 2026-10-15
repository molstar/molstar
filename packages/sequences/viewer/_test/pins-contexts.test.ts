/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// VL-09: a boundary pin over letters reads its font's ink from its own context's text manager. Two viewers whose
// letters' fonts share an id (the same order of fonts asked) at different sizes each pin by their own font, and a
// text manager's `invalidate` (a web font loaded) drops what was read.

import { describe, expect, it } from 'vitest';
import { Dirty, SectionComposer } from '../src/draw/index.js';
import type { TextManager, TextMeasurer } from '../src/text.js';
import { pinsDoc } from './boundary-pins-fixtures.js';
import { recording } from './composer-fixtures.js';
import { lit } from './overlay-fixtures.js';

/** The pins' paths (triangles) of a frame at 16 px a position, as coordinate lists. */
async function viewer(maxLetterPx: number) {
  const t = await lit(pinsDoc(1), { width: 900, height: 400, dpr: 2, render: { maxLetterPx } });
  t.vp?.setView(16, 0, 'api');
  // A fresh composer each time: the whole frame drawn.
  return Object.assign(
    () => {
      const c = new SectionComposer(t.ctx, t.section);
      const r = recording();
      c.frame(r, Dirty.all);
      c.dispose();
      return Array.from(r.pathItems.coords.subarray(0, r.pathItems.nc));
    },
    { text: t.ctx.text },
  );
}

describe("VL-09: the pin's ink per text manager", () => {
  it("two viewers, one font id at two sizes: each viewer's pins as when drawn alone", async () => {
    const small = await viewer(9);
    const alone = small();
    const large = await viewer(16);
    large();
    expect(small()).toEqual(alone);
    expect(alone).not.toEqual(large());
  });

  it('a new measurer (its version bumped) is read again: a smaller capital, a taller triangle (LC-09)', async () => {
    const pins = await viewer(16);
    const before = pins();
    const text = pins.text as TextManager;
    const m = text.port;
    const smaller: TextMeasurer = { measure: (f: string, s: string) => m.measure(f, s) };
    text.setMeasurer(Object.assign(smaller, { metrics: () => ({ ascent: 0.9, descent: 0.2, cap: 0.5 }) }));
    const after = pins();
    expect(after.length).toBe(before.length);
    expect(after).not.toEqual(before);
  });
});
