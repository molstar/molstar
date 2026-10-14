/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The batches (seqspec:viewer#the-render-backend): they grow, keep their room across resets, carry colors per item,
// and append ranges (the recording's stores and its replay).

import { describe, expect, it } from "vitest";
import {
  GlyphBatch,
  NO_COLOR,
  PathBatch,
  PathMode,
  PathVerb,
  PathVerbCoords,
  RectBatch,
  SectorBatch,
  SkinToken,
  TextAlign,
  TextBaseline,
  TextBatch,
} from "../../src/draw/api.ts";

describe("RectBatch", () => {
  it("grows past its capacity with every item kept, and keeps its arrays across resets", () => {
    const b = new RectBatch(4);
    b.begin(0x112233, 128, 2);
    for (let i = 0; i < 10; i++) b.push(i, i * 2, 3, 4);
    expect(b.n).toBe(10);
    expect(b.x.length).toBeGreaterThanOrEqual(10);
    expect([...b.y.subarray(0, 10)]).toEqual([0, 2, 4, 6, 8, 10, 12, 14, 16, 18]);
    expect(b.color[9]).toBe(0x112233);
    expect(b.alpha[9]).toBe(128);
    expect(b.stroke[9]).toBe(2);
    const x = b.x;
    b.reset();
    b.begin(SkinToken.grid);
    for (let i = 0; i < 10; i++) b.push(i, 0, 1, 1);
    expect(b.x).toBe(x);
    expect(b.color[0]).toBe(SkinToken.grid);
    expect(b.alpha[0]).toBe(255);
    expect(b.stroke[0]).toBe(0);
  });

  it("colors follow `begin`; appended ranges carry their own colors", () => {
    const a = new RectBatch(2);
    a.begin(1);
    a.push(0, 0, 1, 1);
    a.begin(2, 50);
    a.push(1, 0, 1, 1);
    a.begin(NO_COLOR);
    a.push(2, 0, 1, 1);
    const b = new RectBatch(1);
    b.append(a, 1, 2);
    b.append(a, 0, 1);
    expect(b.n).toBe(3);
    expect([...b.x.subarray(0, 3)]).toEqual([1, 2, 0]);
    expect([...b.color.subarray(0, 3)]).toEqual([2, NO_COLOR, 1]);
    expect([...b.alpha.subarray(0, 3)]).toEqual([50, 255, 255]);
  });
});

describe("double batches", () => {
  it("write through the slot `add` returns, and grow", () => {
    const s = new SectorBatch(1);
    s.begin(7, 9);
    for (let i = 0; i < 5; i++) s.push(i, 0, 1, 2, 0.5, 1.5);
    expect(s.n).toBe(5);
    expect(s.cx[4]).toBe(4);
    expect(s.a1[4]).toBe(1.5);
    expect(s.alpha[4]).toBe(9);

    const g = new GlyphBatch(1);
    g.begin(3);
    const k = g.add(5);
    g.x[k] = 1.25;
    g.push(2.5, 3.5, 4, 2);
    expect(g.n).toBe(2);
    expect([g.glyph[0], g.glyph[1], g.x[0], g.radius[1], g.color[1]]).toEqual([5, 2, 1.25, 4, 3]);
  });

  it("text holds references to the strings it was given, one style per batch", () => {
    const t = new TextBatch(1);
    t.style(3, TextAlign.right, TextBaseline.top, Math.PI / 2);
    const label = "Kinase domain";
    t.begin(0xff0000, 100);
    t.push(1.5, 2.5, label);
    t.begin(SkinToken.inkQuiet);
    t.push(3, 4, "x");
    expect(t.n).toBe(2);
    expect(t.text[0]).toBe(label);
    expect([t.font, t.align, t.baseline, t.angle]).toEqual([3, TextAlign.right, TextBaseline.top, Math.PI / 2]);
    expect([t.color[0], t.alpha[0], t.color[1], t.alpha[1]]).toEqual([0xff0000, 100, SkinToken.inkQuiet, 255]);
    const copy = new TextBatch(1);
    copy.append(t, 1, 1);
    expect(copy.text[0]).toBe("x");
  });
});

describe("PathBatch", () => {
  it("keeps subpaths with their styles, verbs and coordinates, and appends a range with its offsets", () => {
    const p = new PathBatch(1);
    p.begin(1, 255, PathMode.stroke, 3);
    p.moveTo(0, 0);
    p.lineTo(10, 0);
    p.begin(2);
    p.moveTo(1, 1);
    p.quadTo(2, 2, 3, 3);
    p.close();
    p.begin(3, 128);
    const o = p.verb(PathVerb.roundRect);
    p.coords.set([5, 6, 7, 8, 2], o);
    expect(p.n).toBe(3);
    expect(p.nv).toBe(6);
    expect(p.nc).toBe(4 + 6 + PathVerbCoords[PathVerb.roundRect]);
    expect([...p.start.subarray(0, 3)]).toEqual([0, 2, 5]);
    expect([...p.cstart.subarray(0, 3)]).toEqual([0, 4, 10]);
    expect([p.mode[0], p.width[0], p.mode[1], p.width[1]]).toEqual([PathMode.stroke, 3, PathMode.fill, 2]);

    const q = new PathBatch(1);
    q.begin(9);
    q.moveTo(-1, -1);
    q.append(p, 1, 2);
    expect(q.n).toBe(3);
    expect([...q.start.subarray(0, 3)]).toEqual([0, 1, 4]);
    expect([...q.cstart.subarray(0, 3)]).toEqual([0, 2, 8]);
    expect([...q.verbs.subarray(0, q.nv)]).toEqual([
      PathVerb.move,
      PathVerb.move,
      PathVerb.quad,
      PathVerb.close,
      PathVerb.roundRect,
    ]);
    expect([...q.coords.subarray(0, q.nc)]).toEqual([-1, -1, 1, 1, 2, 2, 3, 3, 5, 6, 7, 8, 2]);
    expect([...q.color.subarray(0, 3)]).toEqual([9, 2, 3]);
  });
});
