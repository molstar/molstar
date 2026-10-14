/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A `ranges` boundary over its track's letters (seqspec:viewer#the-render-backend): a triangle pinned to the row's top
// edge pointing down at the boundary, and a hairline down it, never over a letter's ink (the advance wide, a capital's
// ink high, a descender's room below) at any zoom; where no room is left above the text, the hairline alone; at a
// block's ends inside the data area. Below the letters LOD, its glyph on the lane as before. It hits within a few px of
// the boundary over the whole row.

import { describe, expect, it } from "vitest";
import { Hit, HoverKind } from "../src/index.js";
import { Op, type RecordingBackend } from "../src/testing/index.js";
import { pinsDoc as doc, S } from "./boundary-pins-fixtures.js";
import { recording } from "./composer-fixtures.js";
import { lit } from "./overlay-fixtures.js";

type Box = [number, number, number, number];
const meets = (a: Box, b: Box) => a[0] < b[2] && b[0] < a[2] && a[1] < b[3] && b[1] < a[3];

/** The marks' boxes (device px: x0, y0, x1, y1), triangles and hairlines, by the mark's color; and the triangles. */
function markBoxes(r: RecordingBackend, color: number) {
  const boxes: Box[] = [];
  let triangles = 0;
  const p = r.pathItems;
  for (let i = 0; i < p.n; i++) {
    if (p.color[i] !== color) continue;
    const c0 = p.cstart[i];
    const c1 = i + 1 < p.n ? p.cstart[i + 1] : p.nc;
    const xs: number[] = [];
    const ys: number[] = [];
    for (let c = c0; c < c1; c += 2) {
      xs.push(p.coords[c]);
      ys.push(p.coords[c + 1]);
    }
    boxes.push([Math.min(...xs), Math.min(...ys), Math.max(...xs), Math.max(...ys)]);
    triangles++;
  }
  const q = r.rectItems;
  for (let i = 0; i < q.n; i++)
    if (q.color[i] === color) boxes.push([q.x[i], q.y[i], q.x[i] + q.w[i], q.y[i] + q.h[i]]);
  return { boxes, triangles };
}

/** The letters' ink boxes (device px): the advance wide, a capital's ink above the baseline, a descender's below. */
function inkBoxes(r: RecordingBackend, ink: (font: number) => [number, number]): Box[] {
  const out: Box[] = [];
  const t = r.textItems;
  for (let k = 0; k < r.opCount; k++) {
    if (r.ops[4 * k] !== Op.text) continue;
    const font = r.args[r.ops[4 * k + 3]];
    const [adv, cap] = ink(font);
    for (let i = r.ops[4 * k + 1]; i < r.ops[4 * k + 1] + r.ops[4 * k + 2]; i++) {
      if (t.text[i].length !== 1 || !/[A-Z]/.test(t.text[i])) continue;
      out.push([t.x[i] - adv / 2, t.y[i] - cap, t.x[i] + adv / 2, t.y[i] + Math.ceil(cap / 4)]);
    }
  }
  return out;
}

async function open(height: number, dpr: number) {
  const t = await lit(doc(height), { width: 900, height: 400, dpr });
  // The mark's packed color, from the frame's own paths and rects: the color no letter or cell is drawn in.
  const ink = (font: number): [number, number] => [
    t.ctx.text.cell(font as never, null) * dpr,
    t.ctx.text.ink(font as never, Math.round(dpr * 64))[0],
  ];
  return { ...t, ink };
}

/** The color a frame draws its triangles in (the mark's). */
function triangleColor(r: RecordingBackend): number {
  const p = r.pathItems;
  for (let i = 0; i < p.n; i++) if (p.mode[i] === 0) return p.color[i];
  return -1;
}

describe("a boundary over its track's letters", () => {
  for (const dpr of [1, 2])
    for (const cell of [9, 16, 30])
      it(`never covers a letter's ink: ${cell} px a position, DPR ${dpr}`, async () => {
        const t = await open(1, dpr);
        // The view from the graft's boundary on: the start, the middle (20) and, at 9 px, the end in view.
        t.vp?.setView(cell, 0, "api");
        const r = t.frame(recording());
        const color = triangleColor(r);

        expect(color).not.toBe(-1);
        const { boxes, triangles } = markBoxes(r, color);
        const inks = inkBoxes(r, t.ink);
        expect(inks.length).toBeGreaterThan(20);
        // Both rows' marks at 0 and 20 at least: a triangle and a hairline each.
        expect(triangles).toBeGreaterThanOrEqual(4);
        for (const m of boxes) for (const k of inks) expect(meets(m, k), `${m} over ${k}`).toBe(false);
      });

  it("where no room is left above the text, the hairline alone, still clear of the ink", async () => {
    const t = await open(0.55, 2);
    t.vp?.setView(16, 0, "api");
    const r = t.frame(recording());
    const inks = inkBoxes(r, t.ink);
    expect(inks.length).toBeGreaterThan(20);
    // No triangle; the mark's hairlines (the only rects of 2 device px wide over the rows) never meet the ink.
    expect(Array.from(r.pathItems.mode.subarray(0, r.pathItems.n)).filter((m) => m === 0).length).toBe(0);
    const q = r.rectItems;
    const lines: Box[] = [];
    for (let i = 0; i < q.n; i++) if (q.w[i] === 2) lines.push([q.x[i], q.y[i], q.x[i] + q.w[i], q.y[i] + q.h[i]]);
    expect(lines.length).toBeGreaterThan(0);
    for (const m of lines) for (const k of inks) expect(meets(m, k), `${m} over ${k}`).toBe(false);
  });

  it("at a block's ends the mark sits inside the data area, and it hits within a few px over the row", async () => {
    const t = await open(1, 2);
    // Fitted: the block spans the data area, its first and last boundaries on its edges.
    const r = t.frame(recording());
    const l = t.layout(S);
    const { boxes, triangles } = markBoxes(r, triangleColor(r));
    expect(triangles).toBe(6);
    for (const m of boxes) {
      expect(m[0]).toBeGreaterThanOrEqual(Math.round(l.dataLeft * 2));
      expect(m[2]).toBeLessThanOrEqual(Math.round(l.dataRight * 2));
    }
    // The hit: a few px each side of the middle boundary, top to bottom of the row, is the insertion.
    const a = t.ctx.layout.axis(S);
    const x = a ? a.x(0, a.shownOf(0, 20)) : 0;
    const y = t.yOf(S, "t-of");
    const hit = Hit.create();
    const at = (px: number, py: number) => {
      t.ctx.layout.hitTest(S, px, py, hit);
      return hit.kind === HoverKind.item ? hit.item : hit.kind === HoverKind.seq ? -1 : -2;
    };
    const rowTop = l.scrolledTop + l.entryTop[l.entryTrack.indexOf(t.trackIndex("t-of"))];
    for (const dx of [-3, 0, 3]) for (const py of [rowTop + 1, y, rowTop + 18]) expect(at(x + dx, py)).toBe(1);
    // Beside it, the letters' residues.
    expect(at(x - 7, y)).toBe(-1);
    expect(at(x + 7, y)).toBe(-1);
  });
});
