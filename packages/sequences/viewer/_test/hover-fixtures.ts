/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Fixtures for the hover emphasis (HV1; seqspec:interaction#hover), browser-safe (the builder only): one chain with a
// track per emphasized type (ranges with a two-fragment item and an item over the same positions, points, variants
// stacked at a position, secondary structure, values as bars and as a line, a matrix), and readers of what a recorded
// pass drew as boxes: rects (outlines and fills), rings (a stroked arc's outer box) and glyphs.

import { createSeqViewSpec } from "@molstar/seq-core/builder";
import { Op, type RecordingBackend } from "../src/testing/index.js";

const PROTEIN = "MSEEKLKELAEKAIRLAKEGNPEEALKLLEEALRLAKESGDPELIARVLELLAEVTKQSVHAYRLDDGKWFLQPGGRVSLEMAPTDLQAAH";

/**
 * The hover document: `n` residues, a track per type, each by ref (`t.<type>`); `lineColor`: the line's color theme
 * (none: one color).
 */
export function hoverDoc(n = 100, o: { style?: Record<string, unknown>; lineColor?: unknown } = {}) {
  const b = createSeqViewSpec({ title: "Hover" });
  const a = b.data.sequence({
    ref: "a",
    label: "A",
    value: Array.from({ length: n }, (_, i) => PROTEIN[(i * 7) % PROTEIN.length]).join(""),
  });
  const ann = (ref: string, geometry: string, value: unknown, more: object = {}) =>
    b.data.annotation({ ref, geometry, target: a, value, ...more } as never);
  // g1 in two fragments, g2 over g1's positions (another lane), g3 alone.
  const rg = ann("rg", "ranges", {
    id: ["g1", "g1", "g2", "g3"],
    start: [10, 30, 14, 60],
    stop: [20, 45, 40, 80],
    category: ["x", "x", "y", "z"],
  });
  const pt = ann("pt", "ranges", { start: [12, 50, 70], stop: [13, 51, 71], category: ["p", "q", "p"] });
  const vr = ann("vr", "ranges", { start: [12, 12, 50], stop: [13, 13, 51], alt: ["A", "C", "D"] });
  const ss = ann("ss", "ranges", { start: [5, 25, 55], stop: [15, 35, 65], category: ["helix", "strand", "turn"] });
  const index = Array.from({ length: n }, (_, i) => i);
  const vv = ann("vv", "values", { index, value: index.map((i) => 0.2 + ((i * 31) % 80) / 100) });
  const mx = ann(
    "mx",
    "matrix",
    { index, value: index.map((i) => [(i % 7) / 7, (i % 11) / 11, (i % 3) / 3]) },
    { states: ["x", "y", "z"] },
  );
  const v = b.view({ ref: "v", label: "Hover", height: "box", ...(o.style ? { style: o.style } : {}) } as never);
  v.block({ sequence: a });
  v.track({ ref: "t.letters", label: "Letters" }).letters();
  v.track({ ref: "t.ranges", label: "Ranges" }).ranges({ data: rg });
  v.track({ ref: "t.points", label: "Points" }).points({ data: pt });
  v.track({ ref: "t.variants", label: "Variants" }).variants({ data: vr });
  v.track({ ref: "t.secondary", label: "Secondary" }).secondary({ data: ss });
  v.track({ ref: "t.values", label: "Values" }).values({ data: vv });
  const lineColor = o.lineColor ? { color: o.lineColor } : {};
  v.track({ ref: "t.line", label: "Line" }).values({ data: vv, encoding: "line", ...lineColor } as never);
  v.track({ ref: "t.matrix", label: "Matrix" }).matrix({ data: mx });
  return b.getDocument();
}

/**
 * A drawn box: x0, y0, x1, y1 (device px), its color, and its kind (0 rect fill, 1 rect outline, 2 ring, 3 glyph, 4 a
 * path's bounding box).
 */
export type Box = [number, number, number, number, number, number];
export const BoxKind = { fill: 0, outline: 1, ring: 2, glyph: 3, shape: 4 } as const;

/** What pass `pass` drew, as boxes, in order. */
export function passBoxes(r: RecordingBackend, pass: number): Box[] {
  const out: Box[] = [];
  const box = [Number.POSITIVE_INFINITY, Number.POSITIVE_INFINITY, Number.NEGATIVE_INFINITY, Number.NEGATIVE_INFINITY];
  let arc = false;
  let at = -1;
  for (let k = 0; k < r.opCount; k++) {
    const code = r.ops[4 * k];
    if (code === Op.beginPass) at = r.args[r.ops[4 * k + 3]];
    if (at !== pass) continue;
    const i0 = r.ops[4 * k + 1];
    const i1 = i0 + r.ops[4 * k + 2];
    if (code === Op.rects) {
      const b = r.rectItems;
      for (let i = i0; i < i1; i++)
        out.push([b.x[i], b.y[i], b.x[i] + b.w[i], b.y[i] + b.h[i], b.color[i], b.stroke[i] > 0 ? 1 : 0]);
    } else if (code === Op.glyphs) {
      const g = r.glyphItems;
      for (let i = i0; i < i1; i++)
        out.push([
          g.x[i] - g.radius[i],
          g.y[i] - g.radius[i],
          g.x[i] + g.radius[i],
          g.y[i] + g.radius[i],
          g.color[i],
          3,
        ]);
    } else if (code === Op.paths) {
      const p = r.pathItems;
      for (let q = i0; q < i1; q++) {
        const end = q + 1 < p.n ? p.start[q + 1] : p.nv;
        let c = p.cstart[q];
        for (let v = p.start[q]; v < end; v++) {
          const verb = p.verbs[v];
          if (verb === 4) {
            // An arc: its stroke's outer box (the width is in half device px).
            const R = p.coords[c + 2] + (p.mode[q] === 1 ? p.width[q] / 4 : 0);
            out.push([p.coords[c] - R, p.coords[c + 1] - R, p.coords[c] + R, p.coords[c + 1] + R, p.color[q], 2]);
            arc = true;
          } else
            for (let j = 0; j + 1 < [2, 2, 4, 6, 6, 5, 0][verb]; j += 2) {
              box[0] = Math.min(box[0], p.coords[c + j]);
              box[1] = Math.min(box[1], p.coords[c + j + 1]);
              box[2] = Math.max(box[2], p.coords[c + j]);
              box[3] = Math.max(box[3], p.coords[c + j + 1]);
            }
          c += [2, 2, 4, 6, 6, 5, 0][verb];
        }
        if (!arc && box[0] <= box[2]) out.push([box[0], box[1], box[2], box[3], p.color[q], 4]);
        box.fill(Number.POSITIVE_INFINITY, 0, 2).fill(Number.NEGATIVE_INFINITY, 2, 4);
        arc = false;
      }
    }
  }
  return out;
}
/** Whether box `a` lies within box `b` (a tolerance of `eps` device px). */
export const within = (a: readonly number[], b: readonly number[], eps = 0) =>
  a[0] >= b[0] - eps && a[1] >= b[1] - eps && a[2] <= b[2] + eps && a[3] <= b[3] + eps;
