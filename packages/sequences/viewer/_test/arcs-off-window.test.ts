/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A `pairs` arc whose partner lies beyond the window (seqspec:format#layer-types `pairs`): case 7's contacts join
// barnase (block 0) and barstar (block 1), side by side, so zoomed into one block every partner is off-window in the
// other. Each pair is its own subpath within the data area and the band, on the hit's curve: the parabola over its
// virtual ends (a partner taken at most a margin past the edge), never flat (the first bug: the widest span in view
// set the scale, so a partner a few windows off made every arc a line along the ends), never out of the band (the
// second: an apex past the band, a steep line out of it). A pair spanning the window is the widest arc over it;
// stacked, a partner on another lane is a tick.

import type { SeqViewSpecDoc } from "@molstar/seq-core/format";
import { describe, expect, it } from "vitest";
import { case7 } from "../../core/_test/format/cases.ts";
import { pushArc } from "../src/draw/arcs.ts";
import { PathBatch, PathMode, PathVerb } from "../src/draw/backend.ts";
import type { DrawView } from "../src/draw/view.ts";
import type { PairsState } from "../src/prepare.ts";
import { ARC_MARGIN, ARC_PAST, Arcs } from "../src/renderer-helpers.ts";
import { drawn, layerOf, record } from "./draw-fixtures.ts";

/** Case 7's "Complex" view, its section's blocks side by side (or stacked), `block`'s [start, stop) in view. */
async function zoomed(block: number, start: number, stop: number, doc: SeqViewSpecDoc = case7()) {
  const t = await drawn(doc);
  t.ctx.views.setCurrent("complex-one");
  await t.ctx.graph.update();
  for (const ref of t.ctx.viewports.sections.value) t.ctx.viewports.get(ref)?.setSize(1000, 400, 1);
  await t.ctx.graph.update();
  const { k, section } = layerOf(t.ctx, "Contacts");
  t.ctx.viewports.get(section)?.zoomTo(block, start, stop);
  await t.ctx.graph.update();
  const s = t.ctx.graph.assembled.value?.view?.prepared[k].state as PairsState;
  const { rec, frame } = record(t.ctx, "Contacts");
  return { s, rec, v: frame.v };
}

function stacked(): SeqViewSpecDoc {
  const doc = case7();
  type Node = { kind: string; params?: Record<string, unknown>; children?: Node[] };
  const walk = (n: Node): void => {
    if (n.kind === "section") n.params = { ...n.params, arrangement: "stacked" };
    n.children?.forEach(walk);
  };
  walk(doc.root as Node);
  return doc;
}

/** The stroked subpaths of a batch: per subpath its verbs and coordinates (CSS px at dpr 1). */
function subpaths(b: PathBatch) {
  return Array.from({ length: b.n }, (_, i) => {
    const v1 = i + 1 < b.n ? b.start[i + 1] : b.nv;
    const c1 = i + 1 < b.n ? b.cstart[i + 1] : b.nc;
    return {
      mode: b.mode[i],
      verbs: Array.from(b.verbs.subarray(b.start[i], v1)),
      coords: Array.from(b.coords.subarray(b.cstart[i], c1)),
    };
  });
}

/** Points on a subpath's curve: its move, 11 along each quad, and each line's end (the polyline's vertices). */
function samples(verbs: number[], c: number[]): [number, number][] {
  expect(verbs[0]).toBe(PathVerb.move);
  const out: [number, number][] = [[c[0], c[1]]];
  let [x, y] = [c[0], c[1]];
  for (let i = 1, o = 2; i < verbs.length; i++) {
    if (verbs[i] === PathVerb.quad) {
      for (let j = 1; j <= 10; j++) {
        const t = j / 10;
        const [a, b, q] = [(1 - t) * (1 - t), 2 * t * (1 - t), t * t];
        out.push([a * x + b * c[o] + q * c[o + 2], a * y + b * c[o + 1] + q * c[o + 3]]);
      }
      [x, y] = [c[o + 2], c[o + 3]];
      o += 4;
    } else {
      expect(verbs[i]).toBe(PathVerb.line);
      [x, y] = [c[o], c[o + 1]];
      out.push([x, y]);
      o += 2;
    }
  }
  return out;
}

/** The drawn height at x of the pair in `s.xs` (after `Arcs.apex`): `Arcs.heightAt`, the hit's. */
function heightAt(s: PairsState, x: number): number {
  s.frame[10] = x;
  Arcs.heightAt(s);
  return s.frame[11];
}

/**
 * Every pair with an end in view and its partner beyond it, pushed alone: its own subpath, starting at its end in view
 * on the ends' line, within the data area's x and the band, on `Arcs.heightAt`'s curve (which the hit shares): the
 * parabola over its virtual ends, the partner taken at most `ARC_MARGIN` of the window past the edge. Never flat toward
 * the edge: a partner beyond the margin reaches there at least `(2/3) · band · d / (d + margin)` (d its end's distance
 * from the edge; the scale at most the window and both margins). Its count, and how many partners were beyond the
 * margin.
 */
function checkOffWindow(s: PairsState, v: DrawView, rec: PathBatch, beyond: "left" | "right") {
  const f = s.frame;
  const W = v.right - v.left;
  const M = ARC_MARGIN * W;
  const band = f[2];
  let checked = 0;
  let far = 0;
  let drawnArcs = 0;
  const h = (y: number) => (y - f[0]) * f[1];
  for (let k = 0; k < s.count; k++) {
    const arcs = new PathBatch();
    const wedges = new PathBatch();
    const any = pushArc(v, s, k, 2, 0, arcs, wedges);
    if (any) drawnArcs++;
    const [xa, xb] = [s.xs[0], s.xs[1]];
    const inA = xa >= v.left && xa <= v.right;
    const inB = xb >= v.left && xb <= v.right;
    if (inA === inB) continue;
    const xin = inA ? xa : xb;
    const xout = inA ? xb : xa;
    expect(beyond === "left" ? xout < v.left : xout > v.right).toBe(true);
    // One subpath, a quadratic, its coordinates near the view at any zoom, within the band, on the hit's curve.
    expect(any).toBe(true);
    const [p] = subpaths(arcs);
    expect(arcs.n).toBe(1);
    expect(p.mode).toBe(PathMode.stroke);
    expect(p.verbs).toEqual([PathVerb.move, PathVerb.quad]);
    const pts = samples(p.verbs, p.coords);
    for (const [x, y] of pts) {
      expect(x).toBeGreaterThanOrEqual(v.left - ARC_PAST - 0.01);
      expect(x).toBeLessThanOrEqual(v.right + ARC_PAST + 0.01);
      expect(h(y)).toBeGreaterThanOrEqual(-0.01);
      expect(h(y)).toBeLessThanOrEqual(band + 0.01);
      expect(h(y)).toBeCloseTo(heightAt(s, x), 6);
    }
    // It starts at its end in view, on the ends' line.
    const first = Math.abs(p.coords[0] - xin) < 1e-6 ? 0 : p.coords.length - 2;
    expect(p.coords[first]).toBeCloseTo(xin, 6);
    expect(p.coords[first + 1]).toBeCloseTo(f[0], 6);
    // On the parabola over the end in view and its partner taken in at the margin.
    const vout = beyond === "left" ? Math.max(xout, v.left - M) : Math.min(xout, v.right + M);
    const mid = pts[5];
    const c = h(mid[1]) / ((mid[0] - xin) * (vout - mid[0]));
    expect(c).toBeGreaterThan(0);
    for (const [x, y] of pts) expect(h(y)).toBeCloseTo(c * (x - xin) * (vout - x), 6);
    // Not flat toward the edge.
    const edge = beyond === "left" ? v.left : v.right;
    const d = Math.abs(edge - xin);
    if (Math.abs(xout - edge) >= M) {
      far++;
      expect(heightAt(s, edge)).toBeGreaterThanOrEqual(((2 / 3) * band * d) / (d + M) - 0.01);
    }
    checked++;
  }
  // The draw pushed the same arcs, one subpath each.
  expect(subpaths(rec).filter((p) => p.mode === PathMode.stroke).length).toBe(drawnArcs);
  return { checked, far };
}

describe("pairs: an arc whose partner is beyond the window", () => {
  for (const [start, stop] of [
    [34, 86],
    [36, 60],
  ])
    it(`barstar ${start}–${stop}: each partner left in barnase, the true piece rising to the left edge`, async () => {
      const { s, rec, v } = await zoomed(1, start, stop);
      const { checked, far } = checkOffWindow(s, v, rec.pathItems, "left");
      expect(checked).toBeGreaterThan(5);
      expect(far).toBeGreaterThan(0);
    });

  for (const [start, stop] of [
    [0, 40],
    [20, 45],
  ])
    it(`barnase ${start}–${stop}: each partner right in barstar, the true piece rising to the right edge`, async () => {
      const { s, rec, v } = await zoomed(0, start, stop);
      const { checked, far } = checkOffWindow(s, v, rec.pathItems, "right");
      expect(checked).toBeGreaterThan(2);
      expect(far).toBeGreaterThan(0);
    });

  it("a pair spanning the window is the widest arc over it, its apex the band's top at the window's middle", async () => {
    // Barstar 36–60: the contacts at barstar 73 and 76 (off right) to barnase (off left).
    const { s, v } = await zoomed(1, 36, 60);
    let spanning = 0;
    for (let k = 0; k < s.count; k++) {
      const arcs = new PathBatch();
      const any = pushArc(v, s, k, 2, 0, arcs, new PathBatch());
      if (!(Math.min(s.xs[0], s.xs[1]) < v.left && Math.max(s.xs[0], s.xs[1]) > v.right)) continue;
      spanning++;
      expect(any).toBe(true);
      const [p] = subpaths(arcs);
      expect(p.verbs).toEqual([PathVerb.move, PathVerb.quad]);
      for (const [x, y] of samples(p.verbs, p.coords)) {
        const up = (y - s.frame[0]) * s.frame[1];
        expect(up).toBeGreaterThan(0);
        expect(up).toBeLessThanOrEqual(s.frame[2] + 0.01);
        expect(up).toBeCloseTo(heightAt(s, x), 6);
      }
      expect(heightAt(s, (v.left + v.right) / 2)).toBeCloseTo(s.frame[2], 6);
    }
    expect(spanning).toBeGreaterThan(0);
  });

  it("stacked, a partner on the other block's lane is a tick at its end, its own subpath", async () => {
    for (const [start, stop] of [
      [34, 86],
      [36, 60],
    ]) {
      const { rec } = await zoomed(1, start, stop, stacked());
      const ticks = subpaths(rec.pathItems).filter((p) => p.mode === PathMode.stroke);
      expect(ticks.length).toBeGreaterThan(5);
      for (const p of ticks) {
        expect(p.verbs).toEqual([PathVerb.move, PathVerb.line]);
        expect(p.coords[2]).toBe(p.coords[0]);
      }
    }
  });
});
