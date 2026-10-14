/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// CG2 (seqspec:viewer#the-render-backend): items' labels outside the rings with leaders (placed by angle, a line apart
// per side, dropped by priority), arrowheads by a direction field (on rings and linearly), and pairs as chords through
// the hub.

import { createSeqViewSpec } from "@molstar/seq-core/builder";
import { describe, expect, it } from "vitest";
import { PathVerb, SkinToken } from "../src/draw/index.js";
import { RingLabels } from "../src/draw/ring-labels.js";
import { DefaultPresets } from "../src/presets.js";
import { Op, type RecordingBackend } from "../src/testing/index.js";
import { composed, recording } from "./composer-fixtures.js";
import { opened } from "./layout-fixtures.js";

const registry = [{ presets: DefaultPresets.filter((p) => p.name === "circular") }];
const N = 600;

/** A cyclic 600-nt sequence: three stranded features (one across the origin), two cut sites, a pair. */
function doc(style: "circular" | null) {
  const b = createSeqViewSpec();
  const s = b.data.sequence({ ref: "s", value: "ACGTACGTTA".repeat(N / 10), alphabet: "dna", cyclic: true });
  const f = b.data.annotation({
    ref: "f",
    geometry: "ranges",
    target: s,
    value: {
      start: [100, 300, 550],
      stop: [200, 400, 50],
      label: ["alpha", "beta", "gamma"],
      strand: ["+", "-", "+"],
    },
  } as never);
  const cuts = b.data.annotation({
    ref: "cuts",
    geometry: "ranges",
    target: s,
    value: { boundary: [120, 450], label: ["CutA", "CutB"] },
  } as never);
  const pairs = b.data.annotation({
    ref: "pp",
    geometry: "pairs",
    target: s,
    value: { start: [50], b_target: ["s"], b_start: [350] },
  } as never);
  const v = b.view({ ref: "v", label: "V", ...(style ? { style } : {}) }).block({ ref: "v.b", sequence: s });
  v.track({ label: "" }).ruler();
  v.track({ label: "F", height: 2 }).ranges({ data: f, labels: true, direction_field: "strand" } as never);
  v.track({ label: "C" }).ranges({ data: cuts, labels: true } as never);
  v.track({ label: "P" }).pairs({ data: pairs, display: "arcs" } as never);
  return b.getDocument();
}

/** The recording's texts with their batch's alignment, and its paths' verbs and points (device px). */
function texts(r: RecordingBackend) {
  const out: { text: string; x: number; y: number; color: number }[] = [];
  for (let k = 0; k < r.opCount; k++) {
    if (r.ops[4 * k] !== Op.text) continue;
    for (let i = r.ops[4 * k + 1]; i < r.ops[4 * k + 1] + r.ops[4 * k + 2]; i++)
      out.push({ text: r.textItems.text[i], x: r.textItems.x[i], y: r.textItems.y[i], color: r.textItems.color[i] });
  }
  return out;
}
function paths(r: RecordingBackend) {
  const b = r.pathItems;
  return Array.from({ length: b.n }, (_, q) => {
    const v1 = q + 1 < b.n ? b.start[q + 1] : b.nv;
    const c1 = q + 1 < b.n ? b.cstart[q + 1] : b.nc;
    return {
      verbs: Array.from(b.verbs.subarray(b.start[q], v1)),
      coords: Array.from(b.coords.subarray(b.cstart[q], c1)),
      color: b.color[q],
      mode: b.mode[q],
    };
  });
}

describe("labels on rings", () => {
  it("are set outside the ring past the ruler's numbers, each with a leader from its item, once per item", async () => {
    const t = await composed(doc("circular"), { width: 700, height: 500, registry });
    const l = t.ctx.layout.get(t.section)?.current;
    if (!l) throw new Error("no layout");
    const r = t.frame(recording());
    const labels = texts(r).filter((x) => /^[a-zA-Z]/.test(x.text));
    // gamma across the origin is two runs, one item: labeled once; no label inside the features' ring.
    expect(labels.map((x) => x.text).sort()).toEqual(["CutA", "CutB", "alpha", "beta", "gamma"]);
    for (const x of labels) expect(Math.hypot(x.x - l.ringCx, x.y - l.ringCy)).toBeGreaterThan(l.ringR + 10);
    const leaders = paths(r).filter((p) => p.color === SkinToken.inkQuiet && p.verbs.length === 4);
    expect(leaders).toHaveLength(5);
    // Each leader starts inside the ring (on its item's ring), leaves it at the rings' edge on its item's angle, and ends
    // beside its label; every point is inside the section.
    for (const p of leaders) {
      const at = (j: number) => [p.coords[2 * j] - l.ringCx, p.coords[2 * j + 1] - l.ringCy];
      const [x0, y0] = at(0);
      const [x1, y1] = at(1);
      expect(Math.hypot(x0, y0)).toBeLessThan(l.ringR);
      expect(Math.abs(Math.hypot(x1, y1) - (l.ringR + 2))).toBeLessThan(1e-6);
      expect(Math.abs(Math.atan2(y1, x1) - Math.atan2(y0, x0))).toBeLessThan(1e-9);
      for (let j = 0; j < 4; j++) {
        expect(p.coords[2 * j]).toBeGreaterThanOrEqual(0);
        expect(p.coords[2 * j]).toBeLessThanOrEqual(700);
        expect(p.coords[2 * j + 1]).toBeGreaterThanOrEqual(0);
        expect(p.coords[2 * j + 1]).toBeLessThanOrEqual(l.height);
      }
    }
    for (const x of labels) expect(x.x >= 0 && x.x <= 700 && x.y >= 0 && x.y <= l.height).toBe(true);
    // Room was kept beside the ring for them: the ring is narrower than the section by twice that room at least.
    expect(2 * l.ringR).toBeLessThanOrEqual(700 - 2 * 96);
    expect(t.ctx.issues.current.value.filter((i) => i.code === "seqv.render.labels")).toEqual([]);
    t.ctx.dispose();
  });

  it("follow a rotated origin: an item whose run crosses it on screen is labeled at its middle", async () => {
    const t = await composed(doc("circular"), { width: 700, height: 500, registry });
    await t.ctx.applyPatch({
      format: "seq_view_spec_patch",
      ops: [{ op: "set", ref: "v.b", params: { origin: 350 } }],
    });
    const l = t.ctx.layout.get(t.section)?.current;
    if (!l) throw new Error("no layout");
    // beta, [300, 400), now from 550 to 50 of the ring's 600: its middle (349) just before the top.
    const beta = texts(t.frame(recording())).find((x) => x.text === "beta");
    expect(beta?.y).toBeLessThan(l.ringCy - l.ringR / 2);
    t.ctx.dispose();
  });

  it("are spread a line apart per side, and those a side cannot hold are dropped by priority", async () => {
    const t = await opened(doc(null));
    const items = Array.from({ length: 12 }, (_, i) => ({ angle: 0.1, r0: 100, text: `L${i}`, priority: i }));
    items.push({ angle: Math.PI, r0: 100, text: "left", priority: 0 });
    const frame = { cx: 200, cy: 200, edge: 104, r: 120, left: 0, right: 400, top: 140, bottom: 270, line: 13 };
    const placed = new RingLabels().place(items, frame, t.ctx.text);
    // 130 px hold 10 lines: L0 … L9 on the right, the left side's one.
    expect(placed.text.slice().sort()).toEqual(["L0", "L1", "L2", "L3", "L4", "L5", "L6", "L7", "L8", "L9", "left"]);
    expect(placed.hidden).toBe(2);
    const ys = Array.from({ length: placed.n }, (_, k) => [placed.side[k], placed.pts[10 * k + 9]] as const)
      .filter(([side]) => side > 0)
      .map(([, y]) => y)
      .sort((a, b) => a - b);
    for (let k = 1; k < ys.length; k++) expect(ys[k] - ys[k - 1]).toBeGreaterThanOrEqual(13 - 1e-9);
    expect(ys[0]).toBeGreaterThanOrEqual(140);
    expect(ys.at(-1)).toBeLessThanOrEqual(270);
    t.ctx.dispose();
  });
});

describe("arrowheads by a direction field", () => {
  it("a forward item's head at its right end, a reverse one's at its left, linearly", async () => {
    const t = await composed(doc(null), { width: 1000, height: 300 });
    t.vp?.setView(1.5, 0, "api");
    const r = t.frame(recording());
    const heads = paths(r).filter((p) => p.mode === 0 && p.verbs.length === 4 && p.verbs[3] === PathVerb.close);
    // Each head's tip (its second point) and base: alpha (+) and gamma's first run point right, beta (−) left.
    const dirs = heads.map((p) => Math.sign(p.coords[2] - p.coords[0])).sort();
    expect(dirs).toEqual([-1, 1, 1]);
    t.ctx.dispose();
  });
});

describe("pairs as chords on rings", () => {
  it("a pair's curve bends toward the hub, through it for opposite ends", async () => {
    const t = await composed(doc("circular"), { width: 700, height: 500, registry });
    const l = t.ctx.layout.get(t.section)?.current;
    if (!l) throw new Error("no layout");
    const r = t.frame(recording());
    const chords = paths(r).filter((p) => p.verbs.includes(PathVerb.quad));
    expect(chords).toHaveLength(1);
    const [x0, y0, cx, cy, x1, y1] = chords[0].coords;
    // Ends 300 of 600 apart: opposite, so the control point is the center.
    expect(Math.hypot(cx - l.ringCx, cy - l.ringCy)).toBeLessThan(1);
    // Each end on the pairs ring's outer edge (its ends' line, 3 px in: `display: "arcs"`), against the ring outside
    // it, at its position's angle: the middle of positions 50 and 350's cells, clockwise from the top.
    const av = t.ctx.layout.axis(t.section);
    let e = 0;
    while (e < l.entryCount && l.entryTrack[e] !== l.entryTrack[l.entryCount - 1]) e++;
    const radius = l.ringR - (l.entryTop[e] + 3) * l.ringScale;
    const angle = (q: number) =>
      -Math.PI / 2 + (2 * Math.PI * ((av?.x(0, q as never) ?? 0) + 0.5 * (av?.cell ?? 0) - l.dataLeft)) / l.ringW;
    for (const [x, y, q] of [
      [x0, y0, 50],
      [x1, y1, 350],
    ]) {
      expect(Math.abs(Math.hypot(x - l.ringCx, y - l.ringCy) - radius)).toBeLessThan(0.5);
      const a = Math.atan2(y - l.ringCy, x - l.ringCx);
      // Within a quarter of a cell's angle (the cell's middle as the axis view rounds it).
      expect(Math.abs(Math.atan2(Math.sin(a - angle(q)), Math.cos(a - angle(q))))).toBeLessThan(
        (0.25 * 2 * Math.PI) / N,
      );
    }
    t.ctx.dispose();
  });
});
