/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Rings (CG1a; seqspec:viewer#the-render-backend): the `circular` preset's layout (rings from the outside in, the
// section its diameter), every draw mapped onto the circle (rects as annular sectors at the angles of their positions),
// a feature across the origin as one arc on the lane it shares (packing by arc), polar hits, reductions per angular
// column, and the shortest arc as a focus hull. CG1b: marks as sectors across the rings, the hover as a radial line and
// the hovered item's outline, the ruler's minor ticks, the keyboard cursor going round, and a type with no ring form
// drawn as bars.

import { createSeqViewSpec } from "@molstar/seq-core/builder";
import { Extents } from "@molstar/seq-core/model";
import { describe, expect, it } from "vitest";
import { Dirty, Pass, SkinToken } from "../src/draw/index.ts";
import { Hit, HoverKind, KernelRenderers, type LayerRenderer } from "../src/index.ts";
import type { Features } from "../src/prepare.ts";
import { DefaultPresets } from "../src/presets.ts";
import { Op, type RecordingBackend } from "../src/testing/index.ts";
import { composed, recording } from "./composer-fixtures.ts";
import { driven } from "./input-fixtures.ts";
import { layer } from "./layout-fixtures.ts";

const registry = [{ presets: DefaultPresets.filter((p) => p.name === "circular") }];
const TAU = 2 * Math.PI;

/** A pass's sectors: [r0, r1, a0, a1, color, alpha]; or its paths' points and width: [[x, y]…, color, width]. */
function passSectors(r: RecordingBackend, pass: number): number[][] {
  const b = r.sectorItems;
  const out: number[][] = [];
  let at = -1;
  for (let k = 0; k < r.opCount; k++) {
    const code = r.ops[4 * k];
    if (code === Op.beginPass) at = r.args[r.ops[4 * k + 3]];
    if (code !== Op.sectors || at !== pass) continue;
    for (let i = r.ops[4 * k + 1]; i < r.ops[4 * k + 1] + r.ops[4 * k + 2]; i++)
      out.push([b.r0[i], b.r1[i], b.a0[i], b.a1[i], b.color[i], b.alpha[i]]);
  }
  return out;
}
function passPaths(r: RecordingBackend, pass: number): { points: number[][]; color: number; width: number }[] {
  const b = r.pathItems;
  const out: { points: number[][]; color: number; width: number }[] = [];
  let at = -1;
  for (let k = 0; k < r.opCount; k++) {
    const code = r.ops[4 * k];
    if (code === Op.beginPass) at = r.args[r.ops[4 * k + 3]];
    if (code !== Op.paths || at !== pass) continue;
    for (let q = r.ops[4 * k + 1]; q < r.ops[4 * k + 1] + r.ops[4 * k + 2]; q++) {
      const c0 = b.cstart[q];
      const c1 = q + 1 < b.n ? b.cstart[q + 1] : b.nc;
      const points: number[][] = [];
      for (let c = c0; c + 1 < c1; c += 2) points.push([b.coords[c], b.coords[c + 1]]);
      out.push({ points, color: b.color[q], width: b.width[q] });
    }
  }
  return out;
}
/** Angles as turns from the top, clockwise, in [0, 1). */
const turn = (a: number) => ((((a + Math.PI / 2) / TAU) % 1) + 1) % 1;

/** A cyclic sequence of `n` (240 by default) on rings: a ruler, its letters, and features [200, 30) and [40, 112). */
function plasmid(n = 240, style: string | null = "circular") {
  const b = createSeqViewSpec();
  const s = b.data.sequence({ ref: "s", value: "ACGTACGTTA".repeat(n / 10), alphabet: "dna", cyclic: true });
  const f = b.data.annotation({
    ref: "f",
    geometry: "ranges",
    target: s,
    value: { start: [200, 40], stop: [30, 112] },
  });
  const v = b.view({ ref: "v", label: "V", ...(style ? { style } : {}) }).block({ sequence: s });
  v.track({ label: "" }).ruler();
  v.track({ label: "S" }).letters();
  v.track({ label: "F" }).ranges({ data: f });
  return b.getDocument();
}

async function rings(n = 240) {
  const t = await composed(plasmid(n), { width: 500, height: 500, registry });
  const l = t.ctx.layout.get(t.section)?.current;
  const av = t.ctx.layout.axis(t.section);
  if (!l || !av) throw new Error("no layout");
  /** The angle of shown k (block 0), and the radius at surface y. */
  const angle = (k: number) => -Math.PI / 2 + (TAU * (av.x(0, k as never) - l.dataLeft)) / l.ringW;
  const radius = (y: number) => l.ringR - (y - l.scrolledTop) * l.ringScale;
  const entry = (label: string) => {
    for (let e = 0; e < l.entryCount; e++)
      if (l.entryLayer[e] >= 0 && layer(t.ctx, label).k === l.entryLayer[e]) return e;
    return -1;
  };
  return { t, l, av, angle, radius, entry };
}

describe("rings", () => {
  it("lays the entries out as rings from the outside in, the section its diameter and label ring", async () => {
    const { l } = await rings();
    expect(l.geometry).toBe(2);
    expect([l.ringCx, l.entryCount]).toEqual([250, 3]);
    // Each ring below the last, `ringGap` (2) apart, all of them outside the hub (0.35 of the radius).
    for (let e = 1; e < l.entryCount; e++) expect(l.entryTop[e]).toBe(l.entryTop[e - 1] + l.entryHeight[e - 1] + 2);
    const depth = (l.entryTop[l.entryCount - 1] + l.entryHeight[l.entryCount - 1]) * l.ringScale;
    expect(l.ringR - depth).toBeGreaterThanOrEqual(0.35 * l.ringR - 1e-9);
    // As tall as its width (500, within [320, 560]); the ring inside it with room for its labels.
    expect([l.ringW, l.contentHeight, l.ringCy]).toEqual([TAU * l.ringR, 500, 250]);
    expect(l.ringR).toBeLessThanOrEqual(250 - 8);
  });

  it("fills its box: a section with no height is as tall as wide within [320, 560], the ring the box less its labels", async () => {
    const t = await composed(plasmid(), { width: 1800, height: 200, registry });
    const l0 = t.ctx.layout.get(t.section)?.current;
    // The section asks for 560 px (its content height); sized so, the ring fills it less the label ring.
    expect(t.ctx.layout.get(t.section) && l0?.contentHeight).toBe(560);
    t.vp?.setSize(1800, 560, 1);
    await t.ctx.graph.update();
    const l = t.ctx.layout.get(t.section)?.current;
    if (!l) throw new Error("no layout");
    expect([l.ringCx, l.ringCy]).toEqual([900, 280]);
    // The margin is the widest number ("240") past the ticks less the ruler's room above them: the ring grows into it.
    expect(l.ringR).toBeGreaterThan(260);
    expect(l.ringR).toBeLessThan(280 - 8 + 1e-9);
    // The tracks fill the ring down to the hub (0.35 of the radius), each as wide as its linear height scaled alike.
    const last = l.entryCount - 1;
    expect((l.entryTop[last] + l.entryHeight[last]) * l.ringScale).toBeCloseTo(0.65 * l.ringR, 6);
    for (let e = 0; e < l.entryCount; e++) expect(l.entryHeight[e] * l.ringScale).toBeGreaterThan(3);
    // A narrow section: as tall as its width, at least 320; a `box` section's own height wins.
    const narrow = await composed(plasmid(), { width: 400, height: 0, registry });
    expect(narrow.ctx.layout.get(narrow.section)?.current.contentHeight).toBe(400);
  });

  it("draws a feature as an annular sector at its positions' angles, one across the origin as its arc's two runs", async () => {
    const { t, l, angle, radius, entry } = await rings();
    const r = t.frame(recording());
    const s = r.sectorItems;
    const e = entry("F");
    const [r1, r0] = [radius(l.scrolledTop + l.entryTop[e]), radius(l.scrolledTop + l.entryTop[e] + l.entryHeight[e])];
    // Sectors within the features' ring, by their angles (device px = CSS px at dpr 1; edges rounded to a px).
    const near = (a: number, b: number) => Math.abs(a - b) < TAU / l.ringW + 1e-9;
    const on = Array.from({ length: s.n }, (_, i) => i).filter((i) => s.r0[i] >= r0 - 1 && s.r1[i] <= r1 + 1);
    const has = (k0: number, k1: number) => on.some((i) => near(s.a0[i], angle(k0)) && near(s.a1[i], angle(k1)));
    expect([has(40, 112), has(200, 240), has(0, 30)]).toEqual([true, true, true]);
    // No rects: everything on a ring is a sector, and every sector is centered on the ring's center.
    expect(r.rectItems.n).toBe(0);
    for (let i = 0; i < s.n; i++) expect([s.cx[i], s.cy[i]]).toEqual([l.ringCx, l.ringCy]);
  });

  it("sets the ruler's numbers just past its ticks, aligned away from the center", async () => {
    const { t, l } = await rings();
    const r = t.frame(recording());
    const s = r.sectorItems;
    let ticks = 0;
    for (let i = 0; i < s.n; i++) ticks = Math.max(ticks, s.r1[i]);
    const tx = r.textItems;
    expect(tx.n).toBeGreaterThan(10);
    for (let i = 0; i < tx.n; i++) expect(Math.hypot(tx.x[i] - l.ringCx, tx.y[i] - l.ringCy)).toBeCloseTo(ticks + 5, 6);
    expect(ticks).toBeLessThan(l.ringR);
  });

  it("packs by arc: the feature across the origin shares its lane with one outside its arc", async () => {
    const { t } = await rings();
    const f = layer(t.ctx, "F").state as Features;
    expect([f.lanes, Array.from(f.lane.subarray(0, f.count))]).toEqual([1, [0, 0, 0]]);
    expect(Array.from(f.joinFrom.subarray(0, f.count))).toEqual([-1, -1, -1]);
  });

  it("the story's plasmid: the replication origin, lacZα and AmpR on one lane, the MCS on the next", async () => {
    const b = createSeqViewSpec();
    const p = b.data.sequence({
      ref: "p",
      alphabet: "dna",
      cyclic: true,
      source: { host_file: "pDemo.fasta", format: "fasta" },
    });
    const f = b.data.annotation({
      ref: "f",
      geometry: "ranges",
      target: p,
      source: { host_file: "pDemo.tsv", format: "tsv" },
    });
    b.view({ ref: "v", label: "V", style: "circular" })
      .block({ sequence: p })
      .track({ label: "F" })
      .ranges({ data: f, labels: true });
    const t = await composed(b.getDocument(), { width: 500, height: 500, registry });
    const s = layer(t.ctx, "F").state as Features;
    const lanes = Array.from({ length: s.count }, (_, k) => [s.start[k], s.lane[k]]).sort((x, y) => x[0] - y[0]);
    expect(lanes).toEqual([
      [0, 0],
      [40, 0],
      [58, 1],
      [120, 0],
      [200, 0],
    ]);
  });

  it("hits are polar: the radius picks the ring, the angle the position", async () => {
    const { t, l, angle, radius, entry } = await rings();
    const e = entry("F");
    const at = (k: number) => {
      const a = angle(k + 0.5);
      const r = radius(l.scrolledTop + l.entryTop[e] + l.entryHeight[e] / 2);
      return [l.ringCx + r * Math.cos(a), l.ringCy + r * Math.sin(a)];
    };
    const hit = Hit.create();
    for (const [k, item] of [
      [76, 1],
      [220, 0],
      [10, 0],
    ]) {
      const [x, y] = at(k);
      expect(t.ctx.layout.hitTest(t.section, x, y, hit), `${k}`).toBe(true);
      expect([hit.entry, hit.index, hit.kind, hit.item], `${k}`).toEqual([e, k, HoverKind.item, item]);
      // The point is the surface's again after the hit.
      expect([hit.x, hit.y]).toEqual([x, y]);
    }
    // The hub and outside the rings: no entry.
    expect(t.ctx.layout.hitTest(t.section, l.ringCx, l.ringCy, hit) && hit.entry).toBe(-1);
    expect(t.ctx.layout.hitTest(t.section, 2, 2, hit)).toBe(false);
  });

  it("reduces per angular column: 120,000 positions draw at most a sector per device column of the ring", async () => {
    const { t, l, entry } = await rings(120_000);
    const r = t.frame(recording());
    const e = entry("S");
    expect(e).toBeGreaterThanOrEqual(0);
    expect(r.sectorItems.n).toBeGreaterThan(100);
    expect(r.sectorItems.n).toBeLessThan(4 * l.ringW);
  });

  it("a focus on a whole cyclic unit is its shortest arc (start > stop across the origin)", async () => {
    const { t } = await rings();
    const set = t.ctx.graph.assembled.value?.set;
    const u = set ? set.units.key.findIndex((k) => set.strings.get(k) === "s/unit") : -1;
    const ext = Extents.create();
    Extents.add(ext, { kind: 1, id: u, row: -1 }, 200, 240);
    Extents.add(ext, { kind: 1, id: u, row: -1 }, 0, 30);
    expect(t.ctx.focus.fromExtents(ext, { section: t.section })).toBe(true);
    expect(t.ctx.focus.current.value.get(t.section)).toMatchObject({ block: 0, start: 200, stop: 30 });
  });

  it("marks are sectors across the rings: a selection and a highlight across the origin, from the ruler's ticks down", async () => {
    const { t, l } = await rings();
    t.ctx.selection.set({ systems: ["seqstar:s/unit"], runs: Int32Array.of(0, 225, 240, 0, 0, 12) }, "api");
    t.ctx.highlight.set("host", { systems: ["seqstar:s/unit"], runs: Int32Array.of(0, 230, 240, 0, 0, 6) });
    const r = t.frame(recording(), Dirty.overlay);
    const sel = passSectors(r, Pass.overlay).filter((s) => s[4] === SkinToken.selection);
    const hl = passSectors(r, Pass.overlay).filter((s) => s[4] === SkinToken.highlight);
    // Each the two runs of its range, as its arc: [225, 240) and [0, 12) (a cell is 1/240 of a turn).
    const arcs = (x: number[][]) => x.map((s) => [turn(s[2]), turn(s[3])].map((v) => Math.round(v * 240)));
    expect(arcs(sel).sort((a, b) => a[0] - b[0])).toEqual([
      [0, 12],
      [225, 240],
    ]);
    expect(arcs(hl).sort((a, b) => a[0] - b[0])).toEqual([
      [0, 6],
      [230, 240],
    ]);
    // From the ruler's ticks (not its numbers' room above them) to the last ring's inner edge.
    const last = l.entryCount - 1;
    const outer = l.ringR - (l.ringTop - l.scrolledTop) * l.ringScale;
    const inner = l.ringR - (l.entryTop[last] + l.entryHeight[last]) * l.ringScale;
    for (const s of [...sel, ...hl]) {
      expect(s[1]).toBeCloseTo(outer, 0);
      expect(s[0]).toBeCloseTo(inner, 0);
    }
    expect(outer).toBeLessThan(l.ringR - 10);
    // The selection's sides: radial lines at its arc's ends.
    const sides = passSectors(r, Pass.overlay).filter((s) => s[4] === SkinToken.selectionBorder);
    expect(sides.length).toBe(4);
  });

  it("the hover: a position as a radial line as wide at the hub as at the rim, an item as its sector and outline", async () => {
    // 120,000 positions: a cell is far below a pixel, so a hovered position is a line.
    const big = await rings(120_000);
    const e = big.entry("S");
    const rr = big.radius(big.l.scrolledTop + big.l.entryTop[e] + big.l.entryHeight[e] / 2);
    const a = big.angle(30_000.5);
    const hit = Hit.create();
    big.t.ctx.layout.hitTest(big.t.section, big.l.ringCx + rr * Math.cos(a), big.l.ringCy + rr * Math.sin(a), hit);
    big.t.ctx.hover.setFromHit(hit, "pointer");
    const lines = passPaths(big.t.frame(recording(), Dirty.overlay), Pass.overlay).filter(
      (p) => p.color === SkinToken.hover,
    );
    expect(lines.length).toBe(1);
    const [p0, p1] = lines[0].points;
    const polar = ([x, y]: number[]) => [
      Math.hypot(x - big.l.ringCx, y - big.l.ringCy),
      Math.atan2(y - big.l.ringCy, x - big.l.ringCx),
    ];
    // Radial: both ends at the position's angle (a quarter turn), from the ticks to the hub's ring; 2 px wide.
    expect(turn(polar(p0)[1])).toBeCloseTo(0.25, 3);
    expect(turn(polar(p1)[1])).toBeCloseTo(0.25, 3);
    expect(polar(p0)[0]).toBeGreaterThan(polar(p1)[0] + 50);
    expect(lines[0].width).toBe(4);
    // 240 positions and an item: no line; its range's sector, and its outline (the renderer's overlay) on its ring.
    const { t, l, angle, radius, entry } = await rings();
    const f = entry("F");
    const b = angle(76.5);
    const r = radius(l.scrolledTop + l.entryTop[f] + l.entryHeight[f] / 2);
    t.ctx.layout.hitTest(t.section, l.ringCx + r * Math.cos(b), l.ringCy + r * Math.sin(b), hit);
    t.ctx.hover.setFromHit(hit, "pointer");
    expect([hit.kind, hit.item]).toEqual([HoverKind.item, 1]);
    const rec = t.frame(recording(), Dirty.overlay);
    expect(passPaths(rec, Pass.overlay).filter((p) => p.color === SkinToken.hover)).toEqual([]);
    const lit = passSectors(rec, Pass.overlay).filter((s) => s[4] === SkinToken.hover);
    expect(lit.map((s) => [turn(s[2]), turn(s[3])].map((v) => Math.round(v * 240)))).toEqual([[40, 112]]);
    // The outline: stroked edges in the features' ring, drawn by the overlay pass through the ring map.
    const r0 = radius(l.scrolledTop + l.entryTop[f] + l.entryHeight[f]);
    const r1 = radius(l.scrolledTop + l.entryTop[f]);
    const outline = passSectors(rec, Pass.overlay).filter(
      (s) => s[4] !== SkinToken.hover && s[0] >= r0 - 1 && s[1] <= r1 + 1,
    );
    expect(outline.length).toBeGreaterThanOrEqual(4);
  });

  it("the ruler: unlabeled minor ticks between its numbers on rings, a fifth of the step apart; none linearly", async () => {
    const { t, l } = await rings();
    const r = t.frame(recording());
    // The ticks are the sectors in the ruler's ring; labeled every 10 (the numbers), minor every 2 (1-based: index 9 is
    // 10, index 11 is 12).
    const ruler = passSectors(r, Pass.main).filter(
      (s) => s[1] > l.ringR - l.entryHeight[0] * l.ringScale && s[3] - s[2] < TAU / 240,
    );
    const at = new Set(ruler.map((s) => Math.round(turn((s[2] + s[3]) / 2) * 240 - 0.5)));
    expect([at.has(9), at.has(11), at.has(19), at.has(2)]).toEqual([true, true, true, false]);
    const long = ruler.filter((s) => s[1] - s[0] > 2 * l.ringScale).length;
    expect([long, ruler.length - long]).toEqual([24, 96]);
    // Linear: the labeled ticks only.
    const lin = await composed(plasmid(240, null), { width: 1000 });
    const ticks = lin.frame(recording()).rectItems;
    let n = 0;
    for (let i = 0; i < ticks.n; i++) if (ticks.w[i] === 1 && ticks.h[i] === 4) n++;
    expect(n).toBeLessThan(30);
  });

  it("the keyboard cursor goes round a ring: its cell's sector outlined, the hover on it, past the last onto the first", async () => {
    const t = await driven(plasmid(), {
      renderers: KernelRenderers as LayerRenderer[],
      registry,
      size: { "v.main": [500, 500] },
    });
    const section = t.ctx.viewports.sections.value[0];
    t.focus(section);
    t.key(section, "End");
    expect(t.ctx.ui.keyCursor.value?.index).toBe(239);
    t.key(section, "ArrowRight");
    expect(t.ctx.ui.keyCursor.value).toMatchObject({ index: 0, header: false });
    t.key(section, "ArrowLeft");
    t.key(section, "ArrowLeft");
    expect(t.ctx.ui.keyCursor.value?.index).toBe(238);
    // The hover follows it (the cell's middle on the surface hits it).
    expect(t.ctx.hover.live.current.index).toBe(238);
    t.key(section, "Home");
    expect(t.ctx.ui.keyCursor.value).toMatchObject({ index: 0, header: false });
    const { SectionComposer } = await import("../src/draw/index.ts");
    const c = new SectionComposer(t.ctx, section);
    const rec = recording();
    c.frame(rec, Dirty.all);
    const ring = passSectors(rec, Pass.overlay).filter((s) => s[4] === SkinToken.focusRing);
    // Four edges just outside cell 0 (an arc of 1/240 turn from the top), on the ring of the track it is on (the ruler).
    expect(ring.length).toBe(4);
    const a0 = Math.min(...ring.map((s) => (s[2] + Math.PI / 2) / TAU));
    const a1 = Math.max(...ring.map((s) => (s[3] + Math.PI / 2) / TAU));
    expect([a0 < 0 && a0 > -0.5 / 240, a1 > 1 / 240 && a1 < 1.5 / 240]).toEqual([true, true]);
    c.dispose();
  });

  it("a type with no ring form (a matrix's image tiles) is drawn as bars, and the section says so once", async () => {
    const b = createSeqViewSpec();
    const s = b.data.sequence({ ref: "s", value: "MKTAYIAKQR".repeat(10) });
    const index = Array.from({ length: 100 }, (_, i) => i);
    const value = { index, value: index.map((i) => [i % 3, (i * 7) % 5, 1]) };
    const m = b.data.annotation({ ref: "m", geometry: "matrix", target: s, states: ["x", "y", "z"], value } as never);
    b.view({ ref: "v", label: "V", style: "circular" })
      .block({ sequence: s })
      .track({ label: "M" })
      .matrix({ data: m });
    const t = await composed(b.getDocument(), { width: 400, height: 400, registry });
    const said: string[] = [];
    t.ctx.events.issue.subscribe((i) => said.push(`${i.code} ${i.severity} ${i.message}`));
    const r = t.frame(recording());
    t.frame(recording());
    expect(said).toEqual([
      expect.stringMatching(
        /^svs\.render\.fallback warning layer "v\.matrix[^"]*" \(matrix\) has no ring form in section/,
      ),
    ]);
    expect(passSectors(r, Pass.main).length).toBeGreaterThan(0);
  });
});
