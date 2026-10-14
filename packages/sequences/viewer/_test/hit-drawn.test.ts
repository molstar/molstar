/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// GG-01 and GG-02: hits by drawn geometry (HR2; seqspec:interaction#hover). Every device px of the data area hits the cell drawn there (`DrawView.xd`'s
// device-snapped edges), at the px's left edge (where a mouse reports), its middle and its right edge; linear on plain,
// reversed, cyclic from origin 30 and filtered blocks, and two blocks side by side, at a fractional cell and scroll and
// DPR 1, 1.5 and 2, and wrapped lines with filter breaks (a select drag's position too). A boundary's glyph on a block's
// first edge hits on its outer half too, in the gap before the block, naming the block's first position.

import { createSeqViewSpec } from "@molstar/seq-core/builder";
import { describe, expect, it } from "vitest";
import type { ShownIndex } from "../src/axis.js";
import { wrappedPositionNear } from "../src/hit.js";
import { Hit, HoverKind } from "../src/index.js";
import { recording } from "./composer-fixtures.js";
import { opened } from "./layout-fixtures.js";
import { lit } from "./overlay-fixtures.js";

const S = "v.s";
type Kind = "plain" | "reversed" | "origin30" | "filtered" | "two" | "tworev";
/** 400 positions (150 kept where filtered): wider than the view at the test's zoom, so it scrolls to a fraction. */
function doc(kind: Kind) {
  const b = createSeqViewSpec();
  const p = b.data.sequence({ ref: "p", value: "ACGTACGTTA".repeat(40), alphabet: "dna", cyclic: kind === "origin30" });
  const start = Array.from({ length: 150 }, (_, i) => 2 * i + (i % 3));
  b.data.annotation({ ref: "keep", geometry: "ranges", target: p, value: { start, stop: start.map((s) => s + 1) } });
  const sec = b.view({ ref: "v", label: "V" }).section({
    ref: S,
    ...(kind === "reversed" || kind === "tworev" ? { direction: -1 } : {}),
    ...(kind === "filtered" ? { positions: "keep", direction: -1 } : {}),
  });
  // Two blocks side by side: the second's first edge at a fraction of a device px (LC-26).
  if (kind === "two" || kind === "tworev") {
    sec.block({ sequence: p, range: [0, 37] });
    sec.block({ sequence: p, range: [37, 400] });
  } else sec.block({ sequence: p, ...(kind === "origin30" ? { origin: 30 } : {}) });
  sec.track({ ref: "t" }).letters();
  return b.getDocument();
}

describe("GG-01: a device px hits the cell drawn on it", () => {
  for (const kind of ["plain", "reversed", "origin30", "filtered", "two", "tworev"] as const)
    for (const dpr of [1, 1.5, 2])
      it(`linear, ${kind}, DPR ${dpr}`, async () => {
        const t = await opened(doc(kind));
        const vp = t.vp(S);
        vp.setSize(1000, 300, dpr);
        // A zoom that leaves a fractional cell and scroll, as a wheel zoom does.
        vp.setView(13.37, 101.3, "api");
        const v = t.ctx.layout.axis(S);
        const ax = vp.axis;
        if (!v || !ax) throw new Error("no axis");
        expect(v.scrollDev).toBeGreaterThan(0);
        const y = t.yOf(S, "t");
        const h = Hit.create();
        const bad: string[] = [];
        let total = 0;
        const edge = (b: number, k: number) =>
          Math.round((v.blockX0[b] + k * v.span[b] * v.cell) * v.dpr) - v.scrollDev;
        const lo = Math.ceil(t.layout(S).dataLeft * dpr);
        const hi = Math.floor(t.layout(S).dataRight * dpr);
        for (let b = 0; b < ax.count; b++)
          for (let k = 0; k < ax.shown[b]; k++)
            for (let p = Math.max(lo, edge(b, k)); p < Math.min(hi, edge(b, k + 1)); p++)
              for (const f of [0, 0.5, 0.999]) {
                t.ctx.layout.hitTest(S, (p + f) / dpr, y, h);
                total++;
                if (h.block !== b || h.shown !== k)
                  bad.push(`px ${p}+${f}: drawn ${b}/${k}, hit ${h.block}/${h.shown}`);
              }
        expect(total).toBeGreaterThan(2000);
        expect(bad).toEqual([]);
      });

  it("below a device px a cell, the px's last cell in closed form (LC-29), at DPR 1, 1.25 and 2", async () => {
    const b = createSeqViewSpec();
    const p = b.data.sequence({ ref: "p", value: "ACGT".repeat(250_000), alphabet: "dna" });
    b.view({ ref: "v", label: "V" }).section({ ref: S }).block({ sequence: p }).track({ ref: "t" }).letters();
    const t = await opened(b.getDocument());
    for (const dpr of [1, 1.25, 2]) {
      const vp = t.vp(S);
      vp.setSize(1000, 300, dpr);
      vp.setView(0.0137, 3.3, "api");
      const v = t.ctx.layout.axis(S);
      if (!v) throw new Error("no axis");
      const edge = (k: number) => Math.round((v.blockX0[0] + k * v.span[0] * v.cell) * v.dpr) - v.scrollDev;
      const bad: string[] = [];
      for (let pd = Math.ceil(t.layout(S).dataLeft * dpr); pd < Math.floor(t.layout(S).dataRight * dpr); pd++) {
        // The last cell whose left edge is on or before pd, by bisection.
        let [lo, hi] = [0, 999_999];
        while (lo < hi) {
          const mid = (lo + hi + 1) >> 1;
          if (edge(mid) <= pd) lo = mid;
          else hi = mid - 1;
        }
        // From any hint, the far one too.
        for (const hint of [0, lo, 999_999]) if (v.drawnShown(0, hint, pd) !== lo) bad.push(`pd ${pd} from ${hint}`);
      }
      expect(bad, `DPR ${dpr}`).toEqual([]);
    }
  });

  it("wrapped, filtered and reversed: each line's cells at DPR 1, 1.5 and 2, hit and taken by a select drag", async () => {
    const b = createSeqViewSpec();
    const p = b.data.sequence({ ref: "p", value: "ACGTACGTTA".repeat(30), alphabet: "dna" });
    const start = Array.from({ length: 120 }, (_, i) => 2 * i + (i % 3));
    b.data.annotation({ ref: "keep", geometry: "ranges", target: p, value: { start, stop: start.map((s) => s + 1) } });
    const sec = b.view({ ref: "v", label: "V", layout: "text" }).section({ ref: S, positions: "keep", direction: -1 });
    sec.block({ sequence: p });
    sec.track({ ref: "t" }).letters();
    const t = await opened(b.getDocument());
    const near = new Int32Array(4);
    const at = new Float64Array(2);
    for (const dpr of [1, 1.5, 2]) {
      t.vp(S).setSize(700, 600, dpr);
      await t.ctx.graph.update();
      const l = t.layout(S);
      const v = t.ctx.layout.axis(S);
      if (!v) throw new Error("no axis");
      const h = Hit.create();
      const bad: string[] = [];
      let total = 0;
      for (let e = 0; e < l.entryCount; e++) {
        const line = l.entryRow[e];
        if (l.entryTrack[e] < 0 || line < 0 || line >= l.lines.count) continue;
        const y = l.scrolledTop + l.entryTop[e] + l.entryHeight[e] / 2;
        v.lineStart = l.lines.start[line];
        v.lineStop = l.lines.stop[line];
        for (let k = v.lineStart; k < v.lineStop; k++) {
          const a = Math.round(v.x(0, k as ShownIndex) * dpr);
          const c = Math.round((v.x(0, k as ShownIndex) + (v.runs[0]?.advance ?? 0)) * dpr);
          for (let px = a; px < c; px++)
            for (const f of [0, 0.5, 0.999]) {
              t.ctx.layout.hitTest(S, (px + f) / dpr, y, h);
              total++;
              if (h.shown !== k) bad.push(`line ${line} px ${px}+${f}: drawn ${k}, hit ${h.shown}`);
              // A select drag's press and end alike, on the cell drawn there (LC-27).
              at[0] = (px + f) / dpr;
              at[1] = y;
              if (wrappedPositionNear(l, v, at, near) && (near[2] !== k || near[3] !== 0))
                bad.push(`line ${line} px ${px}+${f}: drawn ${k}, near ${near[2]} past ${near[3]}`);
            }
        }
      }
      expect(total).toBeGreaterThan(300);
      expect(bad, `DPR ${dpr}`).toEqual([]);
    }
  });
});

describe("GG-02: a glyph overhanging a block's first edge into the gap before it hits on both halves", () => {
  const N = 40;
  function ends(direction: 1 | -1) {
    const b = createSeqViewSpec();
    const p = b.data.sequence({ ref: "p", value: "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRV" });
    const data = b.data.annotation({
      ref: "ends",
      geometry: "ranges",
      target: p,
      value: { boundary: [0, N], label: ["start", "end"] },
    });
    const sec = b.view({ ref: "v", label: "V" }).section({ ref: "v.main", direction });
    sec.block({ sequence: p, range: [0, 20] });
    sec.block({ sequence: p });
    sec.track({ ref: "t-ranges", label: "ranges", height: 2 }).ranges({ data });
    sec.track({ ref: "t-points", label: "points", height: 2 }).points({ data });
    sec.track({ ref: "t-variants", label: "variants", height: 2 }).variants({ data });
    return b.getDocument();
  }

  for (const direction of [1, -1] as const)
    it(`direction ${direction}: the boundary on block 1's first edge, 4 and 2 px either side`, async () => {
      const t = await lit(ends(direction), { width: 1400, height: 400 });
      t.vp?.setView(15, 0, "api");
      t.frame(recording());
      const a = t.ctx.layout.axis("v.main");
      if (!a) throw new Error("no axis");
      const x1 = a.x(1, 0 as ShownIndex);
      const h = Hit.create();
      const l = t.layout("v.main");
      // Forward, block 1's first edge holds the N-terminal boundary (item 0); reversed, the C-terminal (item 1).
      const want = direction > 0 ? 0 : 1;
      for (const k of ["ranges", "points", "variants"]) {
        const e = Array.from({ length: l.entryCount }, (_, q) => q).find(
          (q) => l.entryTrack[q] === t.trackIndex(`t-${k}`),
        ) as number;
        const top = l.scrolledTop + l.entryTop[e];
        for (const dx of [-4, -2, 2, 4]) {
          const seen = new Set<number>();
          const where = new Set<string>();
          for (let y = top; y < top + l.entryHeight[e]; y += 0.5) {
            t.ctx.layout.hitTest("v.main", x1 + dx, y, h);
            if (h.kind !== HoverKind.item) continue;
            seen.add(h.item);
            where.add(`${h.block}/${h.shown}/${h.index}`);
          }
          expect([...seen], `${k} at ${dx}`).toEqual([want]);
          // On either half, block 1's first position (LC-28: the outer half named its last).
          expect([...where], `${k} at ${dx}`).toEqual([`1/0/${a.indexOf(1, 0 as ShownIndex)}`]);
        }
      }
    });
});
