/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Hits by drawn geometry (HR2): a glyph wider than its column, a band thinner than its lane and an arc's stroke hit
// where they are drawn: 1 px inside a sized shape's edge is the item, 1 px outside is not (the plain position, the
// column no longer the item's); the tooltip lists the fields the layer's themes encode; and with `UI.OwnRowHover`
// "residue" a bar on an own unit's row hovers that unit's residue, its item kept for the card.

import { describe, expect, it } from "vitest";
import { Hit, HoverKind, SeqViewerConfig } from "../src/index.js";
import { hitGeometryDoc } from "./hit-geometry-fixtures.js";
import { lit } from "./overlay-fixtures.js";
import { CHAIN_AT, splitUnitsDoc } from "./split-units-fixtures.js";

const S = "v.main";

async function open() {
  const t = await lit(hitGeometryDoc(), { width: 900, height: 600 });
  // 3 px a position: an 8 px dot spans five columns.
  t.vp?.setView(3, 0, "api");
  t.frame();
  const ax = () => t.ctx.layout.axis(S);
  /** Position i's center x (CSS px). */
  const cx = (i: number) => {
    const a = ax();
    const k = a?.shownOf(0, i) ?? (0 as never);
    return ((a?.x(0, k) ?? 0) + (a?.x(0, (k + 1) as never) ?? 0)) / 2;
  };
  const hit = Hit.create();
  const at = (x: number, y: number) => {
    t.ctx.layout.hitTest(S, x, y, hit);
    return hit.kind === HoverKind.item ? `item ${hit.item}` : hit.kind === HoverKind.seq ? `seq ${hit.index}` : "none";
  };
  return { ...t, cx, at, hit };
}

describe("hits by drawn geometry", () => {
  it("points: 1 px inside a sized dot's edge is its item, 1 px outside the plain position", async () => {
    const t = await open();
    const y = t.yOf(S, "t-dots");
    // Score 100: radius 8. Its center column and 7 px right (two columns over) are the dot; 9 px right is not.
    expect([t.at(t.cx(60), y), t.at(t.cx(60) + 7, y), t.at(t.cx(60) - 7, y)]).toEqual(["item 1", "item 1", "item 1"]);
    expect(t.at(t.cx(60) + 9, y)).toMatch(/^seq/);
    // Score 0: radius 2; its column still hits it (the fallback), the next column does not.
    expect([t.at(t.cx(20) + 1, y), t.at(t.cx(20) + 3.5, y)]).toEqual(["item 0", "seq 21"]);
  });

  it("ranges: a band's share of its lane is the item, the lane off the band is the position", async () => {
    const t = await open();
    const y = t.yOf(S, "t-bars");
    const l = t.layout(S);
    const e = Array.from({ length: l.entryCount }, (_, k) => k).find((k) => l.entryTrack[k] === t.trackIndex("t-bars"));
    const h = l.entryHeight[e as number];
    // Confidence 1: the whole lane; 0.2: a fifth of it, centered.
    expect([t.at(t.cx(60), y - h / 2 + 4), t.at(t.cx(20), y), t.at(t.cx(20), y - h / 2 + 4)]).toEqual([
      "item 1",
      "item 0",
      "seq 20",
    ]);
  });

  it("ranges: a boundary's sized glyph, 1 px inside and outside its radius", async () => {
    const t = await open();
    const y = t.yOf(S, "t-ticks");
    const edge = (i: number) => t.ctx.layout.axis(S)?.x(0, t.ctx.layout.axis(S)?.shownOf(0, i) ?? (0 as never)) ?? 0;
    // Weight 10: radius 8 (as the lane allows, 8 px here); the glyph sits on the edge before position 60.
    expect([t.at(edge(60) - 7, y), t.at(edge(60) + 7, y)]).toEqual(["item 2", "item 2"]);
    expect(t.at(edge(60) - 9.5, y)).toMatch(/^seq/);
  });

  it("the end boundary's glyph as points and stacked variants: hit from inside the last cell", async () => {
    const t = await open();
    t.vp?.setView(1, 0, "api");
    t.frame();
    const a = t.ctx.layout.axis(S);
    const end = a ? a.x(0, (a.shownOf(0, 399) + 1) as never) : 0;
    for (const track of ["t-tick-points", "t-tick-variants"]) {
      const l = t.layout(S);
      const e = Array.from({ length: l.entryCount }, (_, k) => k).find((k) => l.entryTrack[k] === t.trackIndex(track));
      const top = l.scrolledTop + l.entryTop[e as number];
      const hits = (x: number) => {
        const seen = new Set<string>();
        for (let y = top; y < top + l.entryHeight[e as number]; y += 0.5) seen.add(t.at(x, y));
        return seen;
      };
      expect(hits(end - 1).has("item 3"), track).toBe(true);
    }
  });

  it("pairs: an end's sized glyph, 1 px inside and outside its radius", async () => {
    const t = await open();
    t.frame();
    const l = t.layout(S);
    const e = Array.from({ length: l.entryCount }, (_, k) => k).find(
      (k) => l.entryTrack[k] === t.trackIndex("t-pairs"),
    );
    // The ends' line: a glyph's center on it, found by hitting the end's own column across the entry.
    const top = l.scrolledTop + l.entryTop[e as number];
    let line = -1;
    for (let y = top; y < top + l.entryHeight[e as number] && line < 0; y += 0.5)
      if (t.at(t.cx(40), y) === "item 1") line = y;
    expect(line).toBeGreaterThan(0);
    // Count 10: radius 8, at position 40; 7 px left is the end, 9.5 px left nothing of the pairs.
    for (let y = line; y < line + 8 && t.at(t.cx(40), y + 0.5) === "item 1"; y += 0.5) line = y;
    const mid = line - 3.5;
    expect(t.at(t.cx(40) - 7, mid)).toBe("item 1");
    expect(t.at(t.cx(40) - 9.5, mid)).not.toBe("item 1");
  });

  it("pairs: an arc hits within its stroke's half width (at least 4 px)", async () => {
    const t = await open();
    const l = t.layout(S);
    const e = Array.from({ length: l.entryCount }, (_, k) => k).find((k) => l.entryTrack[k] === t.trackIndex("t-arcs"));
    const top = l.scrolledTop + l.entryTop[e as number];
    // The px of the entry's column at each arc's middle that hit it: about its stroke (16 px wide), at least 8 (1 px).
    const tall = (i: number, item: string) => {
      let n = 0;
      for (let y = top; y < top + l.entryHeight[e as number]; y += 0.5) if (t.at(t.cx(i), y) === item) n += 0.5;
      return n;
    };
    const [thin, wide] = [tall(25, "item 0"), tall(165, "item 1")];
    expect(thin).toBeGreaterThanOrEqual(7);
    expect(thin).toBeLessThanOrEqual(9.5);
    expect(wide).toBeGreaterThanOrEqual(12);
    expect(wide).toBeLessThanOrEqual(18);
  });

  it("the tooltip lists the fields the layer's themes encode", async () => {
    const t = await open();
    t.at(t.cx(60), t.yOf(S, "t-dots"));
    t.ctx.hover.setFromHit(t.hit, "pointer");
    const label = t.ctx.labels.compute(t.ctx.hover.live.current);
    expect(label?.lines.filter((x) => x.key === "score").map((x) => x.value)).toEqual(["100"]);
  });
});

describe("UI.OwnRowHover", () => {
  it('"residue": a bar on an own unit\'s row hovers that unit\'s residue, the item kept; "item" the bar', async () => {
    const t = await lit(splitUnitsDoc(), { width: 1000, height: 400 });
    const hit = Hit.create();
    t.ctx.config.set(SeqViewerConfig.UI.OwnRowHover, "item");
    const x = t.xOf(t.section, 0, CHAIN_AT + 5);
    const y = t.yOf(t.section, "t-chains", 1);
    t.ctx.layout.hitTest(t.section, x, y, hit);
    expect([hit.kind, hit.annotation >= 0]).toEqual([HoverKind.item, true]);
    t.ctx.config.set(SeqViewerConfig.UI.OwnRowHover, "residue");
    t.ctx.layout.hitTest(t.section, x, y, hit);
    const unit = t.set.units.key.findIndex((k) => t.set.strings.get(k) === "B");
    expect([hit.kind, hit.seq.unit, hit.seq.index, hit.split, hit.annotation >= 0]).toEqual([
      HoverKind.seq,
      unit,
      5,
      1,
      true,
    ]);
    t.ctx.hover.setFromHit(hit, "pointer");
    // The residue's title (its unit's label), the item's lines after it: what the card and tooltip name.
    const label = t.ctx.labels.compute(t.ctx.hover.live.current);
    expect([label?.title, label?.lines.find((x) => x.key === "annotation")?.value]).toEqual(["chain B A6", "Observed"]);
  });
});
