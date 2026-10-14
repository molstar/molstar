/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// VL-06: emphasis shows text exactly where the letters under it do (`textShown`, seqspec:viewer#the-render-backend: a
// row never shows some columns' text and not others'): a short track draws no bold letters, only its cells.

import { createSeqViewSpec } from "@molstar/seq-core/builder";
import { describe, expect, it } from "vitest";
import { GRCH38 } from "../../core/_test/format/stream-provider.ts";
import { KernelRenderers } from "../src/index.ts";
import { LayerFrame, RecordingBackend } from "../src/testing/index.ts";
import { drawn, layerOf, texts } from "./draw-fixtures.ts";
import { GENOME, MaxInFlight, openGenome, streamedDoc } from "./stream-fixtures.ts";

/** Inside chr2, clear of the fixture's N runs. */
const AT = 100_040_100;

/** Texts drawn by the track's letters and emphasis layers at 12 px a cell, the track at its height hint. */
async function textsAt(hint: number) {
  const b = createSeqViewSpec();
  const a = b.data.sequence({ ref: "a", value: "MKVLATGHEK".repeat(8) });
  const d = b.data.annotation({ ref: "d", geometry: "ranges", target: a, value: { start: [5], stop: [10] } });
  const s = b.view({ ref: "v", label: "V" }).section({ ref: "v.s", height: "box" }).block({ sequence: a });
  const tr = s.track({ label: "L", ...(hint ? { height: hint } : {}) });
  tr.letters();
  tr.emphasis({ data: d });
  const t = await drawn(b.getDocument(), { width: 1000, height: 300 });
  const { k, section } = layerOf(t.ctx, "L");
  t.ctx.viewports.get(section)?.setView(12, 0, "api");
  await t.ctx.graph.update();
  const out: Record<string, number> = {};
  for (const [name, kk] of [
    ["letters", k],
    ["emphasis", k + 1],
  ] as const) {
    const f = new LayerFrame(t.ctx, t.ctx.registry.get("renderers", name) as never, section, kk);
    const r = new RecordingBackend();
    f.draw(r);
    out[name] = texts(r).length;
    out[`${name}Rects`] = r.rectItems.n;
  }
  return out;
}

describe("VL-06: emphasis text follows the letters' text", () => {
  it("at the default height both draw text: the letters in view, the five emphasized bold", async () => {
    const out = await textsAt(0);
    expect(out.letters).toBeGreaterThan(5);
    expect(out.emphasis).toBe(5);
  });

  it("a track under 10 px (hint 0.4, 8.8 px): neither draws text; emphasis marks its five cells", async () => {
    const out = await textsAt(0.4);
    expect([out.letters, out.emphasis]).toEqual([0, 0]);
    expect(out.emphasisRects).toBeGreaterThan(0);
  });
});

describe("VL-10: emphasis over a streamed block's letters underlines them, never covers them", () => {
  it("at 12 px a cell on chr2, the emphasized residues' marks are underlines only", async () => {
    const doc = streamedDoc({ sequence: true });
    type N = { kind: string; ref?: string; params?: Record<string, unknown>; children?: N[] };
    const root = doc.root as N;
    const data = (root.children ?? []).find((c) => c.kind === "data") as N;
    data.children?.push({
      kind: "annotation",
      ref: "hot",
      params: { geometry: "ranges", target: "hg38/chr2", value: { start: [AT + 10], stop: [AT + 15] } },
    });
    const walk = (n: N) => {
      if (n.kind === "track" && n.params?.label === "Sequence")
        n.children?.push({ kind: "layer", ref: "hot.emphasis", params: { type: "emphasis", data: "hot" } });
      n.children?.forEach(walk);
    };
    walk(root);
    const t = await openGenome({
      chromosomes: GRCH38.slice(0, 3),
      letters: true,
      doc,
      renderers: [...KernelRenderers],
    });
    t.ctx.config.set(MaxInFlight, 64);
    t.vp(GENOME).zoomTo(1, AT, AT + 60);
    for (let k = 0; k < 6; k++) {
      t.layout(GENOME);
      t.manager.pass();
      t.streams.genome.release();
      await new Promise((ok) => setTimeout(ok, 0));
    }
    const { k, section } = layerOf(t.ctx, "Sequence");
    const draw = (kk: number) => {
      const f = new LayerFrame(
        t.ctx,
        t.ctx.registry.get("renderers", kk === k ? "letters" : "emphasis") as never,
        section,
        kk,
      );
      const r = new RecordingBackend();
      f.draw(r);
      return r;
    };
    // The letters layer draws the held letters as text; emphasis, five underlines (2 px at DPR 1).
    expect(texts(draw(k)).length).toBeGreaterThan(20);
    const e = draw(k + 1).rectItems;
    const hs = Array.from(e.h.subarray(0, e.n));
    expect(hs.length).toBe(5);
    expect(new Set(hs)).toEqual(new Set([2]));
  });
});
