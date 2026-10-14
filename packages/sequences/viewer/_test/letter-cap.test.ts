/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Per-position text is never larger than the cap (seqspec:viewer#the-render-backend, per-position text): at the
// closest zoom, case 14's DNA and its protein's residues over their codons, case 15's wild type and its four designs
// placed through their mappings (with their emphases), and an alignment's rows are drawn
// at most at `maxLetterPx`, by default the header labels' size (`labelFontPx`); a host raises the cap, or moves it
// with the labels' size. The text layout's letters are the cap's size by default.

import { createSeqViewSpec } from "@molstar/seq-core/builder";
import type { SeqViewSpecDoc } from "@molstar/seq-core/format";
import { describe, expect, it } from "vitest";
import { case14, case15 } from "../../core/_test/format/cases.ts";
import { DefaultRenderProps, type RenderProps } from "../src/layout.ts";
import { LayerFrame, RecordingBackend } from "../src/testing/index.ts";
import { fontPx } from "../src/text.ts";
import { drawn, texts } from "./draw-fixtures.ts";

/** The renderers that draw per-position text. */
const PER_POSITION = new Set(["letters", "rows", "emphasis"]);

/** An alignment's rows. */
function alignmentRows(): SeqViewSpecDoc {
  const b = createSeqViewSpec();
  const aln = b.data.alignment({
    ref: "aln",
    value: ["MSEEKLKELAEKAIRL", "MSEE-LKELAEKAIRL", "MSDEKLKE-AEKAIRL"],
    flavor: "aligned",
  } as never);
  b.view({ ref: "v", label: "V", height: "box" }).block({ alignment: aln }).track({ label: "Rows" }).rows();
  return b.getDocument();
}

/**
 * Every section at its closest zoom (`maxCell`); the font px of each per-position text item drawn, per renderer, and
 * the text layout's font px of a wrapped section, if any.
 */
async function sizesAtFullZoom(doc: SeqViewSpecDoc, render?: Partial<RenderProps>) {
  const t = await drawn(doc, { width: 1000, height: 600, dpr: 2, render });
  for (const ref of t.ctx.viewports.sections.value) {
    const vp = t.ctx.viewports.get(ref);
    vp?.setView(vp.live.current.maxCell, 0, "api");
  }
  await t.ctx.graph.update();
  const m = t.ctx.graph.assembled.value?.view;
  if (!m) throw new Error("no view");
  const px: number[] = [];
  const renderers = new Set<string>();
  for (let k = 0; k < m.prepared.length; k++) {
    const name = m.prepared[k].renderer;
    if (!PER_POSITION.has(name)) continue;
    const entry = t.ctx.registry.get("renderers", name);
    const section = m.set.strings.get(m.resolved.sections.ref[m.resolved.tracks.section[m.resolved.layers.track[k]]]);
    if (!entry) continue;
    const rec = new RecordingBackend();
    new LayerFrame(t.ctx, entry, section, k).draw(rec);
    for (const d of texts(rec)) {
      // One character a position, a three-letter code three; longer text (a ligand's name) is not per position.
      if (d.text.length > 3) continue;
      px.push(fontPx(t.ctx.text.fontString(d.font)));
      renderers.add(name);
    }
  }
  return { px, renderers, props: t.ctx.layout.props(t.ctx.viewports.sections.value[0]) };
}

describe("per-position text is capped at the header labels' size", () => {
  for (const [name, doc] of [
    ["case 14", case14],
    ["case 15", case15],
    ["an alignment's rows", alignmentRows],
  ] as const)
    it(`${name}: no text above the cap at the closest zoom, and some at it`, async () => {
      const { px, renderers, props } = await sizesAtFullZoom(doc());
      expect(props.maxLetterPx).toBe(DefaultRenderProps.labelFontPx);
      expect(px.length, [...renderers].join()).toBeGreaterThan(10);
      expect(Math.max(...px)).toBe(props.labelFontPx);
    });

  it("a host raises the cap (`maxLetterPx`), or moves it with the labels' size", async () => {
    const raised = await sizesAtFullZoom(case15(), { maxLetterPx: 16 });
    expect(raised.props.maxLetterPx).toBe(16);
    // The row allows 15 (22 px × 0.72): above the labels' 12, under the cap.
    expect(Math.max(...raised.px)).toBe(15);
    const small = await sizesAtFullZoom(case15(), { labelFontPx: 10 });
    expect(small.props.maxLetterPx).toBe(10);
    expect(Math.max(...small.px)).toBe(10);
  });

  it("the text layout's letters are the cap's size unless given", async () => {
    const b = createSeqViewSpec();
    const s = b.data.sequence({ ref: "s", alphabet: "protein", value: "MSEEKLKELAEKAIRL" });
    b.view({ ref: "v", label: "V", layout: "text" }).block({ sequence: s }).track({ label: "P" }).letters();
    const doc = b.getDocument();
    const t = await drawn(doc, { render: { maxLetterPx: 14 } });
    expect(t.ctx.layout.props(t.ctx.viewports.sections.value[0]).wrappedFontPx).toBe(14);
    const own = await drawn(doc, { render: { maxLetterPx: 14, wrappedFontPx: 17 } });
    expect(own.ctx.layout.props(own.ctx.viewports.sections.value[0]).wrappedFontPx).toBe(17);
  });
});
