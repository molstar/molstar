/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A `pairs` layer under a `positions` filter (seqspec:format#layer-types `pairs`). Case 10's "The pocket" filters the
// Chains view to STI's contacts, both sides (a SeqQL query node), so the STI block stays, every contact is placed whole
// and every shown position has one. Filtered to `pocket` alone (the residues within 4 Å of imatinib), the STI block is
// left out: the contacts' kinase ends still draw their glyphs (each pair half-placed), and the poses' contacts, whose
// ligands are in no block of the section, stay out as before.

import type { SeqViewSpecDoc } from "@molstar/seq-core/format";
import { describe, expect, it } from "vitest";
import { case10 } from "../../core/_test/format/cases.js";
import type { PairsState } from "../src/prepare.js";
import { drawn, layerOf, record } from "./draw-fixtures.js";

const CONTACTS = "Contacts (STI)";
async function chains(snapshot: string | null, doc: SeqViewSpecDoc = case10()) {
  const t = await drawn(doc);
  await t.ctx.snapshots.apply(snapshot);
  t.ctx.views.setCurrent("chains");
  await t.ctx.graph.update();
  for (const ref of t.ctx.viewports.sections.value) t.ctx.viewports.get(ref)?.setSize(1000, 400, 1);
  await t.ctx.graph.update();
  const { k } = layerOf(t.ctx, CONTACTS);
  const view = t.ctx.graph.assembled.value?.view;
  const shown = t.ctx.layout.axis("chains.main")?.axis?.shown;
  return { view, shown, state: view?.prepared[k].state as PairsState, rec: record(t.ctx, CONTACTS).rec };
}
/** Case 10 with "The pocket" filtered to `pocket` alone. */
function pocketOnly() {
  const doc = case10();
  const s2 = doc.snapshots?.find((s) => s.ref === "s2");
  if (s2) s2.ops = [{ op: "set", ref: "chains.main", params: { positions: "pocket" } }];
  return doc;
}

describe("pairs under a positions filter", () => {
  it("case 10's pocket keeps the STI block: every contact placed whole, a glyph per kinase end, a glyph per kind at STI", async () => {
    const pocket = await chains("s2");
    expect(pocket.view?.resolved.blocks.count).toBe(2);
    const s = pocket.state;
    expect(s.count).toBe(18);
    for (let k = 0; k < s.count; k++) expect([s.block[2 * k], s.block[2 * k + 1]]).toEqual([0, 1]);
    // All 18 contacts land on STI's one cell, which holds three kinds (hydrophobic, hbond, pi stacking): three glyphs there.
    expect(pocket.rec.glyphItems.n).toBe(18 + 3);
    // Only STI's contacts: every position shown is an end of one (no residue only a docked pose touches).
    const ends = (b: number) =>
      new Set(
        Array.from({ length: s.count * 2 }, (_, e) => e)
          .filter((e) => s.block[e] === b)
          .map((e) => s.start[e]),
      ).size;
    expect([pocket.shown?.[0], pocket.shown?.[1]]).toEqual([ends(0), ends(1)]);
  });

  it("filtered to `pocket` alone: the contacts' kinase ends draw, one glyph each, and no pose contact", async () => {
    const all = await chains("s1");
    const kinaseEnds = Array.from(all.state.block).filter((b, e) => e % 2 === 0 && b === 0).length;
    expect([all.state.count, all.rec.glyphItems.n]).toEqual([18, 36]);

    const pocket = await chains("s2", pocketOnly());
    expect(pocket.view?.resolved.blocks.count).toBe(1);
    const s = pocket.state;
    expect(s.count).toBe(kinaseEnds);
    for (let k = 0; k < s.count; k++) expect([s.block[2 * k], s.block[2 * k + 1]]).toEqual([0, -1]);
    // Only placed ends are indexed; every glyph is drawn.
    expect(s.endOf.length).toBe(s.count);
    expect(pocket.rec.glyphItems.n).toBe(kinaseEnds);
  });
});
