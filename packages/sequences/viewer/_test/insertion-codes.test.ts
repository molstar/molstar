/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */
// Insertion codes in the viewer (seqspec:viewer#labels): a heavy chain in Kabat-style author numbers (52A, 82A-82C,
// 100A-100B) as a document's `labels` numbering. The ruler's ticks fall on plain round numbers only, each once, in
// order; the hover names a coded residue by its code; a document's ranges and a link's refs name residues by codes.

import { createSeqViewSpec } from "@molstar/seq-core/builder";
import { OrderedSet } from "@molstar/seq-core/model";
import { describe, expect, it } from "vitest";
import { drawn, layerOf, record, texts } from "./draw-fixtures.js";
import { driven } from "./input-fixtures.js";

const VH =
  "EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS";
const LABELS: (number | string)[] = [];
for (let n = 1; LABELS.length < VH.length; n++) {
  LABELS.push(n);
  if (n === 52) LABELS.push("52A");
  if (n === 82) LABELS.push("82A", "82B", "82C");
  if (n === 100) LABELS.push("100A", "100B");
}
LABELS.length = VH.length;
const at = (label: string) => LABELS.findIndex((x) => String(x) === label);

function doc() {
  const b = createSeqViewSpec();
  const h = b.data.sequence({ ref: "h", value: VH }).numbering({ name: "auth", method: "labels", labels: LABELS });
  const cdr = b.data.annotation({
    ref: "cdr",
    geometry: "ranges",
    target: h,
    numbering: "auth",
    value: { beg: ["52A", "82A"], end: ["56", "83"] },
  });
  const v = b.view({ ref: "v", label: "V", numbering: "auth" });
  v.track({ label: "" }).ruler();
  v.track({ ref: "v.cdr", label: "CDR" }).ranges({ data: cdr });
  v.track({ ref: "v.letters", label: "Letters" }).letters();
  return b.getDocument();
}

describe("insertion codes in the viewer", () => {
  it.each([4, 8, 30])(
    "at %i px a residue, the ruler's ticks are plain round numbers, once each, in order",
    async (cell) => {
      const t = await drawn(doc(), { width: 1000, height: 200 });
      const { section } = layerOf(t.ctx, "");
      t.ctx.viewports.get(section)?.setView(cell, 0, "api");
      await t.ctx.graph.update();
      const ticks = texts(record(t.ctx, "").rec);
      expect(ticks.length).toBeGreaterThan(3);
      expect(ticks.every((d) => /^\d+$/.test(d.text))).toBe(true);
      const nums = ticks.map((d) => Number(d.text));
      expect(nums).toEqual([...nums].sort((a, b) => a - b));
      expect(new Set(nums).size).toBe(nums.length);
      // Each over its own residue, never a coded one: 55 is two residues further than 53 is (52A between).
      const ax = t.ctx.layout.axis(section);
      if (!ax) throw new Error("no axis");
      for (const d of ticks) {
        const k = ax.shownOf(0, at(d.text) as never);
        expect(Math.abs(d.x - (ax.x(0, k) + ax.x(0, (k + 1) as never)) / 2)).toBeLessThan(1);
      }
    },
  );

  it("the hover names a coded residue by its code", async () => {
    const t = await driven(doc());
    const S = layerOf(t.ctx, "").section;
    const title = (label: string) => {
      t.pointer("move", S, t.xOf(S, 0, at(label)), t.yOf(S, "v.letters"));
      return t.ctx.labels.current.value?.title;
    };
    const labels = ["52", "52A", "53", "82B", "100A"];
    expect(labels.map(title)).toEqual(labels.map((l) => `h ${VH[at(l)]}${l}`));
    expect(title("52A")).toBe("h P52A");
  });

  it("a document's ranges and a link's refs name residues by codes", async () => {
    const t = await driven(doc());
    const set = t.ctx.graph.assembled.value?.set;
    const cdr = set?.annotations.find((a) => a && set.strings.get(a.key).endsWith("cdr"));
    expect([Array.from(cdr?.start ?? []), Array.from(cdr?.stop ?? [])]).toEqual([
      [at("52A"), at("82A")],
      [at("56") + 1, at("83") + 1],
    ]);
    const loci = t.ctx.refs.toLoci("h", { beg: "82A", end: "82C", numbering: "auth" }) as {
      elements: { indices: OrderedSet }[];
    } | null;
    const indices = loci?.elements[0].indices;
    expect(indices && [OrderedSet.min(indices), OrderedSet.max(indices)]).toEqual([at("82A"), at("82C")]);
  });
});
