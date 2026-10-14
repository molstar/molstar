/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The size channels in the format (seqspec:format#layer-types, seqspec:format#themes): each type's channels are layer
// params holding a size theme (`uniform`, `rules`, `field_scale`), the themes resolve to the layer's `themes` after its
// color and glyph, and a `field_scale` over another field with no domain takes its `domain_group`'s.
import { createSeqViewSpec } from "@molstar/seq-core/builder";
import { describe, expect, it } from "vitest";
import { compile, resolveView, type SeqViewSpecDoc, supportedAll, validate } from "../../src/format/index.js";
import { files, fixtureLoaders, registries } from "./load-fixtures.js";

const reg = registries();
type B = ReturnType<typeof createSeqViewSpec>;

/** A sequence `s`, ranges `a` (with a `score`) and pairs `pp`, and what `more` adds. */
function build(more: (b: B, s: ReturnType<B["data"]["sequence"]>) => void) {
  const b = createSeqViewSpec();
  const s = b.data.sequence({ ref: "s", value: "MKVLAAGIVGLLLAQ" });
  b.data.annotation({
    ref: "a",
    geometry: "ranges",
    target: s,
    value: { start: [0, 3, 6], stop: [2, 5, 9], score: [1, 5, 9] },
  } as never);
  b.data.annotation({
    ref: "pp",
    geometry: "pairs",
    target: s,
    value: { start: [0], b_target: ["s"], b_start: [5] },
  } as never);
  more(b, s);
  return b.getDocument() as SeqViewSpecDoc;
}
const codes = (doc: SeqViewSpecDoc) => validate(doc, reg).map((i) => i.code);
const uniform = (size: number) => ({ theme: "uniform", params: { size } });

describe("size channels", () => {
  it("are params of the types that declare them, holding a size theme", () => {
    const ok = build((b) => {
      const t = b.view({ ref: "v", label: "V" }).track({ label: "T" });
      t.ranges({ data: "a", thickness: uniform(0.5), radius: uniform(4) } as never);
      t.pairs({
        data: "pp",
        width: { theme: "field_scale", params: { field: "score", range: [1, 8], scale: "sqrt" } },
      } as never);
      t.points({
        data: "a",
        radius: { theme: "rules", params: { rules: [{ when: [], size: 3, label: "all" }] } },
      } as never);
    });
    expect(codes(ok)).toEqual([]);
  });

  it("are refused where the type has no such channel, or the theme serves another kind, or `range` is missing", () => {
    const bad = (type: "letters" | "ranges", layer: Record<string, unknown>) =>
      codes(
        build((b, s) => {
          const t = b.view({ ref: "v", label: "V" }).block({ sequence: s }).track({ label: "T" });
          if (type === "letters") t.letters(layer as never);
          else t.ranges({ data: "a", ...layer } as never);
        }),
      );
    expect(bad("letters", { thickness: uniform(1) })).toEqual(["svs.params.unknown"]);
    expect(bad("ranges", { thickness: { theme: "letter" } })).toEqual(["svs.theme.channel"]);
    expect(bad("ranges", { thickness: { theme: "field_scale", params: { field: "score" } } })).toHaveLength(1);
    expect(bad("ranges", { color: { theme: "field_scale", params: { range: [1, 2] } } })).toEqual([
      "svs.theme.channel",
    ]);
  });

  it("resolve to the layer's themes after color and glyph; a field_scale takes its group's domain", async () => {
    const doc = build((b, s) => {
      b.data.annotation({
        ref: "o",
        geometry: "ranges",
        target: s,
        value: { start: [1], stop: [2], score: [20] },
      } as never);
      const t = b.view({ ref: "v", label: "V" }).track({ label: "T" });
      const radius = { theme: "field_scale", params: { field: "score", range: [2, 6] } };
      t.points({ ref: "x", data: "a", radius, domain_group: "g" } as never);
      t.points({ ref: "y", data: "o", radius, domain_group: "g" } as never);
      t.points({ ref: "z", data: "a", radius } as never);
    });
    const { set, issues } = await compile(doc, reg, fixtureLoaders(files()));
    expect(issues.filter((i) => i.severity === "error")).toEqual([]);
    const r = resolveView(doc, set, "v", reg, supportedAll(reg));
    const L = r.view.layers;
    const at = (ref: string) => Array.from(L.ref).indexOf(set.strings.find(ref));
    expect(L.themes[at("x")].map((t) => t?.theme)).toEqual(["field_palette", "uniform", "field_scale"]);
    expect(L.themes[at("x")][2]?.params).toMatchObject({ field: "score", domain: [1, 20] });
    expect(L.themes[at("y")][2]?.params).toMatchObject({ domain: [1, 20] });
    expect(L.themes[at("z")][2]?.params).toMatchObject({ domain: null });
    expect(set.strings.get(r.view.domains.field[0])).toBe("score");
  });
});
