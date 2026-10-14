/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// EN1: `domain_group` (seqspec:format#themes, shared scales): the union of the members' domains at resolve.
import { createSeqViewSpec } from "@molstar/seq-core/builder";
import { describe, expect, it } from "vitest";
import {
  compile,
  register,
  resolveView,
  type SeqViewSpecDoc,
  splitRows,
  supportedAll,
  validate,
} from "../../src/format/index.js";
import { files, fixtureLoaders, registries } from "./load-fixtures.js";

const reg = registries();
const all = supportedAll(reg);
const at = (n: number, f: (i: number) => number) => Array.from({ length: n }, (_, i) => f(i));

type Make = (b: ReturnType<typeof createSeqViewSpec>, s: string) => void;
async function resolve(build: Make) {
  const b = createSeqViewSpec();
  b.data.sequence({ ref: "s", value: "MKVLAAGIVGLLLAQ" });
  build(b, "s");
  const doc = b.getDocument() as SeqViewSpecDoc;
  expect(validate(doc, reg)).toEqual([]);
  const { set, issues } = await compile(doc, reg, fixtureLoaders(files()));
  expect(issues.filter((i) => i.severity === "error")).toEqual([]);
  const r = resolveView(doc, set, "v", reg, all);
  const L = r.view.layers;
  const index = (ref: string) => Array.from(L.ref).indexOf(set.strings.find(ref));
  return { ...r, set, L, index };
}
const values = (b: ReturnType<typeof createSeqViewSpec>, s: string, ref: string, hi: number) =>
  b.data.annotation({
    ref,
    geometry: "values",
    target: s,
    value: { index: at(15, (i) => i), value: at(15, (i) => (hi * i) / 14) },
  });
const palette = (l: number, r: Awaited<ReturnType<typeof resolve>>) =>
  (r.L.themes[l][0]?.params?.palette as { value_domain?: (number | null)[] } | undefined)?.value_domain;

describe("EN1: domain_group", () => {
  it("members share the union of their data's domains: the layers' rows, the default palette, the table", async () => {
    const r = await resolve((b, s) => {
      const [small, big] = [values(b, s, "small", 1), values(b, s, "big", 8)];
      const t = b.view({ ref: "v", label: "V" }).track({ label: "T" });
      t.color({ ref: "a", data: small, domain_group: "g" } as never);
      t.color({ ref: "b", data: big, domain_group: "g" } as never);
      t.color({ ref: "c", data: small });
    });
    const [a, bb, c] = ["a", "b", "c"].map(r.index);
    expect([r.L.domainMin[a], r.L.domainMax[a]]).toEqual([0, 8]);
    expect([r.L.domainMin[bb], r.L.domainMax[bb]]).toEqual([0, 8]);
    expect([palette(a, r), palette(bb, r)]).toEqual([
      [0, 8],
      [0, 8],
    ]);
    // A layer outside the group keeps its own domain.
    expect([r.L.domainMin[c], r.L.domainMax[c], r.L.domainGroup[c]]).toEqual([0, 1, -1]);
    expect(palette(c, r)).toEqual([0, 1]);
    expect(r.L.domainGroup[a]).toBe(r.L.domainGroup[bb]);
    expect(r.set.strings.get(r.L.domainGroup[a])).toBe("g");
    expect(r.view.domains.count).toBe(1);
    expect([r.set.strings.get(r.view.domains.field[0]), r.view.domains.min[0], r.view.domains.max[0]]).toEqual([
      "value",
      0,
      8,
    ]);
  });

  it("a signed member makes the default diverging for all; a lone member of a group changes nothing", async () => {
    const r = await resolve((b, s) => {
      const neg = b.data.annotation({
        ref: "neg",
        geometry: "values",
        target: s,
        value: { index: [0, 1], value: [-2, 3] },
      });
      const pos = values(b, s, "pos", 1);
      const t = b.view({ ref: "v", label: "V" }).track({ label: "T" });
      t.color({ ref: "a", data: neg, domain_group: "g" } as never);
      t.color({ ref: "b", data: pos, domain_group: "g" } as never);
      t.color({ ref: "alone", data: pos, domain_group: "solo" } as never);
    });
    const [a, bb, alone] = ["a", "b", "alone"].map(r.index);
    for (const l of [a, bb]) {
      expect(palette(l, r)).toEqual([-3, 3]);
      expect(r.L.themes[l][0]?.params?.palette).toMatchObject({ colors: "RdBu" });
    }
    expect(palette(alone, r)).toEqual([0, 1]);
    expect(r.view.domains.count).toBe(1);
  });

  it("an explicit value_domain wins in its layer, its null ends are the group's", async () => {
    const r = await resolve((b, s) => {
      const [small, big] = [values(b, s, "small", 1), values(b, s, "big", 8)];
      const t = b.view({ ref: "v", label: "V" }).track({ label: "T" });
      const color = (value_domain: (number | null)[]) => ({
        theme: "field_palette" as const,
        params: { palette: { kind: "continuous" as const, colors: "Blues", value_domain } },
      });
      t.color({ ref: "a", data: small, domain_group: "g", color: color([0, 20]) } as never);
      t.color({ ref: "b", data: big, domain_group: "g", color: color([null, null]) } as never);
    });
    // The explicit range stays in its layer; the other's `null` ends come from the viewer's `domain` (the group's).
    expect(palette(r.index("a"), r)).toEqual([0, 20]);
    expect(palette(r.index("b"), r)).toEqual([null, null]);
    expect([r.L.domainMin[r.index("b")], r.L.domainMax[r.index("b")]]).toEqual([0, 8]);
  });

  it("a member's explicit domain is the group's for all of its members", async () => {
    const r = await resolve((b, s) => {
      const [small, big] = [values(b, s, "small", 1), values(b, s, "big", 8)];
      const t = b.view({ ref: "v", label: "V" }).track({ label: "T" });
      t.values({ ref: "a", data: small, domain: [0, 100], domain_group: "g" } as never);
      t.values({ ref: "b", data: big, domain_group: "g" } as never);
    });
    for (const x of ["a", "b"]) expect([r.L.domainMin[r.index(x)], r.L.domainMax[r.index(x)]]).toEqual([0, 100]);
  });

  it("a palette over another field takes the union of that field's range over the members", async () => {
    const r = await resolve((b, s) => {
      const lo = b.data.annotation({
        ref: "lo",
        geometry: "ranges",
        target: s,
        value: { start: [0, 3], stop: [2, 5], score: [1, 2] },
      });
      const hi = b.data.annotation({
        ref: "hi",
        geometry: "ranges",
        target: s,
        value: { start: [0, 3], stop: [2, 5], score: [5, 9] },
      });
      const t = b.view({ ref: "v", label: "V" }).track({ label: "T" });
      const color = {
        theme: "field_palette",
        params: { field: "score", palette: { kind: "continuous", colors: "Blues", value_domain: [null, 12] } },
      };
      t.ranges({ ref: "a", data: lo, domain_group: "g", color } as never);
      t.ranges({ ref: "b", data: hi, domain_group: "g", color } as never);
      t.ranges({ ref: "c", data: lo, color } as never);
    });
    expect(palette(r.index("a"), r)).toEqual([1, 12]);
    expect(palette(r.index("b"), r)).toEqual([1, 12]);
    expect(palette(r.index("c"), r)).toEqual([null, 12]);
    expect([r.set.strings.get(r.view.domains.field[0]), r.view.domains.min[0], r.view.domains.max[0]]).toEqual([
      "score",
      1,
      9,
    ]);
  });
});

describe("EN1: a field provider's name as a split's `by`", () => {
  it("rows by the provider's values; a provider no row has a value of is reported by the layer", async () => {
    const own = registries();
    register(own, [
      {
        fields: [
          {
            name: "lab:parity",
            type: "string",
            applies: () => true,
            read: (a, out, strings) => {
              for (let r = 0; r < a.count; r++) (out as Int32Array)[r] = strings.intern(r % 2 ? "odd" : "even");
            },
          },
        ],
      },
    ]);
    const b = createSeqViewSpec();
    const s = b.data.sequence({ ref: "s", value: "MKVLAAGIVGLLLAQ" });
    const a = b.data.annotation({
      ref: "a",
      geometry: "ranges",
      target: s,
      value: { start: [0, 2, 4], stop: [1, 3, 5] },
    });
    b.view({ ref: "v", label: "V" }).track({ ref: "t" }).ranges({ ref: "l", data: a, split: "lab:parity" });
    const doc = b.getDocument() as SeqViewSpecDoc;
    const { set, issues } = await compile(doc, own, fixtureLoaders(files()));
    expect(issues.map((i) => i.message)).toEqual([]);
    const r = resolveView(doc, set, "v", own, all);
    expect(r.issues.map((i) => i.code)).toEqual([]);
    const ann = set.annotations[r.view.data[0]];
    if (!ann) throw new Error(`no annotation ${r.view.data.length} ${r.view.layers.count} ${set.annotations.length}`);
    const split = { by: "lab:parity", order: null };
    expect(splitRows(set, [ann], [null], split, own).values).toEqual(["even", "odd"]);
    expect(Array.from(splitRows(set, [ann], [null], split, own).rowOf[0])).toEqual([0, 1, 0]);
    // Without the registry the name is a column no row has.
    expect(splitRows(set, [ann], [null], split).values).toEqual(["-"]);
  });
});
