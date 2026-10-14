/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A patch to what a computed node names reruns it, as a fresh open of the patched document gives, one case per ref
// family (seqspec:viewer#updates-and-status; LB-25): a query's names (another provider's answer, LB-02), a `row_stats`'
// `items`, a `paired`'s `pairs`, and a theme's computed inputs (`conservation_shading`'s `data`, `rules`' `profile`),
// whose theme is made again; an unrelated computed node's rerun leaves the theme. The graph's edges come from core's
// `nodeRefs`: a family's edge dropped fails its case.

import { createSeqViewSpec } from "@molstar/seq-core/builder";
import {
  compile,
  RequiredField,
  ref,
  type SeqViewSpecDoc,
  SimpleParamsSchema,
  type StageContext,
} from "@molstar/seq-core/format";
import type { AnnotationInput, SeqSet } from "@molstar/seq-core/model";
import { describe, expect, it } from "vitest";
import { fixtureLoaders } from "../../core/_test/format/load-fixtures.ts";
import { Q } from "../../core/src/query/index.ts";
import { createSeqViewer, DefaultRegistry, DefaultSeqViewerSpec, manualScheduler } from "../src/index.ts";
import { fixtureMeasurer } from "../src/testing/index.ts";
import { patch, refsRan, viewer } from "./graph-fixtures.ts";

const ser = (set: SeqSet | undefined, key: string) => {
  const a = set?.annotations.find((x) => x && set.strings.get(x.key) === key);
  if (!a) return null;
  return {
    count: a.count,
    start: Array.from(a.start),
    stop: Array.from(a.stop),
    sparse: a.sparse ? Array.from(a.sparse) : null,
    dense: a.dense ? a.dense.map((d) => Array.from(d)) : null,
    fields: a.fields.map((f) => JSON.stringify(f, (_k, v) => (ArrayBuffer.isView(v) ? Array.from(v as never) : v))),
  };
};

// A provider whose options name a sequence: a ranges item [start, stop) on it.
const span = {
  name: "test:span",
  options: SimpleParamsSchema({
    of: RequiredField(ref("sequence"), ""),
    start: RequiredField({ name: "int", check: () => {} } as never, ""),
    stop: RequiredField({ name: "int", check: () => {} } as never, ""),
  }),
  geometry: "ranges" as const,
  async compute(_ctx: StageContext, o: unknown) {
    const { of, start, stop } = o as { of: string; start: number; stop: number };
    return { key: "", geometry: "ranges", columns: { target: [of], start: [start], stop: [stop] } } as AnnotationInput;
  },
};
const make = () =>
  createSeqViewer(
    { ...DefaultSeqViewerSpec(), registry: [DefaultRegistry, { svs: { providers: [span] } } as never] },
    { scheduler: manualScheduler(), loaders: fixtureLoaders(), text: fixtureMeasurer() },
  );

describe("a query naming a computed node of another provider", () => {
  const doc = (start: number, stop: number): SeqViewSpecDoc => {
    const b = createSeqViewSpec();
    b.data.sequence({ ref: "p", value: "MACAACAACAAC" });
    b.data.computed({ ref: "sp", provider: "test:span", options: { of: "p", start, stop } } as never);
    b.data.computed({ ref: "q", provider: "query", query: Q.read("(and sp (letters C))") } as never);
    b.view({ ref: "v", label: "V" }).track({ label: "L" }).letters();
    return b.getDocument();
  };
  it("fresh open: the viewer's answer is compile's", async () => {
    const ctx = make();
    const r = await ctx.setDocument(doc(0, 6));
    await ctx.graph.update();
    const compiled = await compile(doc(0, 6), ctx.registry.svs, fixtureLoaders());
    const issues = r.issues.map((i) => `${i.severity} ${i.code} ${i.message}`);
    expect({ v: ser(ctx.graph.assembled.value?.set, "q"), issues }).toEqual({
      v: ser(compiled.set, "q"),
      issues: compiled.issues.map((i) => `${i.severity} ${i.code} ${i.message}`),
    });
  });
  it("a patch to the computed node reruns the query, as a fresh open gives", async () => {
    const ctx = make();
    await ctx.setDocument(doc(0, 6));
    await ctx.graph.update();
    await ctx.applyPatch(patch([{ op: "set", ref: "sp", params: { options: { of: "p", start: 6, stop: 12 } } }]));
    await ctx.graph.update();
    const fresh = make();
    await fresh.setDocument(doc(6, 12));
    await fresh.graph.update();
    expect(ser(ctx.graph.assembled.value?.set, "q")).toEqual(ser(fresh.graph.assembled.value?.set, "q"));
  });
});

describe("computed nodes rerun when what they name changes", () => {
  it("row_stats `items`: a patch to the annotation it counts reruns it", async () => {
    const mk = (start: number) => {
      const b = createSeqViewSpec();
      b.data.alignment({ ref: "aln", value: ["MKVLA", "MKVLA", "MRV-A"], alphabet: "protein" } as never);
      b.data.annotation({
        ref: "site",
        geometry: "ranges",
        target: "aln/0",
        numbering: "ordinal",
        value: { start: [start], stop: [start + 1] },
      } as never);
      b.data.computed({ ref: "st", provider: "row_stats", alignment: "aln", items: "site" } as never);
      const v = b.view({ ref: "v", label: "V", style: "msa" } as never);
      v.block({ alignment: "aln" } as never);
      v.track({ ref: "v.t" }).rows({} as never);
      return b.getDocument();
    };
    const w = viewer();
    expect((await w.ctx.setDocument(mk(1))).issues.filter((i) => i.severity === "error")).toEqual([]);
    const before = ser(w.ctx.graph.assembled.value?.set, "st");
    await w.ctx.applyPatch(patch([{ op: "set", ref: "site", params: { value: { start: [3, 4], stop: [4, 5] } } }]));
    const f = viewer();
    const d = mk(1);
    // the same patch, applied by hand to a fresh document
    const site = (d.root.children?.[0].children ?? []).find((c) => c.ref === "site");
    (site as { params: Record<string, unknown> }).params.value = { start: [3, 4], stop: [4, 5] };
    await f.ctx.setDocument(d);
    const got = ser(w.ctx.graph.assembled.value?.set, "st");
    expect({ ran: refsRan(w.ctx).includes("st"), got }).toEqual({
      ran: true,
      got: ser(f.ctx.graph.assembled.value?.set, "st"),
    });
    expect(got).not.toEqual(before);
  });

  it("paired `pairs`: a patch to the pairs annotation reruns it", async () => {
    const mk = () => {
      const b = createSeqViewSpec();
      b.data.sequence({ ref: "p", value: "MKVLAGWT" });
      b.data.sequence({ ref: "q", value: "MRVLAGWT" });
      b.data.annotation({
        ref: "contacts",
        geometry: "pairs",
        numbering: "ordinal",
        value: { target: ["p"], beg: [1], b_target: ["q"], b_beg: [2] },
      } as never);
      b.data.computed({ ref: "iface", provider: "paired", pairs: ["contacts"], partner: null } as never);
      b.view({ ref: "v", label: "V" }).track({ label: "L" }).letters();
      return b.getDocument();
    };
    const w = viewer();
    const r = await w.ctx.setDocument(mk());
    expect(r.issues.filter((i) => i.severity === "error")).toEqual([]);
    const value = { target: ["p", "p"], beg: [3, 5], b_target: ["q", "q"], b_beg: [4, 6] };
    await w.ctx.applyPatch(patch([{ op: "set", ref: "contacts", params: { value } }]));
    const d = mk();
    const c = (d.root.children?.[0].children ?? []).find((x) => x.ref === "contacts");
    (c as { params: Record<string, unknown> }).params.value = value;
    const f = viewer();
    await f.ctx.setDocument(d);
    expect(ser(w.ctx.graph.assembled.value?.set, "iface")).toEqual(ser(f.ctx.graph.assembled.value?.set, "iface"));
  });
});

describe("themes over a computed node", () => {
  const msa = (theme: Record<string, unknown>) => {
    const b = createSeqViewSpec();
    b.data.alignment({ ref: "msa", value: ["MKV-A", "MKVLA", "MRV-A"], alphabet: "protein" } as never);
    b.data.computed({ ref: "cons", provider: "column_stats", alignment: "msa", measure: "conservation" } as never);
    b.data.computed({ ref: "prof", provider: "profile", alignment: "msa" } as never);
    const v = b.view({ ref: "v", label: "V", style: "msa" } as never);
    v.block({ alignment: "msa" } as never);
    v.track({ ref: "v.t" }).letters({ ref: "v.l", color: theme } as never);
    return b.getDocument();
  };
  it("conservation_shading `data`: remade when the column reruns", async () => {
    const w = viewer();
    const r = await w.ctx.setDocument(msa({ theme: "conservation_shading", params: { data: "cons" } }));
    expect(r.issues.filter((i) => i.severity === "error")).toEqual([]);
    const made = w.created.filter((c) => c.endsWith(":v.l")).length;
    await w.ctx.applyPatch(patch([{ op: "set", ref: "cons", params: { measure: "gap_fraction" } } as never]));
    expect(refsRan(w.ctx)).toContain("cons");
    expect(w.created.filter((c) => c.endsWith(":v.l")).length).toBe(made + 1);
  });
  it("rules `profile`: remade when the profile reruns", async () => {
    const w = viewer();
    const rules = [{ when: [{ measure: { share: "letter" }, min: 0.6 }], color: "#ff0000" }];
    const r = await w.ctx.setDocument(msa({ theme: "rules", params: { rules, profile: "prof" } }));
    expect(r.issues.filter((i) => i.severity === "error")).toEqual([]);
    const made = w.created.filter((c) => c === "rules:v.l").length;
    await w.ctx.applyPatch(patch([{ op: "set", ref: "prof", params: { gaps: "ignore" } } as never]));
    expect(refsRan(w.ctx)).toContain("prof");
    expect(w.created.filter((c) => c === "rules:v.l").length).toBe(made + 1);
  });
  it("an unrelated computed node reruns: the theme is not remade", async () => {
    const w = viewer();
    const rules = [{ when: [{ measure: { column: "cons" }, min: 0.6 }], color: "#ff0000" }];
    await w.ctx.setDocument(msa({ theme: "rules", params: { rules } }));
    const made = w.created.filter((c) => c === "rules:v.l").length;
    await w.ctx.applyPatch(patch([{ op: "set", ref: "prof", params: { gaps: "ignore" } } as never]));
    expect(refsRan(w.ctx)).toContain("prof");
    expect(w.created.filter((c) => c === "rules:v.l").length).toBe(made);
  });
});
