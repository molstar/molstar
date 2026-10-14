/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Data mode (seqspec:viewer#data-mode, seqspec:viewer#the-document; VS1): every validation resolves refs against the
// incoming data, a set or a resolver over the build inputs; the default view validates; one commit validates the
// document once and each snapshot once.

import type { Op, SeqViewSpecDoc } from "@molstar/seq-core/format";
import { type AnnotationInput, type SeqSet, SeqSets, type SetInput } from "@molstar/seq-core/model";
import { describe, expect, it } from "vitest";
import { dataKeys } from "../src/data.js";
import type { SeqIssue, SeqViewerContext } from "../src/index.js";
import { patch, smallDoc, viewer } from "./graph-fixtures.js";

const hide = (ref: string): Op => ({ op: "set", ref, params: { hidden: true } });
const ranges = (key: string, target: string, start = 0): AnnotationInput => ({
  key,
  geometry: "ranges",
  target,
  columns: { start: [start], stop: [start + 2] },
});
const inputs = (...annotations: AnnotationInput[]): SetInput => ({
  molecules: [
    { key: "p", residues: "MKVLATGHEK", ownUnit: true },
    { key: "q", residues: "MKVLA", ownUnit: true },
  ],
  annotations,
});
const errorsOf = (issues: readonly { severity: string; code: string }[]) =>
  issues.filter((i) => i.severity === "error").map((i) => i.code);
function issuesOf(ctx: SeqViewerContext) {
  const out: SeqIssue[] = [];
  ctx.events.issue.subscribe((i) => out.push(i));
  return out;
}
const trackRefs = (ctx: SeqViewerContext) => {
  const a = ctx.graph.assembled.value;
  return Array.from(a?.view?.resolved.tracks.ref ?? [], (r) => a?.set.strings.get(r));
};
const node = (doc: SeqViewSpecDoc | null, ref: string) => {
  const find = (n: SeqViewSpecDoc["root"]): SeqViewSpecDoc["root"] | undefined =>
    n.ref === ref ? n : n.children?.map(find).find((x) => x);
  return doc ? find(doc.root) : undefined;
};
/** Counts full validations: `validate` reads the issue codes once per call (registeredNamespaces). */
function validations(ctx: SeqViewerContext) {
  const svs = ctx.registry.svs as { list: (k: string) => unknown };
  const list = svs.list.bind(svs);
  let n = 0;
  svs.list = (k: string) => {
    if (k === "issueCodes") n++;
    return list(k);
  };
  return () => n;
}

describe("data mode validates against the incoming data", () => {
  it("track.hide, local, preview, undo and a snapshot all work over build inputs", async () => {
    const { ctx } = viewer();
    const issues = issuesOf(ctx);
    const r = await ctx.setData(
      { inputs: inputs(ranges("dom", "p/unit")) },
      {
        views: (b) => {
          const v = b.view({ ref: "v", label: "V" });
          v.block({ unit: "p/unit" });
          v.track({ ref: "v.dom" }).ranges({ data: "dom" });
          v.track({ ref: "v.letters" }).letters();
          b.snapshot({ ref: "s", ops: [hide("v.letters")] });
        },
      },
    );
    expect(errorsOf(r.issues)).toEqual([]);
    expect(ctx.snapshots.list.value.map((s) => s.valid)).toEqual([true]);
    expect(ctx.actions.run("track.hide", "v.dom")).toBe(true);
    await ctx.graph.update();
    expect(node(ctx.document.base.value, "v.dom")?.params?.hidden).toBe(true);
    expect((await ctx.applyPatch(patch([hide("v.letters")]), { layer: "local" })).ok).toBe(true);
    expect((await ctx.document.setPreview(patch([{ op: "set", ref: "v.dom", params: { label: "D" } }]))).ok).toBe(true);
    expect(await ctx.history.undo()).toBe(true);
    expect(node(ctx.document.base.value, "v.dom")?.params?.hidden).toBeUndefined();
    await ctx.snapshots.apply("s");
    expect(ctx.snapshots.current.value).toBe("s");
    expect(node(ctx.document.shown.value, "v.letters")?.params?.hidden).toBe(true);
    expect(issues.map((i) => i.code)).toEqual([]);
  });

  for (const mode of ["inputs", "set"] as const) {
    const data = (...a: AnnotationInput[]) =>
      mode === "inputs" ? { inputs: inputs(...a) } : { set: SeqSets.build(inputs(...a)).set as SeqSet };
    it(`${mode}: a second setData adding a key keeps a local Hide; one removing a key drops one patch, once`, async () => {
      const { ctx } = viewer();
      const issues = issuesOf(ctx);
      expect((await ctx.setData(data(ranges("dom", "p/unit"), ranges("site", "q/unit")))).ok).toBe(true);
      expect((await ctx.applyPatch(patch([hide("view.dom")]), { layer: "local" })).ok).toBe(true);
      expect((await ctx.applyPatch(patch([hide("view.site")]), { layer: "local" })).ok).toBe(true);
      const ruler = node(ctx.document.shown.value, "view.ruler");
      const r = await ctx.setData(data(ranges("dom", "p/unit"), ranges("site", "q/unit"), ranges("more", "p/unit", 4)));
      expect(errorsOf(r.issues)).toEqual([]);
      expect(trackRefs(ctx)).toEqual(["view.ruler", "view.letters", "view.dom", "view.site", "view.more"]);
      expect(node(ctx.document.shown.value, "view.dom")?.params?.hidden).toBe(true);
      expect(node(ctx.document.shown.value, "view.ruler")).toBe(ruler);
      await ctx.setData(data(ranges("dom", "p/unit"), ranges("more", "p/unit", 4)));
      expect(ctx.document.local.value.map((l) => l.patch.ops[0])).toEqual([hide("view.dom")]);
      expect(issues.map((i) => i.code)).toEqual(["seqv.patch.dropped"]);
    });
  }

  it("the resolver over build inputs names units, row refs and kinds as the set does", () => {
    const i: SetInput = {
      molecules: [
        { key: "p", residues: "MKV", ownUnit: true },
        { key: "lig", kind: "ligand", elements: ["C", "N"] },
      ],
      alignments: [{ key: "msa", gapped: ["MKV", "MK-"], flavor: "aligned" }],
    };
    const keys = dataKeys({ inputs: i });
    const set = SeqSets.build(i).set;
    const t = { kind: -1, id: -1, row: -1 };
    const u = { kind: -1, id: -1, row: -1 };
    // A row's index is its digits only, as core's (VS-27): `msa/1e0`, `msa/0x1`, `msa/1.0` name nothing.
    for (const k of ["p", "p/unit", "lig", "msa", "msa/1", "msa/2", "nothing", "msa/1e0", "msa/0x1", "msa/1.0"]) {
      expect(keys.byName(k, t)).toBe(set.byName(k, u));
      if (set.byName(k, u)) expect([t.kind, t.row]).toEqual([u.kind, u.row]);
    }
    expect(keys.molecules.map((m) => m?.kind)).toEqual(set.molecules.map((m) => m?.kind));
  });
});

describe("the default view validates", () => {
  it("a key slugged to a track ref, with a suffix past a collision; a key no ref can name gets no track", async () => {
    const { ctx } = viewer();
    const r = await ctx.setData({
      inputs: inputs(ranges("my feats", "p/unit"), ranges("p/feats", "p/unit"), ranges("ruler", "q/unit")),
    });
    expect(errorsOf(r.issues)).toEqual([]);
    expect(trackRefs(ctx)).toEqual(["view.ruler", "view.letters", "view.p-feats", "view.ruler-2"]);
    // Inputs without molecules (an alignment alone) are valid.
    expect((await viewer().ctx.setData({ inputs: { alignments: [] } as unknown as SetInput })).ok).toBe(true);
  });

  const msa = { key: "msa", gapped: ["MKV-A", "MK-LA"], flavor: "aligned" as const };
  const meta: AnnotationInput = {
    key: "meta",
    geometry: "rows",
    target: "msa",
    columns: { row: [0, 1], species: ["human", "mouse"], score: [1, 2] },
  };

  it("an alignment alone: valid, and one rows track with the rows annotation's fields", async () => {
    const { ctx } = viewer();
    const r = await ctx.setData({ inputs: { alignments: [msa], annotations: [meta] } as unknown as SetInput });
    expect(errorsOf(r.issues)).toEqual([]);
    expect(trackRefs(ctx)).toEqual(["view.ruler", "view.letters", "view.rows"]);
    expect(node(ctx.document.base.value, "view.rows")?.children?.[0]?.params?.data_columns).toEqual([
      { data: "meta", field: "species" },
      { data: "meta", field: "score" },
    ]);
    const set = SeqSets.build({ molecules: [], alignments: [msa], annotations: [meta] }).set;
    const { ctx: c2 } = viewer();
    expect(errorsOf((await c2.setData({ set })).issues)).toEqual([]);
    expect(node(c2.document.base.value, "view.rows")?.children?.[0]?.params?.data_columns).toEqual([
      { data: "meta", field: "species" },
      { data: "meta", field: "score" },
    ]);
  });

  it("an MSA with row molecules and row metadata: no errors and no empty rows track", async () => {
    const { ctx } = viewer();
    const r = await ctx.setData({
      inputs: {
        molecules: [
          { key: "r0", residues: "MKVA", ownUnit: true },
          { key: "r1", residues: "MKLA", ownUnit: true },
        ],
        alignments: [msa],
        annotations: [meta],
      },
    });
    expect(errorsOf(r.issues)).toEqual([]);
    expect(trackRefs(ctx)).toEqual(["view.ruler", "view.letters"]);
  });
});

describe("validation per click", () => {
  it("a committed Hide validates once and each snapshot once; undo and redo the same", async () => {
    const { ctx } = viewer();
    const d = smallDoc();
    d.snapshots = ["v.ruler", "v.letters", "v.dom"].map((t, k) => ({ ref: `s${k}`, ops: [hide(t)] }));
    await ctx.setDocument(d);
    const count = validations(ctx);
    const S = 3;
    expect(ctx.actions.run("track.hide", "v.site", "pointer")).toBe(true);
    await ctx.graph.update();
    expect(node(ctx.document.base.value, "v.site")?.params?.hidden).toBe(true);
    expect(count()).toBeLessThanOrEqual(1 + S);
    expect(await ctx.history.undo()).toBe(true);
    expect(count()).toBeLessThanOrEqual(2 * (1 + S));
    expect(ctx.snapshots.list.value.every((s) => s.valid)).toBe(true);
    expect(count()).toBeLessThanOrEqual(2 * (1 + S));
  });
});
