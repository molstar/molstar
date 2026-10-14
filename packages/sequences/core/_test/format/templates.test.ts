/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// ST1b-b: nested sources and track templates in the format (seqspec:format#placement, seqspec:sources#dynamic-nodes). A
// template's `"$node"` is valid only in its layers; `loadStreamNode` reads a nested source (its manifest asked for, its
// parent named) within `max_depth`, and a node reading it by `source_ref`; `templateNodes` selects by kind, source and
// depth in key order; `resolveView` puts the instances at the template's place; and `ctx.source(ref)` finds the file an
// annotation's `uri` names when a driver loads node by node.

import { createSeqViewSpec } from "@molstar/seq-core/builder";
import { type Issue, type SeqSet, SeqSets } from "@molstar/seq-core/model";
import { describe, expect, it } from "vitest";
import {
  compile,
  type DocIndex,
  type IndexedNode,
  indexDoc,
  type Loaders,
  loadNode,
  loadStreamNode,
  type Manifest,
  MapOver,
  resolveView,
  type SourceResult,
  type StageContext,
  type StreamedData,
  StreamProp,
  sourceParse,
  supportedAll,
  templateNodes,
  validate,
} from "../../src/format/index.ts";
import { fixtureLoaders, registries } from "./load-fixtures.ts";
import { scripted } from "./stream-provider.ts";

const reg = registries();
const codes = (issues: readonly Issue[]) => issues.map((i) => i.code);
const P = "MKVLAAGGHHKKLLMMNNPPQQRRSSTTVVWWYY";
type ForEach = { kind: "annotation"; source: string; depth?: number | null; reaching?: boolean };

/** A protein, a streamed source of genes, and a section of its letters with a template (`for_each`) of ranges. */
function doc(o: { for_each?: ForEach; data?: string; max_depth?: number } = {}) {
  const b = createSeqViewSpec({ title: "t" });
  const p = b.data.sequence({ ref: "p", value: P });
  b.data.sequence({ ref: "q", value: "MKV" });
  b.data.source({ ref: "genes", stream: "genes", data_kind: "items", max_depth: o.max_depth });
  const v = b.view({ ref: "v", label: "v" });
  v.block({ sequence: p });
  v.track({ ref: "seq", label: "Seq" }).letters();
  v.track({ ref: "tpl", for_each: o.for_each ?? { kind: "annotation", source: "genes" } }).ranges({ data: "$node" });
  v.track({ ref: "last", label: "Last" }).ruler();
  if (o.data) v.track({ ref: "other" }).ranges({ data: o.data });
  return b.getDocument();
}
const ann = (key: string, target = "p", label?: string) => ({
  key,
  geometry: "ranges" as const,
  target,
  label,
  columns: { start: [1], stop: [5] },
});

describe("a template's node", () => {
  it('"$node" is valid in a template\'s layers only; its source is a streamed source, or a nested one by ref', () => {
    expect(codes(validate(doc(), reg))).toEqual([]);
    expect(codes(validate(doc({ for_each: { kind: "annotation", source: "genes/g1/tx", depth: 0 } }), reg))).toEqual(
      [],
    );
    expect(codes(validate(doc({ data: "$node" }), reg))).toEqual(["svs.ref.syntax"]);
    expect(codes(validate(doc({ for_each: { kind: "annotation", source: "p" } }), reg))).toEqual([
      "svs.ref.wrong_kind",
    ]);
  });
});

describe("templateNodes", () => {
  const set = SeqSets.build({
    molecules: [
      { key: "p", residues: P, ownUnit: true },
      { key: "genes/s1", residues: "MK", ownUnit: true },
    ],
    annotations: [
      ann("genes/c"),
      ann("genes/a"),
      ann("genes/g1/tx/t1"),
      ann("genes/g1/tx/t2/dom/d1"),
      ann("genes/b"),
      ann("other/x"),
    ],
  }).set;
  const nested = new Set(["genes/g1/tx", "genes/g1/tx/t2/dom"]);
  const pick = (fe: ForEach | { kind: "sequence"; source: string }) =>
    templateNodes("tpl", indexDoc(doc({ for_each: fe as ForEach }), reg), set, nested);

  it("the set's nodes of its kind under its source, at its depth of nested sources, in key order", () => {
    expect(pick({ kind: "annotation", source: "genes", depth: 0 })).toEqual(["genes/a", "genes/b", "genes/c"]);
    expect(pick({ kind: "annotation", source: "genes", depth: 1 })).toEqual(["genes/g1/tx/t1"]);
    expect(pick({ kind: "annotation", source: "genes", depth: 2 })).toEqual(["genes/g1/tx/t2/dom/d1"]);
    expect(pick({ kind: "annotation", source: "genes" })).toHaveLength(5);
    expect(pick({ kind: "annotation", source: "genes/g1/tx", depth: 0 })).toEqual(["genes/g1/tx/t1"]);
    expect(pick({ kind: "sequence", source: "genes" })).toEqual(["genes/s1"]);
  });
});

describe("resolveView with a template", () => {
  async function view(fe: ForEach, extra: ReturnType<typeof ann>[]) {
    const d = doc({ for_each: fe });
    const base = await compile(d, reg, { load: () => Promise.reject(new Error("no files")), stream: () => null });
    const set = SeqSets.withChanges(base.set, { upsert: { annotations: extra } }).set;
    const r = resolveView(indexDoc(d, reg), set, "v", reg, supportedAll(reg));
    const s = (id: number) => (id < 0 ? null : set.strings.get(id));
    const tracks = [...Array(r.view.tracks.count).keys()].map((k) => [
      s(r.view.tracks.ref[k]),
      s(r.view.tracks.label[k]),
    ]);
    const layers = [...Array(r.view.layers.count).keys()].map((k) => [s(r.view.layers.ref[k]), r.view.layers.node[k]]);
    return { tracks, layers, issues: r.issues, set };
  }

  it("an instance per node at the template's place, in key order, labeled by its node; a node placing nothing left out", async () => {
    const r = await view({ kind: "annotation", source: "genes" }, [
      ann("genes/b", "p", "Gene B"),
      ann("genes/a", "p"),
      ann("genes/z", "q"),
    ]);
    expect(r.tracks).toEqual([
      ["seq", "Seq"],
      ["tpl/genes/a", "a"],
      ["tpl/genes/b", "Gene B"],
      ["last", "Last"],
    ]);
    expect(r.layers.filter(([, node]) => node !== null)).toEqual([
      ["v.ranges/genes/a", "genes/a"],
      ["v.ranges/genes/b", "genes/b"],
    ]);
    expect(codes(r.issues)).not.toContain("svs.layer.unplaced");
  });

  it("reaching: false keeps every node; no node shows no track", async () => {
    const all = await view({ kind: "annotation", source: "genes", reaching: false }, [ann("genes/z", "q")]);
    expect(all.tracks.map(([ref]) => ref)).toEqual(["seq", "tpl/genes/z", "last"]);
    const none = await view({ kind: "annotation", source: "genes" }, []);
    expect(none.tracks.map(([ref]) => ref)).toEqual(["seq", "last"]);
  });
});

describe("nested sources", () => {
  const manifest = (nodes: Manifest["nodes"] = []): Manifest => ({
    systems: [{ key: "p", extent: [0, P.length], numbering: null, label: null }],
    levels: [],
    nodes,
  });
  function setup(maxDepth?: number) {
    const tx = scripted(manifest(), () => new Error("no data"));
    const dom = scripted(manifest(), () => new Error("no data"));
    const loaders: Loaders = {
      load: () => Promise.reject(new Error("no files")),
      stream: (s) => (s.stream === "tx-g1" ? tx.provider : s.stream === "dom-t1" ? dom.provider : null),
    };
    const issues: Issue[] = [];
    const index = indexDoc(doc({ max_depth: maxDepth }), reg);
    const nested = new Map<string, IndexedNode>();
    const results = new Map<string, SourceResult>();
    const ctx: StageContext = {
      signal: new AbortController().signal,
      registries: reg,
      set: null,
      sources: new MapOver(new Map(), results),
      unavailable: new Set(),
      source: () => null,
      report: (i) => issues.push(i),
    };
    const at: DocIndex = { nodes: new MapOver(index.nodes, nested), order: index.order };
    const known = (r: string) => r === "p" || nested.has(r);
    /** Reads a nested source node of `source`, holding it (and its result) as the viewer's assembly does. */
    const source = async (parent: string, key: string, stream: string) => {
      const ref = `${parent}/${key}`;
      const params = { stream, data_kind: "items" };
      const r = await loadStreamNode(parent, { kind: "source", key, params }, known, ctx, loaders, at);
      nested.set(ref, { node: { kind: "source", ref, params }, ref, kind: "source", params, parent, path: "" });
      if (r?.source) results.set(ref, r.source);
      return r;
    };
    return { tx, dom, issues, ctx, at, known, loaders, source };
  }

  it("is read as a streamed source: its manifest asked for, its parent named", async () => {
    const t = setup();
    const r = await t.source("genes", "g1/tx", "tx-g1");
    expect(t.issues).toEqual([]);
    expect(r?.source?.kind).toBe("streamed");
    const s = r?.source as Extract<SourceResult, { kind: "streamed" }>;
    expect([s.source.ref, s.source.parent, s.source.data_kind]).toEqual(["genes/g1/tx", "genes", "items"]);
    // A node reading it by `source_ref` is a streamed annotation of it.
    const a = await loadStreamNode(
      "genes",
      { kind: "annotation", key: "g1/items", params: { geometry: "ranges", source_ref: "genes/g1/tx", target: "p" } },
      t.known,
      t.ctx,
      t.loaders,
      t.at,
    );
    const d = a?.annotations?.[0]?.props?.[StreamProp] as StreamedData | undefined;
    expect([a?.annotations?.[0]?.key, d?.source.ref, d?.source.parent]).toEqual([
      "genes/g1/items",
      "genes/g1/tx",
      "genes",
    ]);
  });

  it("deeper than its chain's max_depth is svs.stream.depth, and never requested", async () => {
    const t = setup(2);
    await t.source("genes", "g1/tx", "tx-g1");
    const deep = await t.source("genes/g1/tx", "t1/dom", "dom-t1");
    expect(deep).toBeNull();
    expect(t.issues.map((i) => [i.code, i.severity, i.path])).toEqual([
      ["svs.stream.depth", "warning", "/nodes/t1/dom"],
    ]);
    expect(t.issues[0].message).toContain("genes/g1/tx");
    expect(t.dom.requests).toEqual([]);
    // At the default (4), the same source is read.
    const ok = setup();
    await ok.source("genes", "g1/tx", "tx-g1");
    expect((await ok.source("genes/g1/tx", "t1/dom", "dom-t1"))?.source?.kind).toBe("streamed");
  });
});

describe("ctx.source(ref) of an annotation with a uri", () => {
  it("a driver loading node by node gets the file back as the annotation's source, found by its ref", async () => {
    const b = createSeqViewSpec({ title: "t" });
    b.data.sequence({ ref: "p", value: P });
    b.data.annotation({ ref: "cats", geometry: "ranges", target: "p", uri: "https://x.org/cats.tsv", format: "tsv" });
    const d = b.getDocument();
    const index = indexDoc(d, reg);
    const loaders = fixtureLoaders({ "https://x.org/cats.tsv": "start\tstop\tcategory\n1\t4\thelix\n" });
    const set: SeqSet = (await compile(d, reg, loaders)).set;
    const sources = new Map<string, SourceResult>();
    const ctx: StageContext = {
      signal: new AbortController().signal,
      registries: reg,
      set,
      sources,
      unavailable: new Set(),
      source: (ref) => sourceParse(index, sources, ref),
      report: () => {},
    };
    const r = await loadNode(index.nodes.get("cats")?.node as never, index, ctx, loaders);
    expect(r.annotations?.[0]?.key).toBe("cats");
    expect(r.source?.kind).toBe("loaded");
    sources.set("cats", r.source as SourceResult);
    expect(ctx.source("cats")).toMatchObject({ kind: "table", length: 1 });
  });
});
