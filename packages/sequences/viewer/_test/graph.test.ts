/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The transform graph (seqspec:viewer#updates-and-status): counted runs per change, the axis identity, registry and
// skin inputs, unavailable and failing sources, aborted runs, views, cycles, and equality with compile + resolve.

import { createSeqViewSpec } from "@molstar/seq-core/builder";
import {
  compile,
  createRegistries,
  KernelRegistry,
  MvsRegistry,
  register,
  resolveView,
  type SeqViewSpecDoc,
  SimpleParamsSchema,
  type StreamSource,
} from "@molstar/seq-core/format";
import type { SeqSet } from "@molstar/seq-core/model";
import { tileClient, tilesProvider } from "@molstar/seq-core/stream/tiles";
import { tileServer } from "@molstar/seq-core/stream/tiles/server";
import { describe, expect, it } from "vitest";
import { createSeqViewer, DefaultNodeKinds, DefaultSkin, type NodeKind, type SeqViewerContext } from "../src/index.js";
import { manualLoaders, patch, refsRan, smallDoc, viewer } from "./graph-fixtures.js";

const act = (name: string) => ({ name });

const STAGES = /^(#prepare:|#theme:)|\.l$/;
/** Layer, prepare and theme nodes among the refs run. */
const perLayer = (ctx: SeqViewerContext) => refsRan(ctx).filter((r) => STAGES.test(r));
const set = (ref: string, params: Record<string, unknown>) => patch([{ op: "set", ref, params }]);

async function shown() {
  const v = viewer();
  await v.ctx.setDocument(smallDoc());
  return v;
}

describe("memo and early cutoff", () => {
  it("a set on one layer reruns that layer's stages and the view stage only", async () => {
    const { ctx, prepared } = await shown();
    const before = ctx.graph.runs;
    await ctx.applyPatch(set("v.dom.l", { color: { theme: "uniform", params: { color: "#ff0000" } } }));
    expect(refsRan(ctx)).toEqual(["#prepare:v.dom.l", "#theme:v.dom.l", "#view:v", "v.dom.l"]);
    expect(ctx.graph.runs - before).toBe(4);
    expect(prepared.get("v.dom.l")).toBe(2);
    expect(prepared.get("v.site.l")).toBe(1);
  });

  it("a set on a numbering reruns its parent's load node and the assemblies; layers only where the model rebuilt their data", async () => {
    const { ctx } = await shown();
    const set0 = ctx.graph.assembled.value?.set as SeqSet;
    await ctx.applyPatch(set("a.auth", { first: 201 }));
    const ran = refsRan(ctx);
    expect(ran.filter((r) => !STAGES.test(r))).toEqual([
      "#assemble",
      "#assemble.base",
      "#assemble.mapped",
      "#assemble.stream",
      "#view:v",
      "a",
      "v.main",
    ]);
    // `withChanges` rebuilds the annotation over the changed molecule (`dom`), so its layer reads new data; the
    // other layers, `site`'s on another molecule included, keep every result.
    expect(perLayer(ctx)).toEqual(["#prepare:v.dom.l", "#theme:v.dom.l", "v.dom.l"]);
    const set1 = ctx.graph.assembled.value?.set as SeqSet;
    expect(set1).not.toBe(set0);
    expect(set1.lineage).toBe(set0.lineage); // the set keeps its lineage
    expect(set1.annotations[0]).not.toBe(set0.annotations[0]);
    expect(set1.annotations[1]).toBe(set0.annotations[1]); // and its untouched rows
  });

  it("a filter view, a numbering switch, a track reorder or row sizing rerun no layer, prepare or theme node", async () => {
    const { ctx } = await shown();
    const view0 = ctx.graph.assembled.value?.view;
    for (const p of [
      set("v.main", { positions: ["dom", "site"] }),
      set("v.main", { numbering: "auth" }),
      patch([{ op: "move", ref: "v.site", parent: "v.main", after: "first" }]),
      set("v.main", { row_sizing: "natural" }),
    ]) {
      expect((await ctx.applyPatch(p)).ok).toBe(true);
      expect(refsRan(ctx)).toEqual(["#view:v", "v.main"]);
    }
    const view1 = ctx.graph.assembled.value?.view;
    expect(view1?.resolved.blocks.positions[0]).not.toBeNull(); // the filter is applied
    // The published view carries the new numbering, which a renderer reads live (RN1 draws it).
    expect(view1?.set.strings.get(view1?.resolved.blocks.numbering[0]?.name ?? -1)).toBe("auth");
    // Reused layers take their tracks from the section: the moved track's layer is the first track's now.
    const site = view1?.resolved.layers.ref.findIndex((r) => view1.set.strings.get(r) === "v.site.l") as number;
    expect(view1?.resolved.layers.track[site]).toBe(0);
    const states0 = view0?.prepared.map((p) => p.state) ?? [];
    expect(view1?.prepared.every((p) => states0.includes(p.state))).toBe(true); // the same objects, reordered
  });

  it("a skin change reruns theme nodes only, keeping every prepared state", async () => {
    const { ctx, prepared } = await shown();
    ctx.setSkin({ ...DefaultSkin, dark: true });
    await ctx.graph.update();
    expect(refsRan(ctx).filter((r) => !r.startsWith("#theme:"))).toEqual(["#view:v"]);
    expect(refsRan(ctx).filter((r) => r.startsWith("#theme:")).length).toBe(4);
    expect([...prepared.values()]).toEqual([1, 1, 1, 1]);
    expect(ctx.graph.assembled.value?.view?.prepared[1].themes[0]).toMatchObject({ name: "letter", dark: true });
  });

  it("ST2: a lazy alignment of a stream source reruns when an annotation names another of its rows", async () => {
    // Its named rows come from the document's annotations and mappings (`namedRows`), as a protocol source's do.
    const server = tileServer({ rowCount: 50, columnCount: 12 });
    const stream = (s: StreamSource) =>
      s.stream === "msa" ? tilesProvider(s, tileClient("https://x.test/msa", server.transport)) : null;
    const { ctx } = viewer({ loaders: { load: () => Promise.reject(new Error("no file")), stream } });
    const b = createSeqViewSpec();
    b.data.source({ ref: "msa.src", stream: "msa", data_kind: "alignment" });
    b.data.alignment({ ref: "msa", source_ref: "msa.src", row_count: 50, column_count: 12 } as never);
    b.data.annotation({ ref: "site", geometry: "ranges", target: "msa/3", value: { start: [0], stop: [2] } });
    b.view({ ref: "v", label: "V" }).block({ alignment: "msa" });
    await ctx.setDocument(b.getDocument());
    const rows = () => server.calls.filter((c) => c.startsWith("rows/")).sort();
    expect(rows()).toEqual(["rows/0", "rows/3"]);
    await ctx.applyPatch(set("site", { target: "msa/5" }));
    expect(refsRan(ctx)).toContain("msa");
    expect(rows()).toEqual(["rows/0", "rows/0", "rows/3", "rows/5"]);
  });

  it("a whole-document replacement that changed one node reruns that node only", async () => {
    const { ctx } = await shown();
    const next = smallDoc(); // every object new
    const dom = next.root.children?.find((c) => c.ref === "v")?.children?.[0].children?.[2].children?.[0];
    expect(dom?.ref).toBe("v.dom.l");
    (dom as { params: Record<string, unknown> }).params.color = { theme: "uniform", params: { color: "#00ff00" } };
    await ctx.setDocument(next);
    expect(refsRan(ctx).filter((r) => !r.startsWith("#"))).toEqual(["v.dom.l"]);
    await ctx.setDocument(smallDoc());
    await ctx.setDocument(smallDoc());
    expect(refsRan(ctx)).toEqual([]);
  });

  it("switching views reruns nothing already resolved; only the current view's stages run", async () => {
    const { ctx } = await shown();
    const v = ctx.graph.assembled.value?.view;
    expect(ctx.graph.status("w.main")?.state).toBe("pending");
    ctx.views.setCurrent("w");
    await ctx.graph.update();
    expect(refsRan(ctx)).toEqual(["#prepare:w.letters.l", "#theme:w.letters.l", "#view:w", "w.letters.l", "w.main"]);
    expect(ctx.graph.assembled.value?.view?.ref).toBe("w");
    ctx.views.setCurrent("v");
    await ctx.graph.update();
    expect(refsRan(ctx)).toEqual([]);
    expect(ctx.graph.assembled.value?.view).toBe(v);
  });
});

describe("registry inputs", () => {
  it("a registration that replaces nothing reruns nothing", async () => {
    const { ctx } = await shown();
    const undo = ctx.register({ actions: [act("lab:go")] });
    await ctx.graph.update();
    expect(refsRan(ctx)).toEqual([]);
    undo();
    await ctx.graph.update();
    expect(refsRan(ctx)).toEqual([]);
  });

  it("registering an extension reruns only the nodes whose entry changed", async () => {
    const { ctx } = viewer();
    const bars = {
      name: "lab:bars",
      options: SimpleParamsSchema({}),
      accepts: ["ranges"],
      layouts: ["tracks"],
      axes: ["blocks"],
      sizing: "row",
      fallback: "ranges",
      channels: [],
    } as const;
    ctx.register({ svs: { layerTypes: [bars] } });
    const d = smallDoc();
    const track = d.root.children?.find((c) => c.ref === "v")?.children?.[0].children?.[3];
    track?.children?.push({
      kind: "layer",
      ref: "v.bars.l",
      params: { type: "lab:bars", fallback: "ranges", data: "site" },
    });
    await ctx.setDocument(d);
    const fallback = ctx.graph.status("v.bars.l")?.issues.map((i) => i.code);
    expect(fallback).toEqual(["svs.render.fallback"]);
    const undo = ctx.register({ renderers: [{ name: "lab:bars", geometries: ["linear"], prepare: () => "bars" }] });
    await ctx.graph.update();
    expect(refsRan(ctx)).toEqual(["#prepare:v.bars.l", "#theme:v.bars.l", "#view:v", "v.bars.l"]);
    expect(ctx.graph.assembled.value?.view?.prepared.at(-1)?.state).toBe("bars");
    undo(); // it falls back again
    await ctx.graph.update();
    expect(refsRan(ctx)).toEqual(["#prepare:v.bars.l", "#theme:v.bars.l", "#view:v", "v.bars.l"]);
  });
});

/** A section over an unavailable sequence, one over a failing one, and one over a sequence that loads. */
function threeSources(): SeqViewSpecDoc {
  const b = createSeqViewSpec();
  b.data.sequence({ ref: "x", source: { value: ">x\nMKV\n", format: "genbank" } });
  b.data.sequence({ ref: "y", source: { uri: "https://x.test/bad.fa", format: "fasta" } });
  b.data.sequence({ ref: "z", value: "MKVLAT" });
  for (const s of ["x", "y", "z"]) {
    b.data.annotation({ ref: `on-${s}`, geometry: "ranges", target: s, value: { start: [0], stop: [2] } });
    const v = b.view({ ref: `v${s}`, label: s });
    const sec = v.section({ ref: `v.${s}` }).block({ sequence: s });
    sec.track({ ref: `t.${s}` }).ranges({ ref: `l.${s}`, data: `on-${s}` });
  }
  const d = b.getDocument();
  // One view of the three sections.
  const views = d.root.children?.filter((c) => c.kind === "view") ?? [];
  views[0].children = views.flatMap((v) => v.children ?? []);
  d.root.children = [d.root.children?.[0] as never, views[0]];
  return d;
}

describe("status", () => {
  it("a missing format leaves its readers unavailable and layers unplaced; a failing source omits its layers", async () => {
    const io = manualLoaders();
    const { ctx } = viewer({ loaders: io.loaders });
    const done = ctx.setDocument(threeSources());
    await new Promise((r) => setTimeout(r, 0));
    io.fail("https://x.test/bad.fa");
    await done;
    const st = (ref: string) => {
      const s = ctx.graph.status(ref);
      return [s?.state, s?.cause, ...(s?.issues.map((i) => `${i.code} ${i.severity}`) ?? [])];
    };
    expect(st("x.src")).toEqual(["unavailable", null, "svs.provider.unavailable warning"]);
    expect(st("x")).toEqual(["unavailable", "x.src"]);
    expect(st("on-x")).toEqual(["unavailable", "x"]);
    expect(st("l.x")).toEqual(["ok", null]); // no blocks: nothing more is reported, nothing cascades
    expect(st("y.src")).toEqual(["error", null, "svs.source.load error"]);
    expect(st("y")).toEqual(["unavailable", "y.src"]);
    expect(st("on-y")).toEqual(["unavailable", "y"]);
    expect(st("l.y")).toEqual(["unavailable", "on-y", "svs.render.omitted warning"]);
    for (const r of ["z", "on-z", "l.z", "v.z", "#assemble", "#view:vx"]) expect(st(r)).toEqual(["ok", null]);
    const view = ctx.graph.assembled.value?.view?.resolved;
    expect(view?.sections.count).toBe(3);
    expect(Array.from(view?.sections.blockCount ?? [])).toEqual([0, 0, 1]);
    expect(
      Array.from(view?.layers.ref ?? []).map((r) => view && ctx.graph.assembled.value?.set.strings.get(r)),
    ).toEqual(["l.x", "l.z"]);
    expect(view?.placements.count).toBe(1);
  });

  it("an aborted slow source never publishes", async () => {
    const io = manualLoaders();
    const { ctx } = viewer({ loaders: io.loaders });
    const sets: SeqSet[] = [];
    ctx.graph.assembled.subscribe((a) => a && sets.push(a.set));
    const doc = (uri: string) => {
      const b = createSeqViewSpec();
      b.data.sequence({ ref: "s", source: { uri, format: "fasta" } });
      b.view({ ref: "v", label: "V" }).track().letters();
      return b.getDocument();
    };
    const first = ctx.setDocument(doc("https://x.test/slow.fa"));
    await new Promise((r) => setTimeout(r, 0));
    const second = ctx.setDocument(doc("https://x.test/fast.fa"));
    await new Promise((r) => setTimeout(r, 0));
    expect(io.requests).toEqual(["https://x.test/slow.fa", "https://x.test/fast.fa"]);
    io.resolve("https://x.test/fast.fa", ">s\nMKV\n");
    await Promise.all([first, second]);
    io.resolve("https://x.test/slow.fa", ">s\nAAAAAAAAAA\n"); // after its run was superseded: dropped
    await new Promise((r) => setTimeout(r, 0));
    await ctx.graph.update();
    expect(sets.length).toBe(1);
    expect(sets.map((s) => s.molecules[0]?.length)).toEqual([3]);
    expect(ctx.graph.assembled.value?.set.molecules[0]?.length).toBe(3);
  });

  it("dispose aborts running graph work", async () => {
    let aborted = false;
    const { ctx } = viewer({
      loaders: {
        load: (_r, signal) =>
          new Promise((_, reject) =>
            signal.addEventListener("abort", () => {
              aborted = true;
              reject(signal.reason);
            }),
          ),
      },
    });
    void ctx.setDocument(threeSources());
    await new Promise((r) => setTimeout(r, 0));
    ctx.dispose();
    expect(aborted).toBe(true);
  });

  it("a reference cycle is svs.ref.cycle on every node in it; the rest runs", async () => {
    // A host's node kind that reads another annotation (the format's validation keeps cycles out of documents).
    const base = DefaultNodeKinds.find((k) => k.name === "annotation") as NodeKind;
    const cyclic: NodeKind = {
      ...base,
      inputs(ref, n, doc, out, svs) {
        base.inputs(ref, n, doc, out, svs);
        out.status(ref === "dom" ? "site" : "dom");
      },
    };
    const kinds = [...DefaultNodeKinds.filter((k) => k !== base), cyclic];
    const ctx = createSeqViewer({
      registry: [{ svs: [KernelRegistry, MvsRegistry], nodeKinds: kinds }],
      behaviors: [],
    });
    await ctx.setDocument(smallDoc());
    for (const r of ["dom", "site"]) {
      expect(ctx.graph.status(r)?.state).toBe("error");
      expect(ctx.graph.status(r)?.issues.map((i) => i.code)).toEqual(["svs.ref.cycle"]);
    }
    expect(ctx.graph.status("a")?.state).toBe("ok");
    expect(ctx.graph.assembled.value?.view?.resolved.blocks.count).toBe(2);
  });
});

describe("equal to the format's drivers", () => {
  it("setDocument gives what compile + resolveView give for the same document", async () => {
    const { ctx } = await shown();
    const reg = createRegistries();
    register(reg, [KernelRegistry, MvsRegistry]);
    const check = async (doc: SeqViewSpecDoc) => {
      const a = ctx.graph.assembled.value;
      const fresh = await compile(doc, reg, { load: () => Promise.reject(new Error("none")) });
      const supported = ctx.graph.supported(0);
      const graphSet = a?.set as SeqSet;
      expect(a?.view?.resolved).toEqual(resolveView(doc, graphSet, "v", reg, supported).view);
      const keys = (s: SeqSet) => [s.molecules, s.annotations].map((l) => l.map((x) => x && s.strings.get(x.key)));
      expect(keys(graphSet)).toEqual(keys(fresh.set));
      expect(graphSet.molecules.map((m) => m?.numberings.count)).toEqual(
        fresh.set.molecules.map((m) => m?.numberings.count),
      );
      return fresh;
    };
    // The first document: the same set as compile's, value for value, and the same view.
    const fresh = await check(smallDoc());
    expect(ctx.graph.assembled.value?.view?.resolved).toEqual(
      resolveView(smallDoc(), fresh.set, "v", reg, ctx.graph.supported(0)).view,
    );
    // After changes, the graph's set and view are what the stages give for the shown document.
    await ctx.applyPatch(set("a.auth", { first: 7 }));
    await ctx.applyPatch(set("v.main", { positions: "site", numbering: "auth" }));
    await check(ctx.document.shown.value as SeqViewSpecDoc);
  });
});
