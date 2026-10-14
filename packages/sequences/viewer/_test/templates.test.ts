/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// ST1b-b: track templates and nested sources in the viewer (seqspec:format#placement, seqspec:sources#dynamic-nodes;
// seqspec:viewer#updates-and-status, seqspec:sources#the-viewers-side). A template's instances come from the dynamic
// nodes it selects, at its place in key order; an arrival or eviction under its source reruns its `#template` and the
// view, nothing else, keeping every other instance's prepared state; another source's reruns only `#assemble.stream`.
// Nested sources are asked for only when a template (or a reader) needs them, to `max_depth`. Instances hold nothing (a
// region's node leaves with it); a highlight holds what it names; and `ctx.source(ref)` of an annotation with a `uri`
// works in the graph.

import { createSeqViewSpec } from "@molstar/seq-core/builder";
import {
  RequiredField,
  ref,
  type SeqViewSpecDoc,
  SimpleParamsSchema,
  type StageContext,
} from "@molstar/seq-core/format";
import { type AnnotationInput, type KeyTarget, OrderedSet, type SeqSet } from "@molstar/seq-core/model";
import { describe, expect, it } from "vitest";
import { fixtureLoaders } from "../../core/_test/format/load-fixtures.ts";
import { GRCH38 } from "../../core/_test/format/stream-provider.ts";
import { colorSchemePickers } from "../src/color-schemes.ts";
import { opened } from "./layout-fixtures.ts";
import {
  CacheBytes,
  dynamicDoc,
  GENES,
  GENOME,
  MaxInFlight,
  openGenome,
  PEAK_10MB,
  type StreamedGenome,
  withPeaks,
} from "./stream-fixtures.ts";
import { nested, openTemplates, range } from "./template-fixtures.ts";

const flush = () => new Promise((ok) => setTimeout(ok, 0));
const codes = (t: { ctx: { issues: { current: { value: readonly { code: string }[] } } } }) =>
  t.ctx.issues.current.value.map((i) => i.code);

describe("a template's instances", () => {
  it("one per node its source holds, at the template's place, in key order", async () => {
    const t = await openTemplates({ genes: [range("b", 10, "Gene B"), range("a", 1)] });
    expect(t.tracks()).toEqual(["seq", "tpl/genes/a", "tpl/genes/b", "ruler"]);
    expect(codes(t)).not.toContain("svs.layer.unplaced");
    t.ctx.dispose();
  });

  it("an arrival reruns its template and the view only, every other instance kept; another source's reruns none", async () => {
    const t = await openTemplates({ genes: [range("a", 1), range("b", 10)] });
    const stateOf = (node: string) => {
      const v = t.ctx.graph.assembled.value?.view;
      for (let k = 0; v && k < v.resolved.layers.count; k++)
        if (v.resolved.layers.node[k] === node) return v.prepared[k].state;
      return undefined;
    };
    const a = stateOf("genes/a");
    expect(a).toBeDefined();
    t.manager.nodes.add("genes", [range("c", 20)]);
    await t.frame();
    expect(t.instances()).toEqual(["genes/a", "genes/b", "genes/c"]);
    expect([...t.ctx.graph.ran].sort()).toEqual(["#assemble.stream", "#template:tpl", "#view:v"]);
    expect(stateOf("genes/a")).toBe(a);
    t.manager.nodes.add("notes", [range("x", 30)]);
    await t.frame();
    expect([...t.ctx.graph.ran]).toEqual(["#assemble.stream"]);
    // Its region gone, the node goes: an instance holds nothing.
    t.manager.nodes.release(["genes/c"]);
    await t.frame();
    expect(t.instances()).toEqual(["genes/a", "genes/b"]);
    expect(stateOf("genes/a")).toBe(a);
    t.ctx.dispose();
  });
});

describe("a color scheme picked on an instance", () => {
  it("recolors that instance, its template rerun, the others kept", async () => {
    const t = await openTemplates({ genes: [range("a", 1), range("b", 10)] });
    const at = (node: string) => {
      const v = t.ctx.graph.assembled.value?.view;
      for (let k = 0; v && k < v.resolved.layers.count; k++)
        if (v.resolved.layers.node[k] === node) return { k, track: v.resolved.layers.track[k], p: v.prepared[k] };
      throw new Error(`no instance ${node}`);
    };
    const [a, b] = [at("genes/a"), at("genes/b")];
    const picker = colorSchemePickers(t.ctx, a.track)[0];
    const scheme = picker?.options.find((o) => o.value !== null)?.value;
    if (!picker || !scheme) throw new Error("no scheme to pick");
    t.ctx.ui.setColorScheme(picker.layer, scheme);
    await t.ctx.graph.update();
    expect(t.ctx.graph.ran.has("#template:tpl")).toBe(true);
    expect(at("genes/a").p.themes).not.toBe(a.p.themes);
    expect(at("genes/b").p.themes).toBe(b.p.themes);
    t.ctx.dispose();
  });
});

describe("nested sources", () => {
  const manifests = {
    genes: [range("a", 1), nested("g1/tx", "tx-g1")],
    "tx-g1": [range("t1", 5), nested("t1/dom", "dom-t1")],
    "dom-t1": [range("d1", 9)],
  };
  const open = (depth: number | null, max_depth?: number) =>
    openTemplates(manifests, { for_each: { kind: "annotation", source: "genes", depth }, max_depth });

  it("asked for only when a template reaches their nodes; their manifests' nodes then instances", async () => {
    const own = await open(0);
    expect([own.instances(), own.asked.includes("tx-g1")]).toEqual([["genes/a"], false]);
    own.ctx.dispose();
    const all = await open(null);
    expect(all.instances()).toEqual(["genes/a", "genes/g1/tx/t1", "genes/g1/tx/t1/dom/d1"]);
    expect(all.asked.filter((s) => s !== "genes" && s !== "notes")).toEqual(["tx-g1", "dom-t1"]);
    all.ctx.dispose();
    const one = await open(1);
    expect([one.instances(), one.asked.includes("dom-t1")]).toEqual([["genes/g1/tx/t1"], false]);
    one.ctx.dispose();
  });

  it("their nodes stay through a later arrival, each manifest asked for once", async () => {
    const t = await open(null);
    t.manager.nodes.add("genes", [range("b", 3)]);
    await t.frame();
    expect(t.instances()).toEqual(["genes/a", "genes/b", "genes/g1/tx/t1", "genes/g1/tx/t1/dom/d1"]);
    expect(t.asked.filter((s) => s === "tx-g1" || s === "dom-t1")).toEqual(["tx-g1", "dom-t1"]);
    t.ctx.dispose();
  });

  it("a source deeper than max_depth is svs.stream.depth, never asked for", async () => {
    const t = await open(null, 2);
    expect(t.instances()).toEqual(["genes/a", "genes/g1/tx/t1"]);
    expect(t.asked).not.toContain("dom-t1");
    expect(codes(t)).toContain("svs.stream.depth");
    t.ctx.dispose();
  });
});

// --- over the streamed genome: regions carry the nodes -------------------------------------------------------------

const t0: KeyTarget = { kind: -1, id: -1, row: -1 };
const idOf = (set: SeqSet | undefined, key: string) => (set?.byName(key, t0) ? t0.id : -1);
const setOf = (t: StreamedGenome) => t.ctx.graph.assembled.value?.set;
async function visit(t: StreamedGenome, b: number, start: number, stop: number) {
  t.vp(GENOME).zoomTo(b, start, stop);
  t.layout(GENOME);
  t.manager.pass();
  for (let k = 0; k < 4; k++) {
    await flush();
    t.manager.pass();
  }
  t.scheduler.advance(0);
  await t.ctx.graph.update();
}
/** The dynamic-nodes genome with a template over every peak its regions bring. */
function peaksDoc(): SeqViewSpecDoc {
  const doc = dynamicDoc();
  const section = doc.root.children?.find((c) => c.kind === "view")?.children?.find((c) => c.kind === "section");
  section?.children?.push({
    kind: "track",
    ref: "peaks",
    params: { for_each: { kind: "annotation", source: "phylop_src" } },
    children: [{ kind: "layer", ref: "peaks.l", params: { type: "ranges", data: "$node" } }],
  });
  return doc;
}
const instancesOf = (t: StreamedGenome) => {
  const v = t.ctx.graph.assembled.value?.view;
  const out: string[] = [];
  for (let k = 0; v && k < v.resolved.layers.count; k++) {
    const node = v.resolved.layers.node[k];
    if (node) out.push(node);
  }
  return out;
};
const OTHER = "phylop_src/peak-chr1-2442";

describe("what holds a template's nodes", () => {
  const open = () =>
    openGenome({ chromosomes: GRCH38.slice(0, 3), doc: peaksDoc(), nodes: [GENES], answer: withPeaks() });

  it("nothing: an instance's node leaves with its region; one a document track shows stays", async () => {
    const t = await open();
    t.ctx.config.set(MaxInFlight, 100);
    await visit(t, 0, 10_000_000, 10_010_000);
    expect(instancesOf(t)).toEqual(expect.arrayContaining([PEAK_10MB, OTHER, "phylop_src/genes"]));
    t.ctx.config.set(CacheBytes, 150_000);
    for (const at of [60, 110, 160]) await visit(t, 1, at * 1e6, at * 1e6 + 10_000);
    expect(instancesOf(t)).toContain(PEAK_10MB);
    expect(instancesOf(t)).not.toContain(OTHER);
    expect(idOf(setOf(t), OTHER)).toBe(-1);
    t.ctx.dispose();
  });

  it("a highlight holds the node it names until it is cleared", async () => {
    const t = await open();
    t.ctx.config.set(MaxInFlight, 100);
    await visit(t, 0, 10_000_000, 10_010_000);
    const set = setOf(t) as SeqSet;
    const id = idOf(set, OTHER);
    const items = { kind: "item-loci" as const, set, elements: [{ annotation: id, items: OrderedSet.ofSingleton(0) }] };
    t.ctx.highlight.set("probe", items as never);
    t.ctx.config.set(CacheBytes, 150_000);
    for (const at of [60, 110]) await visit(t, 1, at * 1e6, at * 1e6 + 10_000);
    expect(idOf(setOf(t), OTHER)).toBe(id);
    t.ctx.highlight.clear("probe");
    t.scheduler.advance(0);
    await t.ctx.graph.update();
    expect(idOf(setOf(t), OTHER)).toBe(-1);
    t.ctx.dispose();
  });
});

describe("ctx.source(ref) of an annotation with a uri, in the graph", () => {
  it("a provider reads the file the annotation's uri names, by the annotation's ref", async () => {
    const b = createSeqViewSpec({ title: "t" });
    const p = b.data.sequence({ ref: "p", value: "MKVLAAGGHH" });
    b.data.annotation({ ref: "cats", geometry: "ranges", target: "p", uri: "https://x.org/cats.tsv", format: "tsv" });
    b.data.computed({ ref: "echo", provider: "test:rows_of", options: { of: "cats" } });
    const v = b.view({ ref: "v", label: "v" });
    v.block({ sequence: p });
    v.track({ label: "Echo" }).ranges({ data: "echo" });
    const echo = {
      name: "test:rows_of",
      options: SimpleParamsSchema({ of: RequiredField(ref("annotation"), "") }),
      geometry: "ranges" as const,
      async compute(ctx: StageContext, o: unknown) {
        const t = ctx.source((o as { of: string }).of);
        if (t?.kind !== "table") return null;
        const n = t.length;
        const cols = t.columns as Record<string, ArrayLike<unknown>>;
        const columns = { target: Array(n).fill("p"), start: cols.start, stop: cols.stop };
        return { key: "", geometry: "ranges", columns } as AnnotationInput;
      },
    };
    const loaders = fixtureLoaders({ "https://x.org/cats.tsv": "start\tstop\tcategory\n1\t4\thelix\n" });
    const t = await opened(b.getDocument(), { loaders, registry: [{ svs: { providers: [echo] } } as never] });
    expect(idOf(t.ctx.graph.assembled.value?.set, "echo")).toBeGreaterThanOrEqual(0);
    t.ctx.dispose();
  });
});
