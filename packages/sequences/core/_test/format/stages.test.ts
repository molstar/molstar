/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The stage functions called one node at a time, as the viewer's graph calls them (seqspec:viewer#updates-and-status):
// explicit sources and unavailability, per-node normalization, the scroll-group check in `composeView`, layer reuse.

import { createSeqViewSpec } from "@molstar/seq-core/builder";
import type { Issue, SeqSet } from "@molstar/seq-core/model";
import { describe, expect, it } from "vitest";
import {
  applyPatch,
  assembleSet,
  compile,
  composeView,
  type DocIndex,
  type IndexedNode,
  indexDoc,
  int,
  type LayerTypeEntry,
  type LoadResult,
  loadNode,
  type NormalizedParams,
  normalize,
  OptionalField,
  type Patch,
  type ResolvedLayer,
  register,
  resolveLayer,
  resolveSection,
  resolveView,
  type SeqViewSpecDoc,
  SimpleParamsSchema,
  type SourceResult,
  type StageContext,
  sourceParse,
  supportedAll,
} from "../../src/format/index.ts";
import { checkUnitKinds, remapIssues } from "../../src/format/load/compile.ts";
import { case10 } from "./cases.ts";
import { fixtureLoaders, registries } from "./load-fixtures.ts";

const reg = registries();
const all = supportedAll(reg);
const patched = (doc: SeqViewSpecDoc, ops: Patch["ops"]) => {
  const r = applyPatch(doc, { format: "seq_view_spec_patch", ops }, reg);
  expect(r.issues).toEqual([]);
  return r.doc;
};
const keyed = (issues: Issue[]) => issues.map((i) => `${i.code} ${i.path}`).sort();
const has = (set: SeqSet, key: string) => set.strings.find(key) >= 0 && set.byName(key, { kind: -1, id: -1, row: -1 });

/** A source, a node reading it, an unavailable source (a format no package registered) and its readers. */
function doc() {
  const b = createSeqViewSpec();
  const pasted = b.data.source({ ref: "pasted", value: ">a\nMKV\n>b\nMKVLA\n", format: "fasta" });
  b.data.sequence({ ref: "a", source_ref: pasted, record: "b" });
  const ext = b.data.source({ ref: "ext", value: ">x\nMKV\n", format: "genbank" });
  b.data.sequence({ ref: "x", source_ref: ext });
  b.data.annotation({ ref: "on-x", geometry: "ranges", target: "x", value: { start: [0], stop: [2] } });
  b.data.annotation({ ref: "on-a", geometry: "ranges", target: "a", value: { start: [0], stop: [2] } });
  b.view({ ref: "v", label: "V" }).track({ label: "L" }).letters();
  return b.getDocument();
}

describe("the load stages, one node at a time", () => {
  it("give the set and issues compile gives, sources and unavailability passed in explicitly", async () => {
    const d = doc();
    const index = indexDoc(d, reg);
    const loaders = fixtureLoaders();
    const sources = new Map<string, SourceResult>();
    const unavailable = new Set<string>();
    const status = new Map<string, { unavailable: boolean; issues: Issue[] }>();
    let set: SeqSet | null = null;
    // One graph node: its own context and issues; its result's source and unavailability go to the caller's tables.
    const run = async (n: IndexedNode, ctxSet: SeqSet | null) => {
      const issues: Issue[] = [];
      const ctx: StageContext = {
        signal: new AbortController().signal,
        registries: reg,
        set: ctxSet,
        sources,
        unavailable,
        source: (ref) => sourceParse(index, sources, ref),
        report: (i) => issues.push(i),
      };
      const r: LoadResult = await loadNode(n.node, index, ctx, loaders);
      if (r.source) sources.set(n.ref, r.source);
      if (r.unavailable) unavailable.add(n.ref);
      status.set(n.ref, { unavailable: r.unavailable === true, issues });
      return r;
    };
    const of = (kind: string) => index.order.filter((n) => n.kind === kind);
    for (const n of of("source")) await run(n, null);
    const molecules = [];
    for (const n of of("sequence")) molecules.push(...((await run(n, null)).molecules ?? []));
    set = assembleSet(null, { molecules }, []).set;
    const annotations = [];
    for (const n of of("annotation")) annotations.push(...((await run(n, set)).annotations ?? []));
    set = assembleSet(set, { annotations }, []).set;

    expect(sources.get("pasted")?.kind).toBe("loaded");
    expect(sources.get("ext")?.kind).toBe("unavailable");
    const ext = status.get("ext");
    expect(ext?.unavailable).toBe(true);
    expect(ext?.issues.map((i) => [i.code, i.severity, i.path])).toEqual([
      ["svs.provider.unavailable", "warning", index.nodes.get("ext")?.path],
    ]);
    for (const ref of ["x", "on-x"]) expect(status.get(ref)).toEqual({ unavailable: true, issues: [] });
    for (const ref of ["a", "on-a"]) expect(status.get(ref)).toEqual({ unavailable: false, issues: [] });

    const compiled = await compile(d, reg, fixtureLoaders());
    const perNode = [...status.values()].flatMap((s) => s.issues);
    expect(keyed(perNode)).toEqual(keyed(compiled.issues));
    for (const key of ["a", "on-a", "x", "on-x"]) expect(has(set, key)).toBe(has(compiled.set, key));
    expect([has(set, "a"), has(set, "on-a"), has(set, "x"), has(set, "on-x")]).toEqual([true, true, false, false]);
  });
});

describe("compile's checks, for a driver that assembles itself", () => {
  it("checkUnitKinds leaves out a unit naming a ligand as a sequence; remapIssues gives model paths as the document's", () => {
    const b = createSeqViewSpec();
    b.data.sequence({ ref: "s", value: "MKV" });
    b.data.annotation({ ref: "a", geometry: "ranges", target: "s", value: { start: [0], stop: [9] } });
    const d = b.getDocument();
    d.root.children?.[0]?.children?.push({ kind: "unit", ref: "u", params: { sequence: "lig" } });
    const index = indexDoc(d, reg);
    const issues: Issue[] = [];
    const ctx: StageContext = {
      signal: new AbortController().signal,
      registries: reg,
      set: null,
      sources: new Map(),
      unavailable: new Set(),
      source: () => null,
      report: (i) => issues.push(i),
    };
    const input = {
      molecules: [
        { key: "s", residues: "MKV" },
        { key: "lig", kind: "ligand" as const, elements: ["C"] },
      ],
      units: [
        { key: "u", molecule: "lig" },
        { key: "s/unit", molecule: "s" },
      ],
    };
    checkUnitKinds(input, index, ctx);
    expect(input.units.map((u) => u.key)).toEqual(["s/unit"]);
    expect(issues.map((i) => [i.code, i.path])).toEqual([
      ["svs.ref.wrong_kind", `${index.nodes.get("u")?.path}/params/sequence`],
    ]);
    const upsert = {
      annotations: [{ key: "a", target: "s", geometry: "ranges" as const, columns: { start: [0], stop: [9] } }],
    };
    const r = assembleSet(assembleSet(null, { molecules: [input.molecules[0]] }, []).set, upsert, []);
    expect(r.issues.length).toBeGreaterThan(0);
    const node = index.nodes.get("a")?.path ?? "";
    expect(remapIssues(r.issues, upsert, index).every((i) => i.path.startsWith(`${node}/params/value`))).toBe(true);
  });
});

describe("indexDoc normalizes per node", () => {
  it("keeps an unchanged node's normalized params object across documents, with the caller's cache", () => {
    const cache = new WeakMap<object, NormalizedParams>();
    const d = case10();
    const a = indexDoc(d, reg, cache);
    const changed = patched(d, [{ op: "set", ref: "both.ann", params: { height: "box" } }]);
    const b = indexDoc(changed, reg, cache);
    const params = (x: DocIndex, ref: string) => x.nodes.get(ref)?.params;
    expect(params(b, "both.msa")).toBe(params(a, "both.msa"));
    expect(params(b, "data")).toBe(params(a, "data")); // a node given no params
    expect(params(b, "both.ann")).not.toBe(params(a, "both.ann"));
    expect(params(b, "both.ann")?.height).toBe("box");
    // The same normalized content as a normalized document's index; without a cache, new objects.
    const whole = indexDoc(normalize(changed, reg));
    for (const n of whole.order) {
      expect(b.nodes.get(n.ref)?.params).toEqual(n.params);
      expect(b.nodes.get(n.ref)?.node.custom).toEqual(n.node.custom);
    }
    expect(params(indexDoc(d, reg), "both.msa")).not.toBe(params(a, "both.msa"));
  });

  it("keeps them across a registration that leaves them equal; a node it changes gets new ones", () => {
    const r = registries();
    const cache = new WeakMap<object, NormalizedParams>();
    const b = createSeqViewSpec();
    b.data.sequence({ ref: "s", value: "MKVLA" });
    b.view({ ref: "v", label: "V" }).track({ ref: "t", label: "T" }).letters();
    const d = b.getDocument();
    const track = d.root.children?.find((c) => c.ref === "v")?.children?.[0]?.children?.[0];
    track?.children?.push({ kind: "layer", ref: "ext", params: { type: "org.example:bars", fallback: "omit" } });
    const ext = (name: string, options = SimpleParamsSchema({})): LayerTypeEntry => ({
      ...{ name, options, accepts: "none", layouts: ["tracks"], axes: ["blocks"], sizing: "row" },
      ...{ fallback: "omit", channels: [] },
    });
    const a = indexDoc(d, r, cache);
    const params = (x: DocIndex, ref: string) => x.nodes.get(ref)?.params;
    const version = r.version;
    register(r, { layerTypes: [ext("org.example:other")] });
    const same = indexDoc(d, r, cache);
    expect(r.version).toBeGreaterThan(version);
    for (const n of a.order) expect(params(same, n.ref)).toBe(n.params);
    register(r, { layerTypes: [ext("org.example:bars", SimpleParamsSchema({ k: OptionalField(int, 3, "") }))] });
    const next = indexDoc(d, r, cache);
    expect(params(next, "ext")).not.toBe(params(a, "ext"));
    expect(params(next, "ext")?.options).toEqual({ k: 3 });
    expect(params(next, "s")).toBe(params(a, "s"));
  });
});

describe("composeView", () => {
  it("checks scroll groups: sections of one group with different axes", async () => {
    const bad = patched(case10(), [{ op: "set", ref: "both.ann", params: { positions: "msa-gatekeeper" } }]);
    const { set } = await compile(bad, reg, fixtureLoaders());
    const index = indexDoc(bad, reg);
    const sections = ["both.ann", "both.msa"].map((ref) => resolveSection(index, set, ref, reg).section);
    const issues: Issue[] = [];
    composeView("both", sections, [], set, issues);
    expect(issues.map((i) => [i.code, i.severity, i.path])).toEqual([
      ["svs.section.scroll_group", "error", index.nodes.get("both.msa")?.path],
    ]);
  });

  it("takes each layer's track from the section, so reused layers stay right across track moves and filters", async () => {
    const d = case10();
    const { set } = await compile(d, reg, fixtureLoaders());
    const section = indexDoc(d, reg).order.find((n) => n.kind === "section" && n.parent === "alignment")?.ref ?? "";
    const stages = (x: SeqViewSpecDoc) => {
      const index = indexDoc(x, reg);
      const s = resolveSection(index, set, section, reg).section;
      return { index, s, layers: s.layers.map((l) => resolveLayer(s, l.ref, index, set, reg, all).layer) };
    };
    const before = stages(d);
    const tracks = before.index.order.filter((n) => n.kind === "track" && n.parent === section).map((n) => n.ref);
    const layerOf = (type: string) => before.index.order.find((n) => n.kind === "layer" && n.params.type === type);
    const values = layerOf("values") as IndexedNode;
    const reuse = (after: { s: { layers: readonly { ref: string }[] } }) =>
      after.s.layers.map((l) => before.layers.find((x) => x.ref === l.ref) as ResolvedLayer);

    const changes: Patch["ops"][] = [
      [{ op: "move", ref: tracks[0], parent: section, after: "last" }], // a track reorder
      [
        { op: "move", ref: values.ref, parent: tracks[0] }, // a layer to another track without lanes
        { op: "remove", ref: values.parent as string },
      ],
      [{ op: "set", ref: section, params: { positions: "msa-gatekeeper" } }], // a filter view
    ];
    for (const ops of changes) {
      const x = patched(d, ops);
      const after = stages(x);
      expect(after.s.axisKey).toBe(before.s.axisKey);
      const view = composeView("alignment", [after.s], [reuse(after)], set);
      expect(view).toEqual(resolveView(x, set, "alignment", reg, all).view);
    }
    // The reordered tracks: a reused layer's own track index is stale; the composed view's is the section's.
    const moved = stages(patched(d, changes[0]));
    expect(reuse(moved).map((l) => l.track)).not.toEqual(moved.layers.map((l) => l.track));
    // A layer onto the track with lanes reads them: the key changes, so it is resolved again.
    const kinases = tracks.find((t) => before.index.nodes.get(t)?.params.label === "Kinases") as string;
    const onLanes = stages(
      patched(d, [
        { op: "move", ref: values.ref, parent: kinases },
        { op: "remove", ref: values.parent as string },
      ]),
    );
    expect(onLanes.s.axisKey).not.toBe(before.s.axisKey);
  });
});
