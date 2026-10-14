/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Document layers, patches in and out, the history (seqspec:viewer#the-document, seqspec:viewer#history) and data
// without a document (seqspec:viewer#data-mode). The random test also enters and leaves snapshots, previews and takes
// foreign documents (VS2), and checks the shown document against the layers composed afresh.

import { applyPatch, canonical, type Op, type Patch, type SeqViewSpecDoc, snapshotDoc } from "@molstar/seq-core/format";
import { type SeqSet, SeqSets, type SetInput } from "@molstar/seq-core/model";
import { describe, expect, it } from "vitest";
import type { DocumentPort, PatchEvent, Routed, SeqIssue, SeqViewerContext } from "../src/index.js";
import { patch, refsRan, smallDoc, viewer } from "./graph-fixtures.js";

const hide = (ref: string, hidden = true): Op => ({ op: "set", ref, params: { hidden } });
const nodeOf = (doc: SeqViewSpecDoc | null, ref: string) => {
  const find = (n: SeqViewSpecDoc["root"]): SeqViewSpecDoc["root"] | undefined =>
    n.ref === ref ? n : n.children?.map(find).find((x) => x);
  return doc ? find(doc.root) : undefined;
};
const hidden = (ctx: SeqViewerContext, ref: string, layer: "base" | "shown" = "shown") =>
  nodeOf(ctx.document[layer].value, ref)?.params?.hidden === true;
function issuesOf(ctx: SeqViewerContext) {
  const out: SeqIssue[] = [];
  ctx.events.issue.subscribe((i) => out.push(i));
  return out;
}

describe("document layers", () => {
  it("base → snapshot → local → preview, untouched subtrees shared through every layer", async () => {
    const { ctx } = viewer();
    const d = smallDoc();
    d.snapshots = [{ ref: "quiet", ops: [hide("v.ruler")] }];
    await ctx.setDocument(d);
    const base = ctx.document.base.value as SeqViewSpecDoc;
    await ctx.snapshots.apply("quiet");
    expect((await ctx.applyPatch(patch([hide("v.site")]), { layer: "local" })).ok).toBe(true);
    expect((await ctx.document.setPreview(patch([hide("v.dom")]))).ok).toBe(true);
    const shown = ctx.document.shown.value;
    for (const r of ["v.ruler", "v.site", "v.dom"])
      expect([hidden(ctx, r), hidden(ctx, r, "base")]).toEqual([true, false]);
    expect(shown?.root.children?.[0]).toBe(base.root.children?.[0]); // the data subtree
    expect(nodeOf(shown, "v.letters")).toBe(nodeOf(base, "v.letters"));
    expect(ctx.document.local.value.map((l) => l.snapshot)).toEqual(["quiet"]);
    // The preview is drawn, never committed; clearing layers goes back down.
    await ctx.document.setPreview(null);
    expect(hidden(ctx, "v.dom")).toBe(false);
    ctx.document.clearLocal();
    await ctx.snapshots.apply(null);
    expect(ctx.document.shown.value).toBe(base);
    await ctx.graph.update();
    expect(ctx.graph.assembled.value?.view?.resolved.tracks.hidden.every((h) => h === 0)).toBe(true);
  });

  it("a base patch re-applies local patches; one that no longer applies is dropped and reported", async () => {
    const { ctx } = viewer();
    const issues = issuesOf(ctx);
    await ctx.setDocument(smallDoc());
    await ctx.applyPatch(patch([hide("v.site")]), { layer: "local" });
    await ctx.applyPatch(patch([hide("v.dom")]), { layer: "local" });
    await ctx.applyPatch(patch([hide("v.letters")]));
    expect(["v.site", "v.dom", "v.letters"].map((r) => hidden(ctx, r))).toEqual([true, true, true]);
    await ctx.applyPatch(patch([{ op: "remove", ref: "v.site" }]));
    expect(ctx.document.local.value.map((l) => l.patch.ops[0])).toEqual([hide("v.dom")]);
    expect(hidden(ctx, "v.dom")).toBe(true);
    expect(issues.map((i) => i.code)).toEqual(["seqv.patch.dropped"]);
  });

  it("a failed applyPatch changes nothing; an invalid preview is ok: false and not drawn", async () => {
    const { ctx } = viewer();
    await ctx.setDocument(smallDoc());
    const [base, shown, assembled] = [ctx.document.base.value, ctx.document.shown.value, ctx.graph.assembled.value];
    const r = await ctx.applyPatch(patch([hide("v.ruler"), { op: "remove", ref: "nowhere" }]));
    expect([r.ok, r.issues.map((i) => i.code)]).toEqual([false, ["svs.patch.ref"]]);
    expect(ctx.document.base.value).toBe(base);
    expect(ctx.document.shown.value).toBe(shown);
    expect(ctx.graph.assembled.value).toBe(assembled);
    expect((await ctx.applyPatch(patch([{ op: "remove", ref: "nowhere" }]), { layer: "local" })).ok).toBe(false);
    expect(ctx.document.local.value).toEqual([]);
    expect(ctx.history.canUndo.value).toBe(false);
    for (const p of [
      () => ctx.document.setPreview(patch([{ op: "remove", ref: "nowhere" }])),
      () => ctx.applyPatch(patch([{ op: "set", ref: "v.dom.l", params: { type: "nothing" } }]), { layer: "preview" }),
    ]) {
      const pr = await p();
      expect(pr.ok).toBe(false);
      expect(pr.issues.length).toBeGreaterThan(0);
      expect([ctx.document.preview.value, ctx.document.shown.value]).toEqual([null, shown]);
    }
  });

  it("an invalid document is refused and nothing changes", async () => {
    const { ctx } = viewer();
    await ctx.setDocument(smallDoc());
    const base = ctx.document.base.value;
    const bad = smallDoc();
    (nodeOf(bad, "v.main") as { params?: object }).params = { arrangement: "diagonal" };
    const r = await ctx.setDocument(bad);
    expect(r.ok).toBe(false);
    expect(r.issues[0].code).toMatch(/^svs\.params/);
    expect(ctx.document.base.value).toBe(base);
  });
});

describe("patches out", () => {
  it("the default port commits and keeps history; events carry the routing", async () => {
    const { ctx } = viewer();
    const events: PatchEvent[] = [];
    ctx.events.patch.subscribe((e) => events.push(e));
    await ctx.setDocument(smallDoc());
    const r = await ctx.document.propose(patch([hide("v.site")]), "Hide site", "pointer");
    expect(r).toEqual({ routed: "commit", ok: true, issues: [] });
    expect(hidden(ctx, "v.site", "base")).toBe(true);
    expect(events.map((e) => [e.label, e.routed, e.origin])).toEqual([["Hide site", "commit", "pointer"]]);
    expect(ctx.history.canUndo.value).toBe(true);
    expect(await ctx.history.undo()).toBe(true);
    expect(hidden(ctx, "v.site", "base")).toBe(false);
    expect([ctx.history.canUndo.value, ctx.history.canRedo.value]).toEqual([false, true]);
  });

  it("a patch on the implied `data` or `root` commits too: they are the base's nodes (VS-23)", async () => {
    const { ctx } = viewer();
    await ctx.setDocument(smallDoc());
    const insert: Op = {
      op: "insert",
      parent: "data",
      node: { kind: "sequence", ref: "c", params: { value: "MKV" } },
    } as Op;
    expect((await ctx.document.propose(patch([insert]), "Add", "api")).routed).toBe("commit");
    expect(nodeOf(ctx.document.base.value, "c")).toBeTruthy();
    const root: Op = { op: "set", ref: "root", params: { default_view: "w" } };
    expect((await ctx.document.propose(patch([root]), "Default", "api")).routed).toBe("commit");
    expect(ctx.document.base.value?.root.params?.default_view).toBe("w");
    // A node the base does not have (a local insert's) is still masked.
    const later: Op = { op: "set", ref: "nowhere", params: { label: "x" } };
    expect((await ctx.document.propose(patch([later]), "Set", "api")).routed).toBe("local");
  });

  it("a host's port: local patches stay in the viewer, refusals change nothing, a stale base is rejected", async () => {
    let route: Routed = "local";
    const sent: Patch[] = [];
    const port: DocumentPort = {
      route: () => route,
      commit: async (p) => {
        sent.push(p);
        return { ok: true, issues: [] };
      },
      history: false,
    };
    const { ctx } = viewer({ document: port });
    const issues = issuesOf(ctx);
    await ctx.setDocument({ ...smallDoc(), revision: "r1" });
    expect((await ctx.document.propose(patch([hide("v.site")]), "Hide", "pointer")).routed).toBe("local");
    expect([hidden(ctx, "v.site"), hidden(ctx, "v.site", "base")]).toEqual([true, false]);
    route = "refuse";
    const refused = await ctx.document.propose(patch([hide("v.dom")]), "Hide", "pointer");
    expect([refused.ok, refused.issues.map((i) => i.code)]).toEqual([false, ["seqv.port.refused"]]);
    route = "commit";
    const stale = await ctx.document.propose(patch([hide("v.dom")], { base: "r0" }), "Hide", "pointer");
    expect([stale.ok, stale.issues.map((i) => i.code)]).toEqual([false, ["svs.patch.base"]]);
    route = "local";
    const staleLocal = await ctx.document.propose(patch([hide("v.dom")], { base: "r0" }), "Hide", "pointer");
    expect(staleLocal.issues.map((i) => i.code)).toEqual(["svs.patch.base"]);
    expect(hidden(ctx, "v.dom")).toBe(false);
    route = "commit";
    expect((await ctx.document.propose(patch([hide("v.dom")]), "Hide", "api")).ok).toBe(true);
    // Forwarded with the base's revision; the host passes it back by setDocument.
    expect(sent).toEqual([{ format: "seq_view_spec_patch", base: "r1", ops: [hide("v.dom")] }]);
    expect(ctx.document.local.value.map((l) => l.patch.base)).toEqual([undefined]);
    expect(issues.map((i) => i.code)).toEqual(["seqv.port.refused", "svs.patch.base", "svs.patch.base"]);
    // Only the local patch is in the history (the host owns undo of its commits).
    expect(await ctx.history.undo()).toBe(true);
    expect(await ctx.history.undo()).toBe(false);
    expect(ctx.document.local.value).toEqual([]);
  });
});

/** A small deterministic generator. */
function random(seed: number) {
  let s = seed + 1;
  return () => {
    s = (s * 1103515245 + 12345) % 2147483648;
    return s / 2147483648;
  };
}

/** The shown document composed afresh from the layers: base, snapshot, the local patches it shows, the preview. */
function composed(ctx: SeqViewerContext): string {
  const svs = ctx.registry.svs;
  const d = ctx.document;
  const snap = d.snapshot.value;
  let doc = d.base.value as SeqViewSpecDoc;
  if (snap !== null) doc = snapshotDoc(doc, snap, svs).doc;
  for (const l of d.local.value) {
    const r = l.snapshot === null || l.snapshot === snap ? applyPatch(doc, l.patch, svs) : null;
    if (r?.inverse) doc = r.doc;
  }
  const p = d.preview.value && applyPatch(doc, d.preview.value, svs);
  return canonical((p?.inverse ? p.doc : doc).root);
}

describe("history", () => {
  /** `count` seeded sequences of 14 steps each (from seed 0), every state checked, then undone to the start. */
  async function roundTrips(count: number) {
    const tracks = ["v.ruler", "v.letters", "v.dom", "v.site"];
    const doc = smallDoc();
    doc.snapshots = [
      { ref: "s1", ops: [hide("v.ruler")] },
      { ref: "s2", ops: [{ op: "set", ref: "v.dom", params: { label: "S2" } }] },
    ];
    for (let seed = 0; seed < count; seed++) {
      const rand = random(seed);
      const pick = <T>(xs: readonly T[]) => xs[Math.floor(rand() * xs.length)];
      let route: Routed = "commit";
      const holder: { ctx?: SeqViewerContext } = {};
      const port: DocumentPort = {
        route: () => route,
        commit: (p) => (holder.ctx as SeqViewerContext).applyPatch(p),
        history: true,
      };
      const { ctx } = viewer({ document: port });
      holder.ctx = ctx;
      await ctx.setDocument(doc);
      // What the history restores: the base and the local layer, tags included (not the snapshot shown or the preview).
      const kept = () => canonical([ctx.document.base.value?.root, ctx.document.local.value]);
      const done = [kept()];
      const undone: string[] = [];
      for (let step = 0; step < 14; step++) {
        const r = rand();
        if (r < 0.2 && ctx.history.canUndo.value) {
          expect(await ctx.history.undo()).toBe(true);
          undone.push(done.pop() as string);
        } else if (r < 0.3 && ctx.history.canRedo.value) {
          expect(await ctx.history.redo()).toBe(true);
          done.push(undone.pop() as string);
        } else if (r < 0.4) {
          await ctx.snapshots.apply(pick([null, "s1", "s2"]));
        } else if (r < 0.47) {
          const p = rand() < 0.5 ? null : patch([{ op: "set", ref: pick(tracks), params: { label: `P${step}` } }]);
          expect((await ctx.document.setPreview(p)).ok).toBe(true);
        } else if (r < 0.52) {
          // A foreign document: the history is cleared; the local layer stays.
          const next = structuredClone(ctx.document.base.value as SeqViewSpecDoc);
          const n = nodeOf(next, pick(tracks)) as { params?: Record<string, unknown> };
          n.params = { ...n.params, label: `F${step}` };
          expect((await ctx.setDocument(next)).ok).toBe(true);
          expect([ctx.history.canUndo.value, ctx.history.canRedo.value]).toEqual([false, false]);
          done.splice(0, done.length, kept());
          undone.length = 0;
        } else {
          route = rand() < 0.5 ? "commit" : "local";
          const t = pick(tracks);
          const op: Op =
            rand() < 0.5
              ? { op: "set", ref: t, params: rand() < 0.5 ? { hidden: rand() < 0.5 } : { label: `L${step}` } }
              : {
                  op: "move",
                  ref: t,
                  parent: "v.main",
                  after: pick(["first", "last", ...tracks.filter((x) => x !== t)]),
                };
          expect((await ctx.document.propose(patch([op]), `step ${step}`, "api")).ok).toBe(true);
          done.push(kept());
          undone.length = 0;
        }
        expect(kept()).toBe(done.at(-1));
        expect(canonical(ctx.document.shown.value?.root)).toBe(composed(ctx));
      }
      while (ctx.history.canUndo.value) {
        expect(await ctx.history.undo()).toBe(true);
        done.pop();
        expect(kept()).toBe(done.at(-1));
      }
      expect(done.length).toBe(1);
      ctx.dispose();
    }
  }

  it("undo and redo round-trip 30 random sequences with snapshots, previews and foreign documents", () =>
    roundTrips(30));

  it(
    "undo and redo round-trip 300 random sequences with snapshots, previews and foreign documents",
    { tags: ["slow"] },
    () => roundTrips(300),
  );
});

describe("data without a document", () => {
  const inputs = (): SetInput => ({
    molecules: [
      { key: "p", residues: "MKVLATGHEK", ownUnit: true },
      { key: "q", residues: "MKVLA", ownUnit: true },
    ],
    annotations: [
      { key: "dom", geometry: "ranges", target: "p/unit", columns: { start: [0], stop: [4] } },
      { key: "site", geometry: "ranges", target: "q/unit", columns: { start: [1], stop: [2] } },
    ],
  });

  it("build inputs: one node per entry, a default view; a rebuild with kept entries reruns what changed", async () => {
    const { ctx } = viewer();
    const first = inputs();
    const r = await ctx.setData({ inputs: first });
    expect(r.ok).toBe(true);
    const a0 = ctx.graph.assembled.value;
    const view = a0?.view?.resolved;
    const names = (ids: ArrayLike<number>) => Array.from(ids, (id) => a0?.set.strings.get(id));
    expect(ctx.views.list.value).toEqual([{ ref: "view", label: "View" }]);
    expect(names(view?.sections.ref ?? [])).toEqual(["view.main"]);
    expect(names(view?.tracks.ref ?? [])).toEqual(["view.ruler", "view.letters", "view.dom", "view.site"]);
    expect(view?.blocks.count).toBe(2);
    const next: SetInput = {
      ...first,
      annotations: [first.annotations?.[0], { ...first.annotations?.[1], columns: { start: [0], stop: [3] } }] as never,
    };
    await ctx.setData({ inputs: next });
    // The changed entry, the assemblies over annotations, the section, the site layer's stages and the view.
    expect(refsRan(ctx)).toEqual([
      ...["#assemble", "#assemble.mapped", "#assemble.stream", "#prepare:view.ranges-2", "#theme:view.ranges-2"],
      "#view:view",
      ...["site", "view.main", "view.ranges-2"],
    ]);
    const a1 = ctx.graph.assembled.value;
    expect(a1?.set.lineage).toBe(a0?.set.lineage);
    expect(a1?.set.annotations[0]).toBe(a0?.set.annotations[0]);
  });

  it("a whole set is one node; views by the builder; a new set reruns everything; a document starts over", async () => {
    const { ctx } = viewer();
    const s1 = SeqSets.build(inputs()).set;
    await ctx.setData(
      { set: s1 },
      { views: (b) => b.view({ ref: "mine", label: "Mine" }).track({ ref: "mine.t" }).ranges({ data: "dom" }) },
    );
    expect(ctx.graph.assembled.value?.set).toBe(s1);
    expect(ctx.graph.status("#set")?.state).toBe("ok");
    expect(ctx.graph.assembled.value?.view?.ref).toBe("mine");
    const s2: SeqSet = SeqSets.build(inputs()).set;
    await ctx.setData(
      { set: s2 },
      { views: (b) => b.view({ ref: "mine", label: "Mine" }).track({ ref: "mine.t" }).ranges({ data: "dom" }) },
    );
    expect(refsRan(ctx)).toEqual([
      "#assemble",
      "#assemble.base",
      "#assemble.mapped",
      "#assemble.stream",
      "#prepare:mine.ranges",
      "#set",
      "#theme:mine.ranges",
      "#view:mine",
      "mine.main",
      "mine.ranges",
    ]);
    expect(ctx.graph.assembled.value?.set).toBe(s2);
    await ctx.setDocument(smallDoc());
    const set = ctx.graph.assembled.value?.set as SeqSet;
    expect([set === s2, set.lineage === s2.lineage, ctx.views.current.value]).toEqual([false, false, "v"]);
  });
});
