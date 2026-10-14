/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The graph's lifecycle and robustness (VS4's graph half): a throwing `key` or `inputs` settles, issues are capped per
// node, a superseded run's result is disposed, a cleared preview reruns nothing. Each case also ran against the old
// stack until SW1 (its parity layer, deleted).

import { createSeqViewSpec } from "@molstar/seq-core/builder";
import { describe, expect, it } from "vitest";
import {
  createSeqViewer,
  DefaultNodeKinds,
  DefaultRegistry,
  type NodeKind,
  type SeqViewerContext,
} from "../src/index.js";
import { Issues } from "../src/ui.js";
import { domMoved, manualLoaders, patch, smallDoc } from "./graph-fixtures.js";

/** A context whose `annotation` kind is replaced by `over(base)`. */
function withKind(over: (base: NodeKind) => NodeKind): SeqViewerContext {
  const base = DefaultNodeKinds.find((k) => k.name === "annotation") as NodeKind;
  const kinds = [...DefaultNodeKinds.filter((k) => k !== base), over(base)];
  return createSeqViewer({ registry: [{ ...DefaultRegistry, nodeKinds: kinds }], behaviors: [] });
}
const codes = (ctx: SeqViewerContext, ref: string) => ctx.graph.status(ref)?.issues.map((i) => i.code);

describe("a run always settles", () => {
  for (const part of ["key", "inputs"] as const) {
    it(`a node kind whose \`${part}\` throws is that node's error; the update settles`, async () => {
      const ctx = withKind((base) => ({
        ...base,
        [part]: (...a: Parameters<NonNullable<NodeKind[typeof part]>>) => {
          if ((part === "key" ? (a[0] as { ref: string }).ref : a[0]) === "dom") throw new Error("too old");
          return (base[part] as (...x: unknown[]) => void)(...a);
        },
      }));
      const r = await ctx.setDocument(smallDoc());
      expect(r.ok).toBe(true);
      expect(ctx.graph.status("dom")?.state).toBe("error");
      expect(codes(ctx, "dom")).toEqual(["svs.source.load"]);
      expect(ctx.graph.status("site")?.state).toBe("ok");
      expect(ctx.graph.assembled.value?.view?.ref).toBe("v");
    });
  }
});

describe("issues", () => {
  it("200,000 issues from one node keep 100 and count the rest in one", async () => {
    const ctx = withKind((base) => ({
      ...base,
      run(rc) {
        if (rc.ref === "dom")
          for (let i = 0; i < 200_000; i++)
            rc.stage.report({ code: "svs.source.record", severity: "warning", path: "", message: `${i}` });
        return base.run(rc);
      },
    }));
    await ctx.setDocument(smallDoc());
    const issues = ctx.graph.status("dom")?.issues ?? [];
    expect(issues.length).toBe(101);
    expect(issues[100]).toMatchObject({ code: "seqv.issues.truncated", message: "and 199900 more (199900 warnings)" });
    expect(ctx.graph.issues.value.filter((i) => i.ref === "dom").length).toBe(101);
  });

  it("past the cap, a node keeps its errors before its warnings, and the count says how many of each (GV-08)", async () => {
    const ctx = withKind((base) => ({
      ...base,
      run(rc) {
        if (rc.ref === "dom") {
          for (let i = 0; i < 150; i++)
            rc.stage.report({ code: "svs.source.record", severity: "warning", path: "", message: `${i}` });
          for (let i = 0; i < 3; i++)
            rc.stage.report({ code: "svs.source.load", severity: "error", path: "", message: `e${i}` });
        }
        return base.run(rc);
      },
    }));
    await ctx.setDocument(smallDoc());
    const issues = ctx.graph.status("dom")?.issues ?? [];
    expect(issues.length).toBe(101);
    expect(issues.slice(0, 3).map((i) => i.message)).toEqual(["e0", "e1", "e2"]);
    expect(issues[100]).toMatchObject({ code: "seqv.issues.truncated", message: "and 53 more (53 warnings)" });
  });

  it("the store of issues not from a node evicts the least severe first, the oldest among equals (GV-08)", () => {
    const issues = new Issues();
    const add = (severity: "error" | "warning", message: string) =>
      issues.add({ code: "seqv.state.error", severity, path: "", message, ref: null });
    add("error", "boom");
    for (let i = 0; i < Issues.Limit; i++) add("warning", `w${i}`);
    const held = issues.current.value;
    expect(held.length).toBe(Issues.Limit);
    expect(held.filter((i) => i.severity === "error").map((i) => i.message)).toEqual(["boom"]);
    expect(held.some((i) => i.message === "w0")).toBe(false);
    expect(held.some((i) => i.message === `w${Issues.Limit - 1}`)).toBe(true);
    // Full of errors, a warning is the least severe: it is the one left out.
    const errors = new Issues();
    for (let i = 0; i < Issues.Limit; i++)
      errors.add({ code: "seqv.state.error", severity: "error", path: "", message: `e${i}`, ref: null });
    errors.add({ code: "seqv.state.error", severity: "warning", path: "", message: "late", ref: null });
    expect(errors.current.value.every((i) => i.severity === "error")).toBe(true);
    issues.dispose();
    errors.dispose();
  });
});

describe("superseded work", () => {
  it("a superseded async node's result is disposed, never published", async () => {
    const disposed: unknown[] = [];
    const held: (() => void)[] = [];
    const ctx = withKind((base) => ({
      ...base,
      run: (rc) => (rc.ref === "dom" ? new Promise((done) => held.push(() => done(base.run(rc)))) : base.run(rc)),
      dispose: (r) => void disposed.push(r),
    }));
    const first = ctx.setDocument(smallDoc());
    await new Promise((r) => setTimeout(r, 0));
    const second = ctx.setDocument(domMoved());
    await new Promise((r) => setTimeout(r, 0));
    expect(held.length).toBe(2);
    held[1]();
    await Promise.all([first, second]);
    held[0](); // the first run, superseded: lands after
    await new Promise((r) => setTimeout(r, 0));
    expect(disposed.length).toBe(1);
    expect(disposed[0]).not.toBe(ctx.graph.assembled.value?.set);
  });

  it("a preview cleared before its run settles loads nothing more", async () => {
    const io = manualLoaders();
    const ctx = createSeqViewer(undefined, { loaders: io.loaders });
    const b = createSeqViewSpec();
    b.data.sequence({ ref: "s", source: { uri: "https://x.test/a.fa", format: "fasta" } });
    b.view({ ref: "v", label: "V" }).track().letters();
    const done = ctx.setDocument(b.getDocument());
    await new Promise((r) => setTimeout(r, 0));
    io.resolve("https://x.test/a.fa", ">s\nMKV\n");
    await done;
    const preview = patch([{ op: "set", ref: "s.src", params: { uri: "https://x.test/b.fa" } }]);
    void ctx.applyPatch(preview, { layer: "preview" });
    await new Promise((r) => setTimeout(r, 0));
    expect(io.requests).toEqual(["https://x.test/a.fa", "https://x.test/b.fa"]);
    ctx.document.setPreview(null);
    await ctx.graph.update();
    io.resolve("https://x.test/b.fa", ">s\nMKVLAT\n"); // after it was dropped
    await new Promise((r) => setTimeout(r, 0));
    await ctx.graph.update();
    expect(io.requests.length).toBe(2);
    expect(ctx.graph.assembled.value?.set.molecules[0]?.length).toBe(3);
  });
});
