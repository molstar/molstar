/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Results (LS2; seqspec:interaction#results): sets lit on `results:<name>` (loci as one target per kind, `Coordinates`
// merged), stepped round with the current one focused and accented on `results:<name>:current`, gone to by index,
// carried to a new set, cleared alone; several sets at once; `results.save` as a location set; and a `band` set painted
// as a highlight where the layout is wrapped (no band in the main pass there).

import { createSeqViewSpec } from "@molstar/seq-core/builder";
import { describe, expect, it } from "vitest";
import { itWithin } from "../../core/_test/model/timing.ts";
import { Pass, SkinToken } from "../src/draw/index.ts";
import type { SeqViewerContext } from "../src/index.ts";
import { opened, preparedRenderers } from "./layout-fixtures.ts";
import { lit, passRects, seqLoci, twoBlocks } from "./overlay-fixtures.ts";

const flush = async () => {
  for (let i = 0; i < 5; i++) await new Promise((r) => setTimeout(r, 0));
};
function doc() {
  const b = createSeqViewSpec({ title: "Results" });
  const a = b.data.sequence({ ref: "a", value: "MKCAYCAKQRQISFVKSHFSRQLEERLGLIEVQ" });
  b.view({ ref: "v", label: "V" }).block({ sequence: a }).track({ label: "L" }).letters();
  return b.getDocument();
}
const channels = (ctx: SeqViewerContext) =>
  [...ctx.highlight.channels.value].filter(([c]) => c.startsWith("results:")).map(([c, e]) => [c, e.style]);
const focus = (ctx: SeqViewerContext) => {
  const f = ctx.focus.current.value.get("v.main");
  return f && [f.start, f.stop];
};

describe("ctx.results", () => {
  it("lights a set as one target, steps round focusing and accenting the current one, goes to one, clears alone", async () => {
    const { ctx } = await opened(doc(), { renderers: preparedRenderers() });
    const set = ctx.graph.assembled.value?.set as never;
    const hits = [seqLoci(set, "a", 2, 3), seqLoci(set, "a", 5, 6), seqLoci(set, "a", 20, 24)];
    ctx.results.set("cys", hits, { labels: ["C3", "C6", "LEER"] });
    ctx.results.set("host", [{ systems: ["x"], runs: Int32Array.of(0, 1, 2) }], { style: "accent" });
    expect(channels(ctx)).toEqual([
      ["results:cys", "highlight"],
      ["results:host", "accent"],
    ]);
    expect(ctx.results.current.value).toMatchObject({ name: "host", count: 1, index: -1 });
    ctx.results.step(1, { name: "cys" });
    expect([ctx.results.current.value?.index, focus(ctx)]).toEqual([0, [2, 3]]);
    ctx.results.step(-1);
    expect([ctx.results.current.value?.index, focus(ctx)]).toEqual([2, [20, 24]]);
    ctx.results.step(1);
    expect(ctx.results.current.value?.index).toBe(0);
    ctx.results.go("cys", 1);
    expect([ctx.results.sets.value.get("cys")?.index, focus(ctx)]).toEqual([1, [5, 6]]);
    expect(channels(ctx)).toContainEqual(["results:cys:current", "accent"]);
    ctx.results.clear("cys");
    expect(channels(ctx)).toEqual([["results:host", "accent"]]);
    expect([...ctx.results.sets.value.keys()]).toEqual(["host"]);
    ctx.dispose();
  });

  itWithin(
    "takes thousands of results as it takes a few: one lit target, stepping unchanged, in a bounded time",
    async (within) => {
      const b = createSeqViewSpec({ title: "Big" });
      const a = b.data.sequence({ ref: "a", value: "ACDEFGHIKLMNPQRSTVWY".repeat(1000) });
      b.view({ ref: "v", label: "V" }).block({ sequence: a }).track({ label: "L" }).letters();
      const { ctx } = await opened(b.getDocument(), { renderers: preparedRenderers() });
      const set = ctx.graph.assembled.value?.set as never;
      // Every second residue of 20,000: 10,000 results, as a one-letter search of a long protein makes.
      const hits = Array.from({ length: 10_000 }, (_, k) => seqLoci(set, "a", 2 * k, 2 * k + 1));
      within(1500, () => {
        ctx.results.set("many", hits, { labels: hits.map((_, k) => `hit ${k}`) });
        ctx.results.step(1, { name: "many" });
        ctx.results.go("many", 9_999);
      });
      expect(ctx.results.sets.value.get("many")?.count).toBe(10_000);
      expect(channels(ctx).map(([c]) => c)).toEqual(["results:many", "results:many:current"]);
      expect(focus(ctx)).toEqual([19_998, 19_999]);
      ctx.dispose();
    },
  );

  it("going to the current result again clears it (no accent, the focus stays); stepping starts over", async () => {
    const { ctx } = await opened(doc(), { renderers: preparedRenderers() });
    const set = ctx.graph.assembled.value?.set as never;
    ctx.results.set("cys", [seqLoci(set, "a", 2, 3), seqLoci(set, "a", 5, 6), seqLoci(set, "a", 20, 24)]);
    ctx.results.go("cys", 1);
    expect([ctx.results.sets.value.get("cys")?.index, focus(ctx)]).toEqual([1, [5, 6]]);
    expect(channels(ctx).map(([c]) => c)).toContain("results:cys:current");
    // The current one again: unselected, the focus where it was.
    ctx.results.go("cys", 1);
    expect([ctx.results.sets.value.get("cys")?.index, focus(ctx)]).toEqual([-1, [5, 6]]);
    expect(channels(ctx).map(([c]) => c)).toEqual(["results:cys"]);
    // Previous from none is the last, Next the first; a lone result's step does not clear it.
    ctx.results.step(-1, { name: "cys" });
    expect(ctx.results.sets.value.get("cys")?.index).toBe(2);
    ctx.results.clearCurrent("cys");
    ctx.results.step(1, { name: "cys" });
    expect(ctx.results.sets.value.get("cys")?.index).toBe(0);
    ctx.results.set("one", [seqLoci(set, "a", 7, 8)]);
    ctx.results.step(1, { name: "one" });
    ctx.results.step(1, { name: "one" });
    expect(ctx.results.sets.value.get("one")?.index).toBe(0);
    ctx.dispose();
  });

  it("carries its results to a new set; results.save keeps them as a location set with a highlights entry", async () => {
    const { ctx } = await opened(doc(), { renderers: preparedRenderers() });
    const before = ctx.graph.assembled.value?.set as never;
    ctx.results.set("two", [seqLoci(before, "a", 2, 3), seqLoci(before, "a", 10, 14)], { labels: ["C", "QRQI"] });
    expect(ctx.actions.can("results.save", { section: "v.main", name: "nope" })).toBe(false);
    expect(ctx.actions.run("results.save", { section: "v.main", name: "two" })).toBe(true);
    await flush();
    const after = ctx.graph.assembled.value?.set;
    expect(after).not.toBe(before);
    expect((ctx.results.results("two")[0] as { set: unknown }).set).toBe(after);
    const node = (ref: string) =>
      (ctx.document.shown.value?.root.children ?? [])
        .flatMap((c) => [c, ...(c.children ?? [])])
        .find((c) => c.ref === ref);
    expect(node("v.main.two")?.params?.value).toEqual([
      { target: "a/unit", start: 2, stop: 3, id: "r1", label: "C" },
      { target: "a/unit", start: 10, stop: 14, id: "r2", label: "QRQI" },
    ]);
    expect(node("v.main")?.params?.highlights).toEqual([{ data: "v.main.two" }]);
    // Saved again: a ref of its own.
    ctx.actions.run("results.save", { section: "v.main", name: "two" });
    await flush();
    expect(node("v.main.two-2")).toBeDefined();
    ctx.dispose();
  });

  it("a new set leaves `current` with the set last set or stepped; a set whose count changed loses its accent (VS-09)", async () => {
    const { ctx } = await opened(doc(), { renderers: preparedRenderers() });
    const set = ctx.graph.assembled.value?.set as never;
    ctx.results.set("search", [seqLoci(set, "a", 2, 3), seqLoci(set, "a", 5, 6)]);
    ctx.results.set("other", [seqLoci(set, "a", 10, 14)]);
    ctx.results.step(1, { name: "search" });
    expect(ctx.results.current.value?.name).toBe("search");
    const patch = (value: string) =>
      ctx.applyPatch({ format: "seq_view_spec_patch", ops: [{ op: "set", ref: "a", params: { value } }] });
    await patch("MKCAYCAKQRQISFVKSHFSRQLEERLGLIEVK");
    await flush();
    expect(ctx.results.current.value?.name).toBe("search");
    // An unnamed step steps the current set.
    ctx.results.step(1);
    expect(ctx.results.sets.value.get("search")?.index).toBe(1);
    expect(ctx.results.sets.value.get("other")?.index).toBe(-1);
    // A shorter sequence drops the second result: no current one, and its accent goes with it.
    await patch("MKCA");
    await flush();
    expect(ctx.results.sets.value.get("search")).toMatchObject({ count: 1, index: -1 });
    expect(channels(ctx).map(([c]) => c)).not.toContain("results:search:current");
    ctx.dispose();
  });

  it("a band set in a wrapped section is painted as a highlight on the overlay", async () => {
    const t = await lit(twoBlocks({ wrapped: true }), { width: 600, height: 300 });
    const band = (r: number[][]) => r.filter((x) => x[4] === SkinToken.highlight);
    expect(band(passRects(t.frame(undefined, 7), Pass.overlay))).toEqual([]);
    t.ctx.highlight.set("site", seqLoci(t.set, "a", 20, 30), { style: "band" });
    const r = t.frame(undefined, 7);
    expect(band(passRects(r, Pass.main))).toEqual([]);
    expect(band(passRects(r, Pass.overlay)).length).toBeGreaterThan(0);
  });
});
