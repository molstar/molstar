/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Location sets (LS1; seqspec:format#sections, seqspec:interaction#channels, seqspec:viewer#actions): a section's
// `highlights` as channels `set:<section>/<n>` from the shown, validated document (one per entry, lit in its section
// only, set again only when its entry or the set changes, cleared when it goes, another owner's channel untouched; a
// query node as data), and the bookmark actions (`add` inserts the set and its `highlights` entry, then appends; `step`
// focuses and accents the next item in its section; `remove`).

import { createSeqViewSpec } from "@molstar/seq-core/builder";
import type { SeqViewSpecDoc } from "@molstar/seq-core/format";
import { describe, expect, it } from "vitest";
import { Q } from "../../core/src/query/index.js";
import { bookmarkAt, ChannelMarks, type SeqViewerContext } from "../src/index.js";
import { bookmarkIds } from "../src/location-sets.js";
import { opened, preparedRenderers } from "./layout-fixtures.js";
import { seqLoci } from "./overlay-docs.js";

const flush = async () => {
  for (let i = 0; i < 5; i++) await new Promise((r) => setTimeout(r, 0));
};
const P = { renderers: preparedRenderers() };

/** One sequence shown by two sections; `sites` has items s1 and s2, `cys` is a query; `one` lights them. */
function doc(highlights: unknown[] | null = [{ data: "sites" }, { data: "sites", items: ["s2"], style: "accent" }]) {
  const b = createSeqViewSpec({ title: "Sets" });
  const a = b.data.sequence({ ref: "a", value: "MKCAYCAKQRQISFVKSHFSRQLEERLGLIEVQ" });
  b.data.annotation({
    ref: "sites",
    geometry: "ranges",
    target: a,
    value: { start: [2, 9], stop: [4, 12], id: ["s1", "s2"] },
  });
  b.data.computed({ ref: "cys", provider: "query", query: Q.read("(and a (letters C))") });
  const v = b.view({ ref: "v", label: "V" });
  v.section({ ref: "one", ...(highlights ? { highlights } : {}) } as never)
    .block({ sequence: a })
    .track({ label: "L" })
    .letters();
  v.section({ ref: "two" }).block({ sequence: a }).track({ label: "L" }).letters();
  return b.getDocument() as SeqViewSpecDoc;
}
const sets = (ctx: SeqViewerContext) =>
  [...ctx.highlight.channels.value].filter(([c]) => c.startsWith("set:")).map(([c, e]) => [c, e.style, e.section]);
/** Section s's runs of one mark set, as [start, stop] pairs of its first block. */
function runs(ch: ChannelMarks, s: number, which: "band" | "accent" | "highlight") {
  const m = ch.marks(s)?.[which];
  const out: number[][] = [];
  for (let q = m ? m.from(0) : 0; m && q < m.to(0); q++) out.push([m.start[q], m.stop[q]]);
  return out;
}

describe("a section's highlights as channels", () => {
  it("one channel per entry, `band` by default, lit in its own section only", async () => {
    const { ctx } = await opened(doc(), P);
    expect(sets(ctx)).toEqual([
      ["set:one/0", "band", "one"],
      ["set:one/1", "accent", "one"],
    ]);
    const ch = new ChannelMarks(ctx);
    expect(runs(ch, 0, "band")).toEqual([
      [2, 4],
      [9, 12],
    ]);
    expect(runs(ch, 0, "accent")).toEqual([[9, 12]]);
    expect([runs(ch, 1, "band"), runs(ch, 1, "accent")]).toEqual([[], []]);
    ch.dispose();
    ctx.dispose();
  });

  it("a change sets its own channel alone; a removed entry is cleared; another owner's channel stays", async () => {
    const { ctx } = await opened(doc(), P);
    ctx.highlight.set("host", seqLoci(ctx.graph.assembled.value?.set as never, "a", 0, 1), { style: "accent" });
    const [first, second] = [
      ctx.highlight.channels.value.get("set:one/0"),
      ctx.highlight.channels.value.get("set:one/1"),
    ];
    const set = ctx.graph.assembled.value?.set;
    await ctx.applyPatch({
      format: "seq_view_spec_patch",
      ops: [
        { op: "set", ref: "one", params: { highlights: [{ data: "sites" }, { data: "cys", style: "highlight" }] } },
      ],
    });
    await flush();
    // The set is the same (a section param changed): the unchanged entry keeps its channel's object.
    expect(ctx.graph.assembled.value?.set).toBe(set);
    expect(ctx.highlight.channels.value.get("set:one/0")).toBe(first);
    expect(ctx.highlight.channels.value.get("set:one/1")).not.toBe(second);
    expect(sets(ctx)).toEqual([
      ["set:one/0", "band", "one"],
      ["set:one/1", "highlight", "one"],
    ]);
    // A query node is data as any annotation: the cysteines at 2 and 5.
    const ch = new ChannelMarks(ctx);
    expect(runs(ch, 0, "highlight")).toEqual([
      [2, 3],
      [5, 6],
    ]);
    await ctx.applyPatch({
      format: "seq_view_spec_patch",
      ops: [{ op: "set", ref: "one", params: { highlights: [] } }],
    });
    await flush();
    expect(sets(ctx)).toEqual([]);
    expect(ctx.highlight.channels.value.has("host")).toBe(true);
    ch.dispose();
    ctx.dispose();
  });

  it("an entry's items changed to none, or to one id holding a newline, lights as a fresh open does (GV-03)", async () => {
    const patched = async (from: unknown[], to: unknown[]) => {
      const { ctx } = await opened(doc(from), P);
      await ctx.applyPatch({
        format: "seq_view_spec_patch",
        ops: [{ op: "set", ref: "one", params: { highlights: to } }],
      });
      await flush();
      const ch = new ChannelMarks(ctx);
      const out = runs(ch, 0, "band");
      ch.dispose();
      ctx.dispose();
      const fresh = await opened(doc(to), P);
      const ch2 = new ChannelMarks(fresh.ctx);
      expect(out).toEqual(runs(ch2, 0, "band"));
      ch2.dispose();
      fresh.ctx.dispose();
      return out;
    };
    expect(await patched([{ data: "sites" }], [{ data: "sites", items: [] }])).toEqual([]);
    expect(await patched([{ data: "sites", items: ["s1", "s2"] }], [{ data: "sites", items: ["s1\ns2"] }])).toEqual([]);
  });

  it("a document without them has no set channels", async () => {
    const { ctx } = await opened(doc(null), P);
    expect(sets(ctx)).toEqual([]);
    ctx.dispose();
  });
});

describe("bookmarks", () => {
  it("add inserts the set and its highlights entry, then appends; step focuses and accents; remove drops one", async () => {
    const { ctx } = await opened(doc(null), P);
    const set = () => ctx.graph.assembled.value?.set as never;
    expect(ctx.actions.can("bookmark.add", { section: "one", loci: null })).toBe(false);
    expect(ctx.actions.can("bookmark.step", { section: "one", dir: 1 })).toBe(false);
    expect(ctx.actions.run("bookmark.add", { section: "one", loci: seqLoci(set(), "a", 10, 14), label: "Ten" })).toBe(
      true,
    );
    await flush();
    const node = () =>
      (ctx.document.shown.value?.root.children ?? [])
        .find((c) => c.kind === "data")
        ?.children?.find((c) => c.ref === "one.bookmarks");
    expect(node()?.params?.value).toEqual([{ target: "a/unit", start: 10, stop: 14, id: "b1", label: "Ten" }]);
    expect(sets(ctx)).toEqual([["set:one/0", "band", "one"]]);
    ctx.actions.run("bookmark.add", { section: "one", loci: seqLoci(set(), "a", 20, 22) });
    await flush();
    expect(bookmarkIds(ctx, "one")).toEqual(["b1", "b2"]);
    expect(node()?.params?.value).toHaveLength(2);

    const focus = () => {
      const f = ctx.focus.current.value.get("one");
      return f && [f.start, f.stop];
    };
    expect(ctx.actions.run("bookmark.step", { section: "one", dir: 1 })).toBe(true);
    expect([bookmarkAt(ctx, "one"), focus()]).toEqual(["b1", [10, 14]]);
    ctx.actions.run("bookmark.step", { section: "one", dir: 1 });
    expect([bookmarkAt(ctx, "one"), focus()]).toEqual(["b2", [20, 22]]);
    ctx.actions.run("bookmark.step", { section: "one", dir: 1 });
    expect(bookmarkAt(ctx, "one")).toBe("b1");
    ctx.actions.run("bookmark.step", { section: "one", dir: -1 });
    expect(bookmarkAt(ctx, "one")).toBe("b2");
    const ch = new ChannelMarks(ctx);
    expect([runs(ch, 0, "accent"), runs(ch, 1, "accent")]).toEqual([[[20, 22]], []]);
    // Going to one makes it the current one: focused, accented, and where Next and Previous step from.
    expect(ctx.actions.can("bookmark.go", { section: "one", item: "nope" })).toBe(false);
    expect(ctx.actions.run("bookmark.go", { section: "one", item: "b1" })).toBe(true);
    expect([bookmarkAt(ctx, "one"), focus(), runs(ch, 0, "accent")]).toEqual(["b1", [10, 14], [[10, 14]]]);
    ctx.actions.run("bookmark.step", { section: "one", dir: 1 });
    expect(bookmarkAt(ctx, "one")).toBe("b2");

    expect(ctx.actions.can("bookmark.remove", { section: "one", item: "nope" })).toBe(false);
    expect(ctx.actions.run("bookmark.remove", { section: "one", item: "b1" })).toBe(true);
    await flush();
    expect(bookmarkIds(ctx, "one")).toEqual(["b2"]);
    ch.dispose();
    ctx.dispose();
  });
});
