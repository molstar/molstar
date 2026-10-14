/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A header's hover, as the `headerHover` behavior describes it (seqspec:interaction#hover): a feature track's header its
// shown items (its filter kept), a continuous track's (letters, values) the row band alone, a split row's header that
// row's items, a group's header its shown tracks' (a hidden one's not, a collapsed group's yes); `hover.toLoci()` the
// items as item loci, resolved when read and kept for the version, `hover.size()` known before; past `maxItems` the row
// band alone and nothing forwarded; a filter patched under a resting header, by the pointer or the keyboard, describes
// it again; and a viewer spec without the behavior lights nothing on a header.

import { createSeqViewSpec } from "@molstar/seq-core/builder";
import { OrderedSet } from "@molstar/seq-core/model";
import { describe, expect, it } from "vitest";
import { case13 } from "../../core/_test/format/cases.ts";
import { DefaultHeaderHoverParams, describeHeader, headerLoci, headerSet, lociSize } from "../src/header-hover.ts";
import {
  createSeqViewer,
  DefaultSeqViewerSpec,
  EntryKind,
  type HeaderHoverParams,
  Hit,
  HoverKind,
  SeqViewerBehaviors,
  type SeqViewerContext,
} from "../src/index.ts";
import { opened, preparedRenderers } from "./layout-fixtures.ts";

const S = "v.main";
/** The header's lit runs as "start-stop", sorted. */
const runs = (ctx: SeqViewerContext) => {
  const c = ctx.hover.live.current.header.extents;
  const out = new Set<string>();
  for (let i = 0; i < c.count; i++) out.add(`${c.start[i]}-${c.stop[i]}`);
  return [...out].sort((a, b) => Number.parseInt(a, 10) - Number.parseInt(b, 10));
};
/** `toLoci()`'s items, as "annotation:item". */
const items = (ctx: SeqViewerContext) => {
  const l = ctx.hover.toLoci();
  if (l?.kind !== "item-loci") return l === null ? null : l.kind;
  return l.elements.flatMap((e) =>
    Array.from({ length: OrderedSet.size(e.items) }, (_, k) => `${e.annotation}:${OrderedSet.getAt(e.items, k)}`),
  );
};
/** The pointer's (or the keyboard's) hover at (x, y), as the hover behavior sets it. */
const hoverAt = (ctx: SeqViewerContext, x: number, y: number, owner = "pointer") => {
  const hit = Hit.create();
  ctx.layout.hitTest(S, x, y, hit);
  ctx.hover.setFromHit(hit, owner);
  return hit;
};

/** Feature tracks in a group (one in an inner group, one hidden), one outside it, one filtered, letters and values. */
function doc() {
  const b = createSeqViewSpec({ title: "Headers" });
  const a = b.data.sequence({ ref: "a", value: "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ" });
  const ranges = (ref: string, start: number[], stop: number[], label: string[]) =>
    b.data.annotation({ ref, geometry: "ranges", target: a, value: { start, stop, label } });
  const d1 = ranges("d1", [1, 6], [3, 8], ["p", "q"]);
  const d2 = ranges("d2", [12], [15], ["r"]);
  const d3 = ranges("d3", [20, 25], [22, 28], ["s", "t"]);
  const d4 = ranges("d4", [30], [32], ["u"]);
  const vals = b.data.annotation({ ref: "vals", geometry: "values", target: a, value: { value: [0.5, 1, 0.2] } });
  const v = b.view({ ref: "v", label: "V" }).block({ sequence: a });
  v.track({ ref: "v.letters", label: "Letters" }).letters();
  v.track({ ref: "v.vals", label: "Values" }).values({ data: vals } as never);
  v.track({ ref: "v.outside", label: "Outside" }).ranges({ data: d3 });
  v.track({ ref: "v.only", label: "Only t" }).ranges({ data: d3, filter: { field: "label", values: ["t"] } } as never);
  const g = v.group({ ref: "v.g", label: "Group", collapsed: false } as never);
  g.track({ ref: "v.g.one", label: "One" }).ranges({ data: d1 });
  g.track({ ref: "v.g.hidden", label: "Hidden", hidden: true } as never).ranges({ data: d4 });
  g.group({ ref: "v.g.inner", label: "Inner" } as never)
    .track({ ref: "v.g.two", label: "Two" })
    .ranges({ data: d2 });
  return b.getDocument();
}
const open = () => opened(doc(), { size: { [S]: [1000, 700] } });

describe("a header's hover (the headerHover behavior)", () => {
  it("a feature track's header: its shown items, its filter kept; item loci when read, the size before", async () => {
    const t = await open();
    const hit = hoverAt(t.ctx, 10, t.yOf(S, "v.outside"));
    expect(hit.kind).toBe(HoverKind.header);
    const h = t.ctx.hover.live.current.header;
    expect([h.on, h.items, h.length, h.light, h.forward]).toEqual([true, 2, 5, true, true]);
    expect(t.ctx.hover.size()).toEqual({ items: 2, length: 5 });
    expect(runs(t.ctx)).toEqual(["20-22", "25-28"]);
    const d3 = t.ctx.graph.assembled.value?.set.annotations.find((x) => x && x.count === 2 && x.start[0] === 20)?.id;
    expect(items(t.ctx)).toEqual([`${d3}:0`, `${d3}:1`]);
    // Resolved once a version: the same loci again until the hover changes.
    expect(t.ctx.hover.toLoci()).toBe(t.ctx.hover.toLoci());
    hoverAt(t.ctx, 10, t.yOf(S, "v.only"));
    expect(runs(t.ctx)).toEqual(["25-28"]);
    expect(items(t.ctx)).toEqual([`${d3}:1`]);
  });

  it("a continuous track's header (letters, values): the row band alone, nothing positional, nothing forwarded", async () => {
    const t = await open();
    for (const track of ["v.letters", "v.vals"]) {
      hoverAt(t.ctx, 10, t.yOf(S, track));
      const h = t.ctx.hover.live.current.header;
      expect([track, h.items, h.light, h.extents.count, t.ctx.hover.toLoci(), t.ctx.hover.size()]).toEqual([
        track,
        0,
        false,
        0,
        null,
        null,
      ]);
      expect(t.ctx.hover.live.current.track).toBeGreaterThanOrEqual(0);
    }
  });

  it("a group's header: its shown tracks, an inner group's too, a hidden track's not", async () => {
    const t = await open();
    const l = t.layout(S);
    const view = t.ctx.graph.assembled.value?.view;
    let y = -1;
    for (let e = 0; e < l.entryCount; e++)
      if (
        l.entryKind[e] === EntryKind.group &&
        view?.set.strings.get(view.resolved.groups.ref[l.entryGroup[e]]) === "v.g"
      )
        y = l.scrolledTop + l.entryTop[e] + l.entryHeight[e] / 2;
    const hit = hoverAt(t.ctx, 10, y);
    expect([hit.kind, hit.group >= 0]).toEqual([HoverKind.header, true]);
    expect(runs(t.ctx)).toEqual(["1-3", "6-8", "12-15"]);
  });

  it("a group whose tracks share an annotation counts its items and positions once, as its loci resolve (VS-10)", async () => {
    const b = createSeqViewSpec({ title: "G" });
    const a = b.data.sequence({ ref: "a", value: "MKCAYCAKQRQISFVKSHFSRQLEERLGLIEVQ" });
    const sites = { start: [2, 9], stop: [4, 12], id: ["s1", "s2"] };
    b.data.annotation({ ref: "sites", geometry: "ranges", target: a, value: sites });
    b.data.annotation({ ref: "other", geometry: "ranges", target: a, value: { start: [20], stop: [22], id: ["o1"] } });
    const sec = b.view({ ref: "v", label: "V" }).section({ ref: "one" } as never);
    sec.block({ sequence: a });
    const g = sec.group({ ref: "g", label: "G" } as never);
    g.track({ ref: "t1", label: "T1" } as never).ranges({ data: "sites" } as never);
    g.track({ ref: "t2", label: "T2" } as never).ranges({ data: "other" } as never);
    g.track({ ref: "t3", label: "T3" } as never).ranges({ data: "sites" } as never);
    const { ctx } = await opened(b.getDocument(), { renderers: preparedRenderers() });
    const model = ctx.graph.assembled.value?.view;
    const out = headerSet();
    out.group = model ? model.resolved.groups.ref.indexOf(model.set.strings.find("g")) : -1;
    describeHeader(model, DefaultHeaderHoverParams, out);
    const size = lociSize(headerLoci(model, DefaultHeaderHoverParams, out));
    expect({ items: out.items, length: out.length }).toEqual(size);
    expect(size).toEqual({ items: 3, length: 7 });
    // Described again: the same count (what the last pass saw is cleared).
    describeHeader(model, DefaultHeaderHoverParams, out);
    expect([out.items, out.length]).toEqual([3, 7]);
    ctx.dispose();
  });

  it("a split track's row header: that row's items; its label's header the whole track", async () => {
    const t = await opened(case13(), { size: { "features.main": [1000, 900] }, renderers: preparedRenderers() });
    const NV = "features.natural-variants";
    const ROWS = [..."GAVLISTCMDNEQRKHFYWP", "-", "*"];
    const at = (row: number) => {
      const hit = Hit.create();
      t.ctx.layout.hitTest("features.main", 10, t.yOf("features.main", NV, row), hit);
      t.ctx.hover.setFromHit(hit, "pointer");
      return runs(t.ctx).map((r) => Number.parseInt(r, 10));
    };
    expect(at(ROWS.indexOf("H"))).toEqual([174, 272]);
    expect(at(-2).length).toBeGreaterThan(2);
  });

  it("past `maxItems`: the row band alone and nothing forwarded; the count stays on the hover", async () => {
    const t = await open();
    t.ctx.behaviors.get<HeaderHoverParams>("headerHover")?.setParams({ maxItems: 1 });
    hoverAt(t.ctx, 10, t.yOf(S, "v.outside"));
    const h = t.ctx.hover.live.current.header;
    expect([h.items, h.light, h.forward, h.extents.count]).toEqual([2, false, false, 0]);
    expect([t.ctx.hover.toLoci(), t.ctx.hover.size()]).toEqual([null, null]);
    // `forward: false`: lit, but no linked view gets it; a type set to the row band lights nothing.
    t.ctx.behaviors.get<HeaderHoverParams>("headerHover")?.setParams({ maxItems: 5_000, forward: false });
    hoverAt(t.ctx, 10, t.yOf(S, "v.only"));
    expect([runs(t.ctx), t.ctx.hover.toLoci()]).toEqual([["25-28"], null]);
    t.ctx.behaviors.get<HeaderHoverParams>("headerHover")?.setParams({ forward: true, types: { ranges: "rowBand" } });
    hoverAt(t.ctx, 10, t.yOf(S, "v.outside"));
    expect([t.ctx.hover.live.current.header.items, t.ctx.hover.toLoci()]).toEqual([0, null]);
  });

  it("a filter patched under a resting header is read again, for the pointer's hover and the keyboard's", async () => {
    const t = await open();
    const view = t.ctx.graph.assembled.value?.view;
    const rv = view?.resolved;
    let layer = "";
    for (let i = 0; rv && view && i < rv.layers.count; i++)
      if (view.set.strings.get(rv.tracks.ref[rv.layers.track[i]]) === "v.only")
        layer = view.set.strings.get(rv.layers.ref[i]);
    for (const owner of ["pointer", "keyboard"]) {
      hoverAt(t.ctx, 10, t.yOf(S, "v.only"), owner);
      expect([owner, runs(t.ctx)]).toEqual([owner, ["25-28"]]);
      const version = t.ctx.hover.live.version.value;
      await t.ctx.applyPatch({
        format: "seq_view_spec_patch",
        ops: [{ op: "set", ref: layer, params: { filter: { field: "label", values: ["s"] } } }],
      });
      await t.ctx.graph.update();
      expect([owner, runs(t.ctx), t.ctx.hover.live.version.value > version]).toEqual([owner, ["20-22"], true]);
      await t.ctx.applyPatch({
        format: "seq_view_spec_patch",
        ops: [{ op: "set", ref: layer, params: { filter: { field: "label", values: ["t"] } } }],
      });
      await t.ctx.graph.update();
    }
  });

  it("a viewer spec without headerHover lights nothing on a header; removing the behavior turns it off", async () => {
    const spec = DefaultSeqViewerSpec();
    const ctx = createSeqViewer({
      ...spec,
      behaviors: spec.behaviors.filter((b) => b.behavior !== SeqViewerBehaviors.HeaderHover),
    });
    expect(ctx.behaviors.get("headerHover")).toBeNull();
    ctx.dispose();
    const t = await open();
    t.ctx.behaviors.get("headerHover")?.remove();
    hoverAt(t.ctx, 10, t.yOf(S, "v.outside"));
    const h = t.ctx.hover.live.current.header;
    expect([t.ctx.hover.live.current.kind, h.on, h.extents.count, t.ctx.hover.toLoci()]).toEqual([
      HoverKind.header,
      false,
      0,
      null,
    ]);
  });
});
