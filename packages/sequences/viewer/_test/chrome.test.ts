/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The default UIs' shared helpers (seqspec:interaction#the-default-ui, seqspec:interaction#slots-and-chrome; UX1):
// the chrome's compositions (`placedParts`) with a part placed under `clean` by its slot, and the slots
// (`resolveSlots`), which the chrome never changes, `undefined` skipped and an invalid slot reported once; the track menu from `UI.TrackActions` with a host's own actions, `none` and "Retry" only on a failed
// region; the section header, hidden tracks, the plain `Description`, and the track the header menu sits on.

import { createSeqViewSpec } from "@molstar/seq-core/builder";
import type { Loaders, SeqViewSpecDoc } from "@molstar/seq-core/format";
import { OrderedSet } from "@molstar/seq-core/model";
import { NotServed, tilesStream } from "@molstar/seq-core/stream/tiles";
import { tileServer } from "@molstar/seq-core/stream/tiles/server";
import { describe, expect, it } from "vitest";
import { case9, case12, case13 } from "../../core/_test/format/cases.js";
import { files, fixtureLoaders } from "../../core/_test/format/load-fixtures.js";
import {
  headerLine,
  headerTrack,
  hiddenTracks,
  isDomPart,
  isReactPart,
  issueBadge,
  issueList,
  issuesByRef,
  labelText,
  placedParts,
  plainDescription,
  resolveSlots,
  sectionLabel,
  trackLabel,
  trackLine,
  trackMenu,
} from "../src/chrome.js";
import {
  ChromeCompositions,
  CleanParts,
  descriptionOf,
  type KeyCursor,
  SeqViewerConfig,
  type SeqViewerContext,
  sectionInfo,
  shownSections,
} from "../src/index.js";
import { driven } from "./input-fixtures.js";
import { opened, preparedRenderers } from "./layout-fixtures.js";

const flush = async () => {
  for (let i = 0; i < 5; i++) await new Promise((r) => setTimeout(r, 0));
};
const S = "v.main";

/** One sequence, a letters track and a ranges track with its data. */
function tracksDoc(): SeqViewSpecDoc {
  const b = createSeqViewSpec({ title: "Menu" });
  const a = b.data.sequence({ ref: "a", label: "A", value: "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ" });
  const d = b.data.annotation({
    ref: "d",
    geometry: "ranges",
    target: a,
    value: { start: [2, 9], stop: [4, 12], label: ["x", "y"] },
  });
  const v = b.view({ ref: "v", label: "V" }).block({ sequence: a });
  v.track({ ref: "v.letters", label: "Letters" }).letters();
  v.track({ ref: "v.feat", label: "Features" }).ranges({ data: d });
  return b.getDocument();
}
const P = { renderers: preparedRenderers() };
const items = (ctx: SeqViewerContext, track: string) => trackMenu(ctx, S, track).map((m) => m.label);

/** The parts a default UI has, each a function; the host's replace them. */
const DEFAULTS = {
  ViewTabs: () => "tabs",
  SectionHeader: () => "header",
  TrackActions: () => "actions",
  HiddenChips: () => "chips",
  Legend: () => "legend",
  Export: () => "export",
  Issues: () => "issues",
  Tooltip: () => "tooltip",
  SnapshotPlayer: () => "player",
  Description: () => "description",
  Hint: () => "hint",
  ScrollArea: () => "scroll",
  Menu: () => "menu",
};
type Slots = typeof DEFAULTS;
const nulls = (s: Slots) =>
  Object.entries(s)
    .filter(([, v]) => v === null)
    .map(([k]) => k);

const sorted = (s: ReadonlySet<string>) => [...s].sort();
describe("the chrome's compositions and the slots", () => {
  it("full places its whole composition; clean leaves `CleanParts` out and keeps the canvas parts, the player and the states", async () => {
    const { ctx } = await opened(tracksDoc(), P);
    expect(sorted(placedParts(ctx, "full", undefined))).toEqual([...ChromeCompositions.full].sort());
    const clean = placedParts(ctx, "clean", undefined);
    expect(sorted(clean)).toEqual(
      ChromeCompositions.full.filter((k) => !(CleanParts as readonly string[]).includes(k)).sort(),
    );
    expect(["Tooltip", "HoverCard", "SnapshotPlayer", "Description", "Search"].every((k) => clean.has(k))).toBe(true);
    // The chrome never changes a slot: a part a host places itself is the same part under either chrome.
    expect(nulls(resolveSlots(ctx, DEFAULTS, undefined, isReactPart))).toEqual([]);
    ctx.dispose();
  });

  it("the chrome is the UI's option: unset it reads the config (clean by default), and the config is never written", async () => {
    const { ctx } = await opened(tracksDoc(), P);
    expect(ctx.config.get(SeqViewerConfig.UI.Chrome)).toBe("clean");
    expect(placedParts(ctx, undefined, undefined).has("Legend")).toBe(false);
    ctx.config.set(SeqViewerConfig.UI.Chrome, "full");
    expect(placedParts(ctx, undefined, undefined).has("Legend")).toBe(true);
    ctx.config.set(SeqViewerConfig.UI.Chrome, "clean");
    expect(placedParts(ctx, undefined, undefined).has("Legend")).toBe(false);
    // One clean UI does not turn another on the same context clean, nor the reverse.
    expect(placedParts(ctx, "full", undefined).has("Legend")).toBe(true);
    expect(ctx.config.get(SeqViewerConfig.UI.Chrome)).toBe("clean");
    ctx.config.set(SeqViewerConfig.UI.Chrome, "full");
    placedParts(ctx, "clean", undefined);
    expect(ctx.config.get(SeqViewerConfig.UI.Chrome)).toBe("full");
    ctx.dispose();
  });

  it("a host's slot places its part under clean (a part or `default`); `none` leaves it out wherever it is placed", async () => {
    const { ctx } = await opened(tracksDoc(), P);
    const mine = () => "mine";
    const own = placedParts(ctx, "clean", { Legend: mine });
    expect([own.has("Legend"), own.has("ViewTabs")]).toEqual([true, false]);
    expect(resolveSlots(ctx, DEFAULTS, { Legend: mine }, isReactPart).Legend).toBe(mine);
    const back = placedParts(ctx, "clean", { Legend: "default", Issues: "default" });
    expect([back.has("Legend"), back.has("Issues"), back.has("ViewTabs")]).toEqual([true, true, false]);
    // `none`: not placed, and null as a slot (a header leaves its chips out); a primitive keeps its default.
    expect(placedParts(ctx, "full", { Legend: "none" }).has("Legend")).toBe(false);
    const none = resolveSlots(
      ctx,
      DEFAULTS,
      { HiddenChips: "none", Hint: "none", Menu: "none", ScrollArea: "none" },
      isReactPart,
    );
    expect([none.HiddenChips, none.Hint, none.Menu, none.ScrollArea]).toEqual([
      null,
      DEFAULTS.Hint,
      DEFAULTS.Menu,
      DEFAULTS.ScrollArea,
    ]);
    ctx.dispose();
  });

  it("a slot set to `undefined` (a part not loaded yet) is not set: the chrome's choice stands", async () => {
    const { ctx } = await opened(tracksDoc(), P);
    const slots = resolveSlots(ctx, DEFAULTS, { Menu: undefined, Legend: undefined }, isReactPart);
    expect([slots.Menu, slots.Legend]).toEqual([DEFAULTS.Menu, DEFAULTS.Legend]);
    expect(placedParts(ctx, "clean", { Legend: undefined }).has("Legend")).toBe(false);
    ctx.dispose();
  });

  it("a value that is not a part falls back to the default and reports `seqv.ui.slot` once per slot", async () => {
    const { ctx } = await opened(tracksDoc(), P);
    const issues: string[] = [];
    ctx.events.issue.subscribe((i) => issues.push(`${i.code} ${i.message}`));
    const severities = new Set<string>();
    ctx.events.issue.subscribe((i) => severities.add(i.severity));
    const bad = { Menu: 42, Legend: "bogus", Hint: true } as unknown as Partial<Record<keyof Slots, "none">>;
    for (let i = 0; i < 3; i++) {
      const out = resolveSlots(ctx, DEFAULTS, bad, isReactPart);
      expect([out.Menu, out.Legend, out.Hint]).toEqual([DEFAULTS.Menu, DEFAULTS.Legend, DEFAULTS.Hint]);
    }
    expect(issues.map((m) => m.split(" ")[0])).toEqual(["seqv.ui.slot", "seqv.ui.slot", "seqv.ui.slot"]);
    expect(issues.map((m) => /"(\w+)"/.exec(m)?.[1])).toEqual(["Menu", "Legend", "Hint"]);
    // At its registered severity: the default is used, a warning (GV-05; `report` decides no severity of its own).
    expect([...severities]).toEqual(["warning"]);
    // Another context reports on its own; `memo` and `forwardRef` parts (objects) are parts.
    const other = await opened(tracksDoc(), P);
    const memo = { $$typeof: Symbol.for("react.memo") };
    const out = resolveSlots(other.ctx, DEFAULTS, { Menu: memo as never }, isReactPart);
    expect(out.Menu).toBe(memo);
    // The DOM calls a part: an object there is no part, its default is used (VB-12).
    expect(resolveSlots(other.ctx, DEFAULTS, { Legend: {} as never }, isDomPart).Legend).toBe(DEFAULTS.Legend);
    ctx.dispose();
    other.ctx.dispose();
  });

  it("`UI.TrackActions: none` and a hidden legend leave their parts out of the composition; a slot places them again", async () => {
    const { ctx } = await opened(tracksDoc(), P);
    ctx.config.set(SeqViewerConfig.UI.TrackActions, "none");
    expect(placedParts(ctx, "full", undefined).has("TrackActions")).toBe(false);
    expect(placedParts(ctx, "full", undefined, "hidden").has("Legend")).toBe(false);
    expect(placedParts(ctx, "full", { Legend: "default" }, "hidden").has("Legend")).toBe(true);
    ctx.dispose();
  });
});

describe("trackMenu: UI.TrackActions in order, a host's own actions", () => {
  it("the default list: each built-in when it applies, in order", async () => {
    const { ctx } = await opened(tracksDoc(), P);
    expect(items(ctx, "v.feat")).toEqual(["Show only what this track covers", "Focus", "Hide track"]);
    // A track with no data has no filter and no focus.
    expect(items(ctx, "v.letters")).toEqual(["Hide track"]);
    expect(trackMenu(ctx, S, "v.feat").map((m) => [m.key, !!m.danger])).toEqual([
      ["track.filterView", false],
      ["focus", false],
      ["track.hide", false],
    ]);
    // "Remove track" is built in but off by default: a host adds `remove` to the list.
    ctx.config.set(SeqViewerConfig.UI.TrackActions, ["hide", "remove"]);
    expect(trackMenu(ctx, S, "v.feat").map((m) => [m.key, !!m.danger])).toEqual([
      ["track.hide", false],
      ["track.remove", true],
    ]);
    ctx.dispose();
  });

  it("Focus moves the section's viewport onto what the track covers", async () => {
    const { ctx, scheduler } = await opened(tracksDoc(), P);
    const vp = ctx.viewports.get(S);
    if (!vp) throw new Error("no viewport");
    vp.setSize(800, 200, 1);
    await ctx.graph.update();
    const logs: string[] = [];
    ctx.events.issue.subscribe((i) => logs.push(`${i.code} ${i.message}`));
    const before = vp.live.current.cell;
    trackMenu(ctx, S, "v.feat")
      .find((m) => m.key === "focus")
      ?.run();
    // The menu closes and gives the surface its focus back, as the DOM does.
    ctx.input.focus(ctx.layout.sectionIndex(S), true);
    for (let i = 0; i < 5; i++) await new Promise((r) => setTimeout(r, 0));
    scheduler.advance(2000);
    await ctx.graph.update();
    for (let i = 0; i < 5; i++) await new Promise((r) => setTimeout(r, 0));
    const s = vp.live.current;
    // The features cover positions 2..12 of 33: the view zooms in to them.
    expect({ zoomed: s.cell > before, logs }).toEqual({ zoomed: true, logs: [] });
    ctx.dispose();
  });

  it("Focus covers what the track's filter keeps, not all of the data its layer reads", async () => {
    const b = createSeqViewSpec({ title: "Filtered" });
    const a = b.data.sequence({ ref: "a", label: "A", value: "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ" });
    const d = b.data.annotation({
      ref: "d",
      geometry: "ranges",
      target: a,
      value: { start: [0, 9, 20], stop: [33, 12, 24], category: ["Domain", "Motif", "Motif"] },
    });
    const v = b.view({ ref: "v", label: "V" }).block({ sequence: a });
    v.track({ ref: "v.motif", label: "Motif" }).ranges({ data: d, filter: { field: "category", values: ["Motif"] } });
    v.track({ ref: "v.all", label: "All" }).ranges({ data: d });
    const { ctx } = await opened(b.getDocument(), P);
    const vp = ctx.viewports.get(S);
    if (!vp) throw new Error("no viewport");
    vp.setSize(800, 200, 1);
    await ctx.graph.update();
    const focused = async (track: string) => {
      ctx.focus.set(S, null);
      trackMenu(ctx, S, track)
        .find((m) => m.key === "focus")
        ?.run();
      for (let i = 0; i < 5; i++) await new Promise((r) => setTimeout(r, 0));
      const f = ctx.focus.current.value.get(S);
      return f && [f.block, f.start, f.stop];
    };
    expect([await focused("v.motif"), await focused("v.all")]).toEqual([
      [0, 9, 24],
      [0, 0, 33],
    ]);
    ctx.dispose();
  });

  it("orders the menu as the config does, with a host's actions; a `when` that fails leaves one out", async () => {
    const { ctx } = await opened(tracksDoc(), P);
    const ran: string[] = [];
    ctx.config.set(SeqViewerConfig.UI.TrackActions, [
      "hide",
      "rows",
      { id: "log", label: "Log it", run: (_c, section, track) => void ran.push(`${section} ${track}`) },
      { id: "never", label: "Never", when: () => false, run: () => {} },
      { id: "danger", label: "Drop", danger: true, when: (_c, _s, t) => t === "v.feat", run: () => {} },
    ]);
    expect(items(ctx, "v.feat")).toEqual(["Hide track", "Log it", "Drop"]);
    expect(items(ctx, "v.letters")).toEqual(["Hide track", "Log it"]);
    trackMenu(ctx, S, "v.feat")[1].run();
    expect(ran).toEqual(["v.main v.feat"]);
    expect(trackMenu(ctx, S, "v.feat").map((m) => m.danger)).toEqual([undefined, undefined, true]);
    ctx.dispose();
  });

  it('"none", or a list with nothing that applies, gives no items', async () => {
    const { ctx } = await opened(tracksDoc(), P);
    ctx.config.set(SeqViewerConfig.UI.TrackActions, "none");
    expect(trackMenu(ctx, S, "v.feat")).toEqual([]);
    ctx.config.set(SeqViewerConfig.UI.TrackActions, [
      { id: "never", label: "Never", when: () => false, run: () => {} },
    ]);
    expect(trackMenu(ctx, S, "v.feat")).toEqual([]);
    ctx.config.set(SeqViewerConfig.UI.TrackActions, ["rows", "retry"]);
    expect(trackMenu(ctx, S, "v.feat")).toEqual([]);
    ctx.dispose();
  });

  it("running an item proposes the action's patch: Hide hides the track, and `hiddenTracks` lists it", async () => {
    const { ctx } = await opened(tracksDoc(), P);
    expect(hiddenTracks(ctx, S)).toEqual([]);
    trackMenu(ctx, S, "v.feat")
      .find((m) => m.key === "track.hide")
      ?.run();
    await flush();
    expect(hiddenTracks(ctx, S)).toEqual([{ ref: "v.feat", label: "Features" }]);
    expect(ctx.actions.run("track.show", "v.feat")).toBe(true);
    await flush();
    expect(hiddenTracks(ctx, S)).toEqual([]);
    ctx.dispose();
  });

  it("`rows`: the toggle of a split track, named by its state; the group action and the menu agree (ported)", async () => {
    const t = await opened(case13(), { size: { "features.main": [1000, 900] }, renderers: preparedRenderers() });
    const NV = "features.natural-variants";
    const menu = () => trackMenu(t.ctx, "features.main", NV).find((m) => m.key === "rows");
    const rows = () => t.entries("features.main").filter(([, name]) => name === NV).length;
    expect(menu()?.label).toBe("Collapse rows");
    const open = rows();
    menu()?.run();
    await flush();
    expect([menu()?.label, rows()]).toEqual(["Expand rows", 1]);
    expect(t.ctx.actions.run("group.toggle", NV)).toBe(true);
    await flush();
    expect([menu()?.label, rows()]).toEqual(["Collapse rows", open]);
    // Only a split track has the toggle.
    expect(t.ctx.actions.can("group.toggle", "features.main")).toBe(false);
    t.ctx.dispose();
  });
});

describe('"Retry" in the track menu', () => {
  const URI = "https://example.org/lazy";
  function lazy() {
    const b = createSeqViewSpec({ title: "lazy" });
    const aln = b.data.alignment({
      ref: "lazy",
      row_count: 2000,
      column_count: 400,
      source: { uri: URI, protocol: "seqstar-tiles/1" },
    });
    const v = b.view({ ref: "v", label: "A", style: "bare", height: "box" });
    v.block({ alignment: aln });
    v.track({ ref: "v.rows", label: "Rows" }).rows();
    return b.getDocument();
  }
  /** Tiles fail `n` times (NotServed when asked), then answer. */
  function loaders(n: number, notServed = false): Loaders {
    const server = tileServer({ rowCount: 2000, columnCount: 400 });
    const seen = new Map<string, number>();
    return {
      load: () => Promise.reject(new Error("no fixture")),
      stream: tilesStream((u) =>
        u === URI
          ? (path, signal) => {
              const k = seen.get(path) ?? 0;
              seen.set(path, k + 1);
              return path === "tiles/0/0/0/0" && k < n
                ? Promise.reject(notServed ? new NotServed("404") : new Error("down"))
                : server.transport(path, signal);
            }
          : null,
      ),
    };
  }
  const pass = (ctx: SeqViewerContext) => (ctx.streams as unknown as { pass(): void }).pass();

  it("is in the menu only while a layer of the track has a failed region, and ends the cooldown", async () => {
    const t = await opened(lazy(), {
      ...{ renderers: preparedRenderers() },
      loaders: loaders(4),
      size: { "v.main": [800, 400] },
    });
    const rows = () => trackMenu(t.ctx, S, "v.rows").map((m) => m.label);
    pass(t.ctx);
    await flush();
    // The first failure: backing off, and the menu offers Retry.
    expect(rows()).toEqual(["Hide track", "Retry"]);
    expect(t.ctx.actions.can("stream.retry", "v.rows")).toBe(true);
    for (const ms of [1000, 2000, 4000]) {
      t.scheduler.advance(ms);
      await flush();
    }
    const retry = trackMenu(t.ctx, S, "v.rows").find((m) => m.key === "stream.retry");
    expect(retry?.label).toBe("Retry");
    retry?.run();
    t.scheduler.advance(16);
    t.scheduler.advance(16);
    await flush();
    expect([
      trackMenu(t.ctx, S, "v.rows").some((m) => m.key === "stream.retry"),
      t.ctx.actions.can("stream.retry", "v.rows"),
    ]).toEqual([false, false]);
    t.ctx.dispose();
  });

  it("not for a region the source does not serve (it is never retried), nor for another track", async () => {
    const t = await opened(lazy(), {
      renderers: preparedRenderers(),
      loaders: loaders(9, true),
      size: { "v.main": [800, 400] },
    });
    pass(t.ctx);
    await flush();
    expect(trackMenu(t.ctx, S, "v.rows").some((m) => m.key === "stream.retry")).toBe(false);
    expect(t.ctx.actions.can("stream.retry", "nothing")).toBe(false);
    t.ctx.dispose();
  });
});

describe("the section header, the description and the header track", () => {
  it("sectionInfo: a switch whose value removes what offers it stays offered (case 12's heme as atoms, its copies each)", async () => {
    const { ctx } = await opened(case12(), P);
    const S12 = "design-view.main";
    const value = (action: string) => sectionInfo(ctx, S12).switches.find((x) => x.action === action)?.value;
    expect(value("section.ligandForm")).toBe("position");
    for (const [action, v] of [
      ["section.ligandForm", "atoms"],
      ["section.ligandForm", "position"],
      ["section.copies", "each"],
      ["section.copies", "one"],
    ] as const) {
      await ctx.actions.run(action, { section: S12, value: v });
      await ctx.graph.update();
      expect([action, value(action)]).toEqual([action, v]);
    }
    ctx.dispose();
  });

  it("sectionInfo: the numbering switch offers only numberings that label differently (case 9's 1A3N: none)", async () => {
    const numbering = (ctx: SeqViewerContext) =>
      sectionInfo(ctx, "msa.chains")
        .switches.find((x) => x.action === "section.numbering")
        ?.options.map((o) => o.value);
    // 1A3N's chains: auth (1-141, 1-146), label (label_seq_id from 1) and ordinal label every residue alike.
    const plain = await opened(case9(), P);
    expect(numbering(plain.ctx)).toBeUndefined();
    // The section naming one keeps it reachable, and the others it equals.
    await plain.ctx.actions.run("section.numbering", { section: "msa.chains", value: "label" });
    await plain.ctx.graph.update();
    expect(numbering(plain.ctx)).toEqual([null, "auth", "label"]);
    plain.ctx.dispose();
    // Chain B's auth numbers from 101: auth differs; label still equals ordinal.
    const f = files();
    const uri = "https://www.ebi.ac.uk/pdbe/entry-files/download/1a3n_updated.cif";
    const rec = JSON.parse(f[uri]);
    for (const u of rec.units) if (u.key === "B") u.numberings[0].segments = [[0, 101, 146]];
    f[uri] = JSON.stringify(rec);
    const offset = await opened(case9(), { ...P, loaders: fixtureLoaders(f) });
    expect(numbering(offset.ctx)).toEqual([null, "auth", "ordinal"]);
    offset.ctx.dispose();
  });

  it("sectionInfo: the label, the style switch over the presets of the section's geometry, a clear-filter chip", async () => {
    const { ctx } = await opened(tracksDoc(), P);
    const info = sectionInfo(ctx, S);
    expect(info.label).toBe(null);
    const style = info.switches.find((s) => s.action === "section.style");
    expect(style?.options.map((o) => o.label)).toEqual(["auto", "tracks", "msa", "msa-fixed", "bare", "railed"]);
    expect(info.filter).toBe(null);
    ctx.actions.run("track.filterView", "v.feat");
    await flush();
    expect(sectionInfo(ctx, S).filter).toBe("d");
    // The style switch is a runtime action (`ui.style`), not a patch.
    expect(ctx.actions.run("section.style", { section: S, value: "railed" })).toBe(true);
    expect(sectionInfo(ctx, S).switches.find((s) => s.action === "section.style")?.value).toBe("railed");
    expect(ctx.actions.run("section.style", { section: S, value: null })).toBe(true);
    expect(sectionInfo(ctx, S).switches.find((s) => s.action === "section.style")?.value).toBe(null);
    ctx.dispose();
  });

  it('"Show only what this track covers" on a filtered track keeps the rows its filter keeps, not all of the data', async () => {
    const b = createSeqViewSpec({ title: "Filtered" });
    const a = b.data.sequence({ ref: "a", label: "A", value: "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ" });
    const d = b.data.annotation({
      ref: "d",
      geometry: "ranges",
      target: a,
      value: { start: [0, 9, 20], stop: [33, 12, 24], category: ["Domain", "Motif", "Motif"] },
    });
    const v = b.view({ ref: "v", label: "V" }).block({ sequence: a });
    v.track({ ref: "v.motif", label: "Motif" }).ranges({ data: d, filter: { field: "category", values: ["Motif"] } });
    v.track({ ref: "v.all", label: "All" }).ranges({ data: d });
    const { ctx } = await opened(b.getDocument(), P);
    const shown = async (track: string) => {
      ctx.actions.run("section.clearFilter", S);
      await flush();
      ctx.actions.run("track.filterView", track);
      await flush();
      const set = ctx.graph.assembled.value?.view?.resolved.blocks.positions[0];
      return set ? Array.from({ length: OrderedSet.size(set) }, (_, i) => OrderedSet.getAt(set, i)) : null;
    };
    expect(await shown("v.motif")).toEqual([9, 10, 11, 20, 21, 22, 23]);
    expect(sectionInfo(ctx, S).filter).toBe("Motif");
    expect(await shown("v.all")).toHaveLength(33);
    ctx.dispose();
  });

  it("issuesByRef groups neighbors of one ref in the manager's order; the document's own is null", () => {
    const i = (ref: string | null, code: string) => ({
      code,
      severity: "warning" as const,
      path: "",
      message: code,
      ref,
    });
    const groups = issuesByRef([i("a", "x"), i("a", "y"), i("b", "z"), i(null, "w")]);
    expect(groups.map((g) => [g.ref, g.issues.map((x) => x.code)])).toEqual([
      ["a", ["x", "y"]],
      ["b", ["z"]],
      [null, ["w"]],
    ]);
  });

  it("the issue badge counts errors and warnings, marked by the most severe; notes alone are a quiet count (VB-14)", () => {
    const i = (severity: "error" | "warning" | "info") => ({ code: "x", severity, path: "", message: "", ref: null });
    expect(issueBadge([])).toBeNull();
    expect(issueBadge([i("info"), i("info")])).toEqual({ severity: "info", text: "2 notes" });
    expect(issueBadge([i("info"), i("warning")])).toEqual({ severity: "warning", text: "1 issue" });
    expect(issueBadge([i("warning"), i("error"), i("info")])).toEqual({ severity: "error", text: "2 issues" });
  });

  it("the issue list cuts by severity, keeps the manager's order within it, and lists notes last (VB-14)", () => {
    const i = (ref: string | null, severity: "error" | "warning" | "info", code: string) => ({
      code,
      severity,
      path: "",
      message: code,
      ref,
    });
    // In the manager's order (by ref): a's notes and warnings, then z's error, after the cut's worth of warnings.
    const many = [i("a", "info", "n1"), i("a", "warning", "w1"), i("a", "warning", "w2"), i("z", "error", "e1")];
    const cut = issueList(many, 2);
    expect([cut.problems.map((g) => [g.ref, g.issues.map((x) => x.code)]), cut.notes, cut.more]).toEqual([
      [
        ["a", ["w1"]],
        ["z", ["e1"]],
      ],
      [],
      2,
    ]);
    const all = issueList(many);
    expect([
      all.problems.map((g) => g.issues.map((x) => x.code)),
      all.notes.map((g) => [g.ref, g.issues.length]),
    ]).toEqual([[["w1", "w2"], ["e1"]], [["a", 1]]]);
    expect(all.more).toBe(0);
  });

  it("labels: tracks, sections and the shown sections of the current view", async () => {
    const { ctx } = await opened(tracksDoc(), P);
    expect([trackLabel(ctx, "v.feat"), sectionLabel(ctx, S), shownSections(ctx)]).toEqual(["Features", "V", [S]]);
    expect(labelText({ title: "A T315", lines: [{ key: "uniprot", value: "334" }] })).toBe("A T315, uniprot 334");
    expect(labelText(null)).toBe("");
    ctx.dispose();
  });

  it("the default Description shows a `plaintext` description as paragraphs and a `markdown` one not at all", async () => {
    const doc = tracksDoc();
    doc.metadata = {
      description: "First line\nstill first.\n\nSecond [link](https://x.test) *plain*.",
      description_format: "plaintext",
    };
    const { ctx } = await opened(doc, P);
    expect(plainDescription(ctx)).toEqual(["First line still first.", "Second [link](https://x.test) *plain*."]);
    expect(descriptionOf(ctx)?.format).toBe("plaintext");
    ctx.dispose();
    const md = tracksDoc();
    md.metadata = {
      description: "# Title\n\nSome **bold** [T315](ctx:site.gatekeeper).",
      description_format: "markdown",
    };
    const m = await opened(md, P);
    expect([descriptionOf(m.ctx)?.format, plainDescription(m.ctx)]).toEqual(["markdown", null]);
    m.ctx.dispose();
    const none = await opened(tracksDoc(), P);
    expect([descriptionOf(none.ctx), plainDescription(none.ctx)]).toEqual([null, null]);
    none.ctx.dispose();
  });

  it("the snapshot's description wins over the document's", async () => {
    const doc = tracksDoc();
    doc.metadata = { description: "The document.", description_format: "plaintext" };
    doc.snapshots = [{ ref: "s1", ops: [], description: "The snapshot.", description_format: "plaintext" }];
    const { ctx } = await opened(doc, P);
    expect(plainDescription(ctx)).toEqual(["The document."]);
    await ctx.snapshots.apply("s1");
    expect(plainDescription(ctx)).toEqual(["The snapshot."]);
    ctx.dispose();
  });

  it("headerTrack: the keyboard cursor's track in the section, else none; a cursor of another section is not this one's", async () => {
    const { ctx } = await opened(tracksDoc(), P);
    // Without the keyboard behavior, whose cursor hovers the header it stands on.
    for (const b of ctx.behaviors.list.value) ctx.behaviors.get(b.name)?.remove();
    expect(headerTrack(ctx, S)).toBe(null);
    const cursor = (section: string, track: string | null): KeyCursor => ({
      section,
      track,
      header: true,
      block: 0,
      index: 0,
      row: -1,
      displayRow: 0,
    });
    ctx.ui.setKeyCursor(cursor(S, "v.feat"));
    expect(headerTrack(ctx, S)).toBe("v.feat");
    ctx.ui.setKeyCursor(cursor("other.main", "v.feat"));
    expect(headerTrack(ctx, S)).toBe(null);
    ctx.ui.setKeyCursor(cursor(S, null));
    expect(headerTrack(ctx, S)).toBe(null);
    ctx.dispose();
  });
});

describe("trackLine: the header line the menu sits on (ported; allocation-free)", () => {
  it("the line under the pointer, else the track's first in view; none for a track not placed", async () => {
    const t = await driven(tracksDoc(), { size: { [S]: [1000, 400] } });
    const out = headerLine();
    const feat = t.yOf(S, "v.feat");
    const letters = t.yOf(S, "v.letters");
    // The pointer elsewhere: the track's own entry, at its place.
    expect(trackLine(t.ctx, S, "v.feat", out)).toBe(true);
    const l = t.layout(S);
    expect([l.entryTrack[out.entry], out.row]).toEqual([t.trackIndex("v.feat"), 0]);
    expect([out.top + out.height / 2]).toEqual([feat]);
    // The pointer over the other track's header does not move this track's line.
    t.pointer("move", S, 10, letters);
    expect(trackLine(t.ctx, S, "v.feat", out)).toBe(true);
    expect(out.top + out.height / 2).toBe(feat);
    t.pointer("move", S, 10, feat + 1);
    expect(trackLine(t.ctx, S, "v.feat", out)).toBe(true);
    expect(trackLine(t.ctx, S, "v.nothing", out)).toBe(false);
    expect(out.entry).toBe(-1);
    t.ctx.dispose();
  });
});
