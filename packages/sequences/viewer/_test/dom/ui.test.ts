/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The vanilla UI (UX2; seqspec:interaction#slots-and-chrome, seqspec:interaction#accessibility,
// seqspec:interaction#the-dom-element): the element connects, disconnects, moves and reconnects (a context it made
// outlives a move and is disposed a task after a disconnect, a host's is not), each member of the table changes only
// what it names, its events are frozen copies; every part mounts and disposes with no subscription or node left; a
// state change updates only the nodes it owns (a mutation count); slots replace parts; menus close on Escape, on a
// choice and on an outside press and give focus back; a rebuilt node keeps focus by key; the default Description shows
// plain text only. React/DOM parity is _test/react/parity.test.tsx (UX3).

import { afterEach, describe, expect, it, vi } from "vitest";
import * as D from "../../src/dom/index.ts";
// The skin, for the legend's grid (its layout is CSS).
import "../../src/skin.css";
import {
  type ActionEntry,
  createSeqViewer,
  Hit,
  type ItemTarget,
  manualScheduler,
  SeqViewerConfig,
  type SeqViewerContext,
} from "../../src/index.ts";
import { fixtureMeasurer } from "../../src/testing/index.ts";
import {
  at,
  flush,
  headlessSpec,
  keys,
  lines,
  openedContext,
  pointer,
  stubSurface,
  subscriptions,
  uiDoc,
  WIDTH,
  xOf,
} from "./fixtures.ts";

D.defineSeqViewer();
const cleanups: (() => void)[] = [];
afterEach(() => {
  for (const c of cleanups.splice(0)) c();
});

function container() {
  const el = document.createElement("div");
  el.style.width = `${WIDTH}px`;
  document.body.appendChild(el);
  cleanups.push(() => el.remove());
  return el;
}
async function opened() {
  const t = await openedContext();
  cleanups.push(() => t.ctx.dispose());
  return t;
}
async function mounted(spec: D.DomUISpec = {}) {
  const t = await opened();
  const host = container();
  const m = D.mountSeqViewer(host, t.ctx, { surface: stubSurface, ...spec });
  cleanups.unshift(() => m.dispose());
  await flush();
  return { ...t, host, m };
}

/** A made context's ports: the manual clock and the fixture measurer, the default spec, the canvas stub. */
function madeElement(tag = "seq-viewer") {
  const el = document.createElement(tag) as D.SeqViewerElement;
  const scheduler = manualScheduler();
  el.spec = headlessSpec();
  el.ports = { scheduler, text: fixtureMeasurer() };
  el.ui = { surface: stubSurface };
  return Object.assign(el, { advance: (ms: number) => scheduler.advance(ms) });
}
/** A task: what a deferred disposal waits for. */
const task = () => new Promise((r) => setTimeout(r, 5));

describe("the <seq-viewer> element (UX2; UI1's rows)", () => {
  it("makes a context on connect, disposes it a task after disconnect, makes a new one on reconnect", async () => {
    const el = madeElement();
    el.document = uiDoc();
    expect(el.context).toBeNull();
    const host = container();
    const told: unknown[] = [];
    host.addEventListener("seq-context", (e) => told.push((e as CustomEvent).detail));
    host.append(el);
    const first = el.context as SeqViewerContext;
    await flush();
    const shadow = el.shadowRoot as ShadowRoot;
    expect(first.document.shown.value?.metadata?.title).toBe("UI fixture");
    expect(shadow.querySelectorAll(".seqv-tab").length).toBe(2);
    expect(shadow.adoptedStyleSheets.length).toBe(1);
    // The `--seq-*` defaults on the document, so a host themes it from outside.
    expect(getComputedStyle(el).getPropertyValue("--seq-accent").trim()).not.toBe("");
    // `seq-context` told it once, frozen.
    expect(told).toEqual([{ context: first }]);
    expect(Object.isFrozen(told[0])).toBe(true);
    el.remove();
    // The UI goes at once; the context stays for a task.
    expect([shadow.childNodes.length, first.disposed]).toEqual([0, false]);
    await task();
    expect([first.disposed, el.context]).toEqual([true, null]);
    host.append(el);
    await flush();
    const second = el.context as SeqViewerContext;
    expect(second).not.toBe(first);
    expect(second.disposed).toBe(false);
    expect(shadow.querySelectorAll(".seqv-tab").length).toBe(2);
    expect(told.length).toBe(2);
    el.remove();
    await task();
    expect(second.disposed).toBe(true);
  });

  it("a story switch: an element removed with its page, then a new one (its chrome set after connect) renders", async () => {
    const first = madeElement();
    first.document = uiDoc();
    const page = container();
    page.append(first);
    await flush();
    page.remove();
    // The next story: a new host and a new element, its options set as a story's effects do, after it connected.
    const next = madeElement();
    next.document = uiDoc();
    const host = container();
    host.append(next);
    next.chrome = "clean";
    await flush();
    const shadow = next.shadowRoot as ShadowRoot;
    expect(shadow.querySelectorAll(".seqv-section").length).toBeGreaterThan(0);
    expect(next.mountedUI?.placed.has("Legend")).toBe(false);
    // A part the host mounts beside it shares its UI.
    const tabs = D.ViewTabs(host, next.context as SeqViewerContext, { ui: next.mountedUI ?? undefined });
    expect(host.querySelectorAll(".seqv-tab").length).toBe(2);
    tabs.dispose();
    next.remove();
    await task();
  });

  it("a spec from another copy of the viewer in the page (an HMR update) mounts: the viewer's entries are its own", async () => {
    vi.resetModules();
    const other = await import("../../src/index.ts");
    const spec = other.DefaultSeqViewerSpec();
    const el = madeElement();
    el.spec = { ...spec, config: [...(spec.config ?? []), [SeqViewerConfig.UI.Chrome, "full"]] };
    el.document = uiDoc();
    const errors: unknown[] = [];
    const host = container();
    host.addEventListener("seq-error", (e) => errors.push((e as CustomEvent).detail.error));
    host.append(el);
    await flush();
    expect(errors).toEqual([]);
    expect((el.shadowRoot as ShadowRoot).querySelectorAll(".seqv-section").length).toBeGreaterThan(0);
    el.remove();
    await task();
  });

  it("a mount that fails shows its error in the viewer's place and fires `seq-error`; a new spec mounts", async () => {
    const el = madeElement();
    const bad = headlessSpec();
    el.spec = {
      ...bad,
      registry: [...(bad.registry ?? []), { svs: { issueCodes: [{ code: "seqv.mine", severity: "error" }] } }],
    };
    el.document = uiDoc();
    const errors: unknown[] = [];
    const host = container();
    host.addEventListener("seq-error", (e) => errors.push((e as CustomEvent).detail.error));
    host.append(el);
    const shadow = el.shadowRoot as ShadowRoot;
    expect(errors.length).toBe(1);
    expect(shadow.querySelector("[role=alert].seqv-error")?.textContent).toContain('"seqv.mine" is reserved');
    expect([el.context, el.mountedUI]).toEqual([null, null]);
    el.spec = headlessSpec();
    await flush();
    expect(shadow.querySelector("[role=alert]")).toBeNull();
    expect(shadow.querySelectorAll(".seqv-section").length).toBeGreaterThan(0);
    el.remove();
    await task();
  });

  it("`mountedUI` is replaced on a remount, announced by `seq-ui`; its surfaces map stays the element's", async () => {
    const el = madeElement();
    el.document = uiDoc();
    const uis: unknown[] = [];
    const host = container();
    host.addEventListener("seq-ui", (e) => uis.push((e as CustomEvent).detail.ui));
    host.append(el);
    await flush();
    const first = el.mountedUI;
    expect(uis).toEqual([first]);
    el.chrome = "full";
    await flush();
    const second = el.mountedUI;
    expect(second).not.toBe(first);
    expect(uis).toEqual([first, second]);
    // A host's part on the first UI still reads the surfaces the second mount filled.
    expect(first?.surfaces).toBe(second?.surfaces);
    expect(second?.surfaces.size).toBeGreaterThan(0);
    el.remove();
    await task();
  });

  it("a moved element keeps the context it made, with its state, and does not load again", async () => {
    const el = madeElement();
    el.document = uiDoc();
    const a = container();
    const b = container();
    a.append(el);
    await flush();
    const ctx = el.context as SeqViewerContext;
    const load = vi.spyOn(ctx, "setDocument");
    ctx.selection.set(ctx.refs.toLoci("a", { start: 3, stop: 9 }), "api");
    ctx.views.setCurrent("msa");
    await ctx.actions.run("track.hide", "t.rows");
    await flush();
    // appendChild disconnects, then connects, in one task.
    b.append(el);
    await task();
    expect([el.context, ctx.disposed]).toEqual([ctx, false]);
    expect(el.shadowRoot?.querySelectorAll(".seqv-section").length).toBe(1);
    expect([ctx.views.current.value, ctx.selection.current.value.seq !== null]).toEqual(["msa", true]);
    expect(el.shadowRoot?.querySelector(".seqv-chips")?.textContent).toBe("Show Rows");
    expect(load).not.toHaveBeenCalled();
    el.remove();
    await task();
    expect(ctx.disposed).toBe(true);
  });

  it("`context`: a host's outlives the element; a new one remounts the UI on it and releases the made one", async () => {
    const { ctx } = await opened();
    const el = document.createElement("seq-viewer");
    el.context = ctx;
    el.ui = { surface: stubSurface };
    const host = container();
    host.append(el);
    await flush();
    const shadow = el.shadowRoot as ShadowRoot;
    expect(shadow.querySelector(".seqv-legend")).not.toBeNull();
    el.remove();
    await task();
    expect(ctx.disposed).toBe(false);
    host.append(el);
    await flush();
    expect(el.context).toBe(ctx);
    // A spec set later does not replace a host's context.
    el.spec = headlessSpec();
    el.ports = {};
    expect(el.context).toBe(ctx);
    // Light-DOM slots for host content.
    expect([...shadow.querySelectorAll("slot")].map((s) => s.name)).toEqual(["header", "toolbar"]);
    // Another host context: the UI is on it, and `document` loads into it.
    const other = await openedContext();
    cleanups.push(() => other.ctx.dispose());
    const doc = { ...uiDoc(), metadata: { title: "Other" } };
    el.document = doc;
    await flush();
    expect(ctx.document.shown.value?.metadata?.title).toBe("Other");
    el.context = other.ctx;
    await flush();
    expect(other.ctx.document.shown.value?.metadata?.title).toBe("Other");
    expect(shadow.querySelectorAll(".seqv-section").length).toBe(1);
    expect(subscriptions(ctx)).toBeLessThan(subscriptions(other.ctx));
    // Back on the first: the document it already has is not applied again (its edits since would be lost).
    const again = vi.spyOn(ctx, "setDocument");
    el.context = ctx;
    await flush();
    expect([el.context, again.mock.calls.length]).toEqual([ctx, 0]);
  });

  it("`spec` and `ports` remake the context it made, remount the UI and load again; the old one is disposed", async () => {
    const el = madeElement();
    el.document = uiDoc();
    container().append(el);
    await flush();
    const first = el.context as SeqViewerContext;
    first.ui.setStyle("chains.main", "compact");
    const text = { measure: () => 7 };
    el.ports = { scheduler: manualScheduler(), text };
    const second = el.context as SeqViewerContext;
    await flush();
    expect(second).not.toBe(first);
    expect(second.ui.style.value.size).toBe(0);
    expect(second.document.shown.value?.metadata?.title).toBe("UI fixture");
    expect(el.shadowRoot?.querySelectorAll(".seqv-tab").length).toBe(2);
    await task();
    expect([first.disposed, second.disposed]).toEqual([true, false]);
    el.spec = headlessSpec();
    await flush();
    expect(el.context).not.toBe(second);
    await task();
    expect(second.disposed).toBe(true);
    el.remove();
  });

  it("`ui`, `slots` and `chrome` remount the UI alone: the document, selection, hidden tracks and view are kept", async () => {
    const { ctx } = await opened();
    const el = document.createElement("seq-viewer");
    el.context = ctx;
    el.ui = { surface: stubSurface };
    container().append(el);
    await flush();
    const shadow = el.shadowRoot as ShadowRoot;
    ctx.selection.set(ctx.refs.toLoci("a", { start: 3, stop: 9 }), "api");
    ctx.actions.run("track.hide", "t.dom");
    await ctx.snapshots.apply("s2");
    const kept = () => [
      ctx.disposed,
      ctx.snapshots.current.value,
      ctx.selection.current.value.seq !== null,
      ctx.document.shown.value?.metadata?.title,
      shadow.querySelector(".seqv-chips")?.textContent,
    ];
    const before = kept();
    el.slots = { Legend: "none" };
    expect(shadow.querySelector(".seqv-legend")).toBeNull();
    expect(shadow.querySelector(".seqv-player")).not.toBeNull();
    el.ui = { surface: stubSurface, className: "host-theme" };
    expect(shadow.querySelector(".seqv-root")?.classList.contains("host-theme")).toBe(true);
    el.chrome = "clean";
    expect(shadow.querySelector(".seqv-tabs")).toBeNull();
    el.chrome = "full";
    await flush();
    expect(kept()).toEqual(before);
    expect(before).toEqual([false, "s2", true, "UI fixture", "Show Domains"]);
  });

  it("`ui`, `slots` and `chrome` keep a context the element made, told once, with its state", async () => {
    const el = madeElement();
    el.document = uiDoc();
    const host = container();
    const told: unknown[] = [];
    host.addEventListener("seq-context", (e) => told.push((e as CustomEvent).detail.context));
    host.append(el);
    await flush();
    const ctx = el.context as SeqViewerContext;
    await ctx.actions.run("track.hide", "t.dom");
    el.slots = { Legend: "none" };
    el.ui = { surface: stubSurface, className: "host-theme" };
    el.chrome = "clean";
    el.chrome = "full";
    await task();
    expect([el.context, ctx.disposed, told]).toEqual([ctx, false, [ctx]]);
    expect(el.shadowRoot?.querySelector(".seqv-chips")?.textContent).toBe("Show Domains");
    el.remove();
  });

  it("`document` and `data` load once per (context, object); a new one of either does not apply the other again", async () => {
    const el = madeElement();
    const doc = uiDoc();
    const data = { sequences: { a: "MKTAYIAKQRQISF" } } as unknown as D.SeqViewerElement["data"];
    el.document = doc;
    container().append(el);
    await flush();
    const ctx = el.context as SeqViewerContext;
    const setDocument = vi.spyOn(ctx, "setDocument");
    const setData = vi.spyOn(ctx, "setData").mockResolvedValue({ ok: true, issues: [] } as never);
    el.document = doc;
    el.document = doc;
    expect(setDocument).not.toHaveBeenCalled();
    el.document = { ...doc };
    expect(setDocument).toHaveBeenCalledTimes(1);
    el.data = data;
    el.data = data;
    expect(setData).toHaveBeenCalledTimes(1);
    // A new data object goes in alone: the document is not applied again, so the reader's edits stay.
    el.data = { ...(data as object) } as typeof data;
    expect([setDocument.mock.calls.length, setData.mock.calls.length]).toEqual([1, 2]);
    // A remount of the UI loads nothing again.
    el.slots = { Legend: "none" };
    el.chrome = "clean";
    expect([setDocument.mock.calls.length, setData.mock.calls.length]).toEqual([1, 2]);
    el.remove();
  });

  it("`legend` as an attribute or property is this UI's option: hidden leaves the part out, expanded opens it", async () => {
    const { ctx } = await opened();
    const el = document.createElement("seq-viewer");
    el.setAttribute("legend", "hidden");
    el.context = ctx;
    el.ui = { surface: stubSurface };
    container().append(el);
    await flush();
    const shadow = el.shadowRoot as ShadowRoot;
    expect([el.legend, shadow.querySelector(".seqv-legend")]).toEqual(["hidden", null]);
    el.legend = "expanded";
    await flush();
    expect((shadow.querySelector(".seqv-legend") as HTMLDetailsElement).open).toBe(true);
    // Not a mode: the UI spec's, else the context's (`UI.Legend`, here the fixtures' "expanded").
    el.setAttribute("legend", "sideways");
    expect(el.legend).toBe(ctx.config.get(SeqViewerConfig.UI.Legend));
    el.remove();
  });

  it("`chrome` as an attribute is this UI's option: the parts left out, the canvas parts kept, the context's config untouched", async () => {
    const { ctx } = await opened();
    ctx.report("seqv.test", "a reported problem");
    const el = document.createElement("seq-viewer");
    el.setAttribute("chrome", "clean");
    el.context = ctx;
    el.ui = { surface: stubSurface };
    container().append(el);
    await flush();
    const shadow = el.shadowRoot as ShadowRoot;
    expect([el.chrome, ctx.config.get(SeqViewerConfig.UI.Chrome)]).toEqual(["clean", "full"]);
    for (const c of [".seqv-tabs", ".seqv-section-header", ".seqv-legend", ".seqv-badge"])
      expect(shadow.querySelector(c)).toBeNull();
    expect([shadow.querySelector(".seqv-player"), shadow.querySelector(".seqv-surface")]).not.toContain(null);
    const h = at(ctx, "chains.main", "t.dom");
    pointer(ctx)("move", "chains.main", h.header, h.y);
    expect(shadow.querySelector(".seqv-track-actions")).toBeNull();
    // A slot's "default" brings one part back; the property returns the rest.
    el.slots = { Legend: "default" };
    expect([shadow.querySelector(".seqv-legend"), shadow.querySelector(".seqv-tabs")]).toEqual([
      expect.anything(),
      null,
    ]);
    el.chrome = "full";
    expect([shadow.querySelector(".seqv-tabs"), shadow.querySelector(".seqv-badge")]).not.toContain(null);
    // One clean UI does not turn another on the same context clean.
    const other = document.createElement("seq-viewer");
    other.context = ctx;
    other.ui = { surface: stubSurface };
    container().append(other);
    await flush();
    el.chrome = "clean";
    expect(other.shadowRoot?.querySelector(".seqv-tabs")).not.toBeNull();
    expect(ctx.config.get(SeqViewerConfig.UI.Chrome)).toBe("full");
  });

  it("takes properties set before it was defined", async () => {
    const el = document.createElement("seq-viewer-late") as D.SeqViewerElement;
    const made = madeElement();
    Object.assign(el, { spec: made.spec, ports: made.ports, ui: made.ui, document: uiDoc() });
    container().append(el);
    D.defineSeqViewer("seq-viewer-late");
    await flush();
    expect(el.context?.document.shown.value?.metadata?.title).toBe("UI fixture");
    el.remove();
  });

  it("re-dispatches selection, hover and snapshot as composed, bubbling events with frozen details; hover once a frame", async () => {
    const { ctx, advance } = await opened();
    const el = document.createElement("seq-viewer");
    el.context = ctx;
    el.ui = { surface: stubSurface };
    const outer = container();
    outer.append(el);
    await flush();
    const got: [string, unknown, boolean][] = [];
    for (const type of ["seq-selection", "seq-hover", "seq-snapshot"])
      outer.addEventListener(type, (e) => got.push([type, (e as CustomEvent).detail, e.composed]));
    // Not the states' replay on connect: a frame with no change tells nothing.
    advance(20);
    expect(got).toEqual([]);
    ctx.selection.set(ctx.refs.toLoci("a", { start: 3, stop: 9 }), "api");
    const res = at(ctx, "chains.main", "t.res");
    const move = pointer(ctx);
    // Three elements hovered between two frames: one event, for the last.
    for (const i of [10, 11, 12]) move("move", "chains.main", xOf(ctx, "chains.main", i), res.y);
    expect(got.filter((g) => g[0] === "seq-hover")).toEqual([]);
    advance(20);
    await ctx.snapshots.apply("s2");
    await flush();
    const hovers = got.filter((g) => g[0] === "seq-hover");
    expect(hovers.length).toBe(1);
    expect(got.every((g) => g[2])).toBe(true);
    const hover = hovers[0][1] as D.SeqHoverDetail;
    expect([hover.owner, hover.section, hover.track, hover.index]).toEqual(["pointer", "chains.main", "t.res", 12]);
    // The loci are the lit element's, closed; as coordinates in the shown set's systems, each with its descriptor.
    expect(hover.loci).not.toBeNull();
    expect(hover.coordinates?.systems).toEqual(["seqstar:a/unit"]);
    expect(Array.from(hover.coordinates?.runs ?? [])).toEqual([0, 12, 13]);
    expect(hover.systems["seqstar:a/unit"]).toMatchObject({ key: "seqstar:a/unit", kind: "unit", owner: "seqstar" });
    expect([Object.isFrozen(hover), Object.isFrozen(hover.coordinates), Object.isFrozen(hover.systems)]).toEqual([
      true,
      true,
      true,
    ]);
    // A hover that did not change the element tells nothing; a later one does not change the kept detail.
    move("move", "chains.main", xOf(ctx, "chains.main", 12) + 1, res.y);
    advance(20);
    expect(got.filter((g) => g[0] === "seq-hover").length).toBe(1);
    move("move", "chains.main", xOf(ctx, "chains.main", 30), res.y);
    advance(20);
    expect(got.filter((g) => g[0] === "seq-hover").length).toBe(2);
    expect([hover.index, Array.from(hover.coordinates?.runs ?? [])]).toEqual([12, [0, 12, 13]]);
    const sel = got.find((g) => g[0] === "seq-selection")?.[1] as D.SeqSelectionDetail;
    expect([sel.origin, sel.coordinates.systems, Array.from(sel.coordinates.runs), Object.keys(sel.systems)]).toEqual([
      "api",
      ["seqstar:a/unit"],
      [0, 3, 9],
      ["seqstar:a/unit"],
    ]);
    expect(Object.isFrozen(sel)).toBe(true);
    expect(got.find((g) => g[0] === "seq-snapshot")?.[1]).toEqual({ ref: "s2" });
    // Disposed with the UI: a removed element tells nothing more, a frame pending at its removal included.
    let own = 0;
    el.addEventListener("seq-hover", () => own++);
    move("move", "chains.main", xOf(ctx, "chains.main", 35), res.y);
    el.remove();
    move("move", "chains.main", xOf(ctx, "chains.main", 40), res.y);
    advance(20);
    expect([got.filter((g) => g[0] === "seq-hover").length, own]).toEqual([2, 0]);
  });

  it("`seq-hover`'s loci, coordinates and systems are a reader's cost: computed on the first read, kept, null once stale", async () => {
    const { ctx, advance } = await opened();
    const el = document.createElement("seq-viewer");
    el.context = ctx;
    el.ui = { surface: stubSurface };
    container().append(el);
    await flush();
    const details: D.SeqHoverDetail[] = [];
    el.addEventListener("seq-hover", (e) => details.push((e as CustomEvent).detail));
    const toLoci = vi.spyOn(ctx.hover, "toLoci");
    const coordinatesOf = vi.spyOn(ctx.systems, "coordinatesOf");
    const res = at(ctx, "chains.main", "t.res");
    const move = pointer(ctx);
    for (const i of [10, 20, 30]) {
      move("move", "chains.main", xOf(ctx, "chains.main", i), res.y);
      advance(20);
    }
    // Three events, nobody read their loci: none was made.
    expect([details.length, toLoci.mock.calls.length, coordinatesOf.mock.calls.length]).toEqual([3, 0, 0]);
    const last = details[2];
    expect(Array.from(last.coordinates?.runs ?? [])).toEqual([0, 30, 31]);
    expect([last.loci, last.systems["seqstar:a/unit"]]).not.toContain(undefined);
    expect([toLoci.mock.calls.length, coordinatesOf.mock.calls.length]).toEqual([1, 1]);
    // Read before the hover moved on, kept; first read after it, null (never another element's).
    move("move", "chains.main", xOf(ctx, "chains.main", 40), res.y);
    expect(Array.from(last.coordinates?.runs ?? [])).toEqual([0, 30, 31]);
    expect([details[1].index, details[1].loci, details[1].coordinates, details[1].systems]).toEqual([
      20,
      null,
      null,
      {},
    ]);
    el.remove();
  });

  it("a menu in the shadow root stays open on a press inside it and closes on one outside (`composedPath`)", async () => {
    const el = madeElement();
    el.document = uiDoc();
    container().append(el);
    await flush();
    const shadow = el.shadowRoot as ShadowRoot;
    (shadow.querySelector(".seqv-export button") as HTMLButtonElement).click();
    const item = shadow.querySelector("[role=menuitem]") as HTMLButtonElement;
    // Seen from the document, the press inside is on the element: only the composed path says it is in the menu.
    item.dispatchEvent(new PointerEvent("pointerdown", { bubbles: true, composed: true }));
    expect(shadow.querySelector("[role=menu]")).not.toBeNull();
    document.body.dispatchEvent(new PointerEvent("pointerdown", { bubbles: true, composed: true }));
    expect(shadow.querySelector("[role=menu]")).toBeNull();
    el.remove();
  });

  it("a made context registers the browser's `view.exportSvg`, so the Export part offers it; a spec's own is kept", async () => {
    const el = madeElement();
    el.document = uiDoc();
    container().append(el);
    await flush();
    const shadow = el.shadowRoot as ShadowRoot;
    const ctx = el.context as SeqViewerContext;
    expect(ctx.registry.get("actions", "view.exportSvg")).toBeDefined();
    expect(shadow.querySelector(".seqv-export button")?.textContent).toContain("Export");
    el.remove();
    // The host's own registration, here recording what it is offered.
    const offered: string[] = [];
    const own = madeElement();
    own.spec = {
      ...headlessSpec(),
      registry: [
        ...(headlessSpec().registry ?? []),
        { actions: [D.browserExportAction({ offer: (_svg, name) => void offered.push(name) })] },
      ],
    };
    own.document = uiDoc();
    container().append(own);
    await flush();
    const ownShadow = own.shadowRoot as ShadowRoot;
    (ownShadow.querySelector(".seqv-export button") as HTMLButtonElement).click();
    (ownShadow.querySelector("[role=menuitem]") as HTMLButtonElement).click();
    for (let i = 0; i < 40 && offered.length === 0; i++) await new Promise((r) => setTimeout(r, 50));
    expect(offered).toEqual(["chains.svg"]);
    own.remove();
  });
});

function RowToggle(host: HTMLElement, _ctx: SeqViewerContext, p: D.RowProps) {
  const b = D.button({ "aria-label": `Row ${p.displayRow}` }, () => {}, "3D");
  host.append(b);
  return { dispose: () => b.remove() };
}

/** The parts the context's view selects; each class marks the part that owns the nodes under it. */
const PARTS = [
  "seqv-tabs",
  "seqv-busy",
  "seqv-issues",
  "seqv-player",
  "seqv-description",
  "seqv-tooltip",
  "seqv-legend",
  "seqv-section-header",
  "seqv-chips",
  "seqv-track-actions",
  "seqv-row-controls",
  "seqv-state",
  "seqv-surface",
];
const partOf = (n: Node | null): string | null => {
  for (let e = n instanceof Element ? n : (n?.parentElement ?? null); e; e = e.parentElement)
    for (const c of PARTS) if (e.classList.contains(c)) return c;
  return null;
};
/** The parts whose nodes `fn` changed (records of the whole tree). */
async function mutated(root: Node, fn: () => unknown): Promise<string[]> {
  const records: MutationRecord[] = [];
  const mo = new MutationObserver((r) => records.push(...r));
  mo.observe(root, { subtree: true, childList: true, attributes: true, characterData: true });
  await fn();
  await flush();
  const out = new Set<string>();
  for (const r of [...records, ...mo.takeRecords()]) {
    // A part's old nodes taken apart after they left the document (a menu out of its hint) are nobody's to see.
    if (!r.target.isConnected) continue;
    const nodes = r.type === "childList" ? [...r.addedNodes, ...r.removedNodes] : [r.target];
    for (const n of nodes) if (n.nodeType !== Node.COMMENT_NODE) out.add(partOf(n) ?? partOf(r.target) ?? "?");
  }
  mo.disconnect();
  return [...out].sort();
}

describe("the DOM parts (UX2)", () => {
  it("a state change updates only the nodes its part owns", async () => {
    const t = await mounted({ slots: { RowControl: RowToggle } });
    const { ctx, host } = t;
    const move = pointer(ctx);
    const res = at(ctx, "chains.main", "t.res");
    // A pan, a selection, a highlight: nothing.
    expect(
      await mutated(host, () => {
        ctx.viewports.get("chains.main")?.panBy(40, 0, "api");
        ctx.selection.set(ctx.refs.toLoci("a", { start: 3, stop: 9 }), "api");
        ctx.highlight.set("host", ctx.refs.toLoci("dom"));
      }),
    ).toEqual([]);
    // Hovering positions: the tooltip only.
    expect(
      await mutated(host, () => {
        for (let i = 20; i < 30; i++) move("move", "chains.main", xOf(ctx, "chains.main", i), res.y);
      }),
    ).toEqual(["seqv-tooltip"]);
    const busy = ctx.tasks.busy as unknown as { set(v: boolean): void };
    expect(await mutated(host, () => busy.set(true))).toEqual(["seqv-busy"]);
    expect(await mutated(host, () => busy.set(false))).toEqual(["seqv-busy"]);
    expect(await mutated(host, () => ctx.report("seqv.test", "a reported problem"))).toEqual(["seqv-issues"]);
    expect(await mutated(host, () => ctx.snapshots.apply("s2"))).toEqual(["seqv-player"]);
    // Snapshot 1 hides a track and has its own description: the header (its row holds the chips), the description and
    // the legend follow, and the tooltip, which is placed again from the surface's offsets.
    expect(await mutated(host, () => ctx.snapshots.apply("s1"))).toEqual([
      "seqv-description",
      "seqv-legend",
      "seqv-player",
      "seqv-section-header",
      "seqv-tooltip",
    ]);
    // A scroll of the alignment's rows moves the row controls' content by a transform, nothing else.
    await ctx.actions.run("view.switch", "msa");
    await flush();
    expect(host.querySelectorAll(".seqv-row-control").length).toBeGreaterThan(0);
  });

  it("every part mounts and disposes with no subscription or node left", async () => {
    const { ctx } = await opened();
    // Warm up what the context makes lazily (per-section states of both views, snapshots' documents).
    await ctx.snapshots.apply("s1");
    await ctx.snapshots.apply("s2");
    ctx.views.setCurrent("msa");
    await flush();
    ctx.views.setCurrent("chains");
    ctx.report("seqv.test", "a reported problem");
    await flush();
    const before = subscriptions(ctx);
    const host = container();
    const ui = D.createDomUI(ctx, { surface: stubSurface });
    const single: [string, D.DomPart<never>, object][] = [
      ["ViewTabs", D.ViewTabs, {}],
      ["SectionHeader", D.SectionHeader, { section: "chains.main" }],
      ["SectionSwitches", D.SectionSwitches, { section: "chains.main" }],
      ["TrackActions", D.TrackActions, { section: "chains.main", track: "t.dom" }],
      ["HiddenChips", D.HiddenChips, { section: "chains.main" }],
      ["Tooltip", D.Tooltip, {}],
      ["Legend", D.Legend, {}],
      ["SnapshotPlayer", D.SnapshotPlayer, {}],
      ["Description", D.Description, {}],
      ["Issues", D.IssueBadge, {}],
      ["Busy", D.Busy, {}],
      ["EmptyState", D.EmptyState, {}],
      ["ErrorState", D.ErrorState, {}],
      ["ViewSurface", D.ViewSurface, {}],
    ];
    for (const [name, Part, props] of single) {
      const p = (Part as D.DomPart<object>)(host, ctx, { ui, ...props });
      expect([name, host.childNodes.length > 0]).toEqual([name, true]);
      expect([name, subscriptions(ctx) > before]).toEqual([name, true]);
      p.dispose();
      expect([name, host.childNodes.length, subscriptions(ctx)]).toEqual([name, 0, before]);
    }
    // The whole UI, with its leaf state in use: a track's menu open, the issue list open, a snapshot, the hover.
    const net: Record<string, number> = {};
    const count = (name: "addEventListener" | "removeEventListener", by: number) => {
      const real = document[name].bind(document) as (...a: unknown[]) => void;
      return vi.spyOn(document, name).mockImplementation(((type: string, ...rest: unknown[]) => {
        net[type] = (net[type] ?? 0) + by;
        real(type, ...rest);
      }) as typeof document.addEventListener);
    };
    const spies = [count("addEventListener", 1), count("removeEventListener", -1)];
    const m = D.mountSeqViewer(host, ctx, { surface: stubSurface, slots: { RowControl: RowToggle } });
    await flush();
    const h = at(ctx, "chains.main", "t.dom");
    pointer(ctx)("move", "chains.main", h.header, h.y);
    (host.querySelector(".seqv-track-actions button") as HTMLButtonElement).click();
    expect(host.querySelectorAll("[role=menuitem]").length).toBe(4);
    (host.querySelector(".seqv-badge") as HTMLButtonElement).click();
    expect(host.querySelector(".seqv-issue-list")).not.toBeNull();
    pointer(ctx)("move", "chains.main", xOf(ctx, "chains.main", 12), at(ctx, "chains.main", "t.res").y);
    expect(host.querySelector(".seqv-tooltip")).not.toBeNull();
    await ctx.snapshots.apply("s1");
    await flush();
    m.dispose();
    expect(host.childNodes.length).toBe(0);
    expect(subscriptions(ctx)).toBe(before);
    // No listener of an open popover or a surface is left on the document (a press outside would reach a disposed part).
    for (const spy of spies) spy.mockRestore();
    expect(Object.entries(net).filter(([, n]) => n !== 0)).toEqual([]);
  });

  it("the issue list closes when the issues go, and the badge comes back with it closed", async () => {
    const { ctx, host } = await mounted();
    ctx.report("seqv.test", "a reported problem");
    await flush();
    (host.querySelector(".seqv-badge") as HTMLButtonElement).click();
    expect(host.querySelector(".seqv-issue-list")).not.toBeNull();
    // A new document drops the reported ones.
    await ctx.setDocument({ ...uiDoc() });
    await flush();
    expect(host.querySelector(".seqv-badge")).toBeNull();
    ctx.report("seqv.test", "another");
    await flush();
    const badge = host.querySelector(".seqv-badge") as HTMLButtonElement;
    expect([badge.getAttribute("aria-expanded"), host.querySelector(".seqv-issue-list")]).toEqual(["false", null]);
  });

  it("the tooltip rewrites a header legend whose swatch changes under the same words, as React's does (VB-11)", async () => {
    const { ctx, host } = await mounted();
    const label = (color: string) => ({
      title: "Domains",
      lines: [],
      legend: { items: [{ label: "kinase", color, glyph: -1 }], ramps: [], more: 0 },
    });
    ctx.labels.current.set(label("#ff0000"));
    const swatch = () => (host.querySelector(".seqv-tooltip .seqv-swatch") as HTMLElement).style.background;
    const red = swatch();
    ctx.labels.current.set(label("#0000ff"));
    expect([red, swatch()]).toEqual(["rgb(255, 0, 0)", "rgb(0, 0, 255)"]);
  });

  it("notes alone make a quiet badge; the list shows errors and warnings first, the notes last under their heading", async () => {
    const { ctx, host } = await mounted();
    ctx.report("seqv.render.labels", "a quiet note");
    await flush();
    const badge = host.querySelector(".seqv-badge") as HTMLButtonElement;
    expect([badge.textContent, badge.dataset.severity]).toEqual(["1 note", "info"]);
    // Quiet by no fill, in the skin's quiet ink: text, never the hairlines' (LB-15, contrast 1.7:1 on white).
    const probe = document.createElement("span");
    probe.style.color = "var(--seq-ink-quiet)";
    badge.parentElement?.appendChild(probe);
    expect(getComputedStyle(badge).color).toBe(getComputedStyle(probe).color);
    probe.remove();
    // An unregistered code is an error.
    ctx.report("seqv.test", "a reported problem");
    await flush();
    expect([badge.textContent, badge.dataset.severity]).toEqual(["1 issue", "error"]);
    badge.click();
    const list = host.querySelector(".seqv-issue-list") as HTMLElement;
    const rows = [...list.querySelectorAll(".seqv-issue-heading, .seqv-issue")].map((e) =>
      e.classList.contains("seqv-issue") ? (e as HTMLElement).dataset.severity : e.textContent,
    );
    expect(rows).toEqual(["error", "Notes", "info"]);
    const note = list.querySelector('.seqv-issue[data-severity="info"]') as HTMLElement;
    const code = note.querySelector("code") as HTMLElement;
    expect(getComputedStyle(code).color).toBe(getComputedStyle(note).color);
  });

  it("slots replace parts, `none` leaves one out; no `title` attributes", async () => {
    const tabs: D.DomPart<D.PartProps> = (host) => {
      const el = D.h("div", { "data-testid": "host-tabs" }, "host tabs");
      host.append(el);
      return { dispose: () => el.remove() };
    };
    const hints: string[] = [];
    const hint: D.DomPart<D.HintProps> = (host, ctx, p) => {
      hints.push(p.label);
      return D.Hint(host, ctx, p);
    };
    const t = await mounted({
      slots: { ViewTabs: tabs, Legend: "none", Description: "none", Hint: hint, Menu: "none" },
    });
    const { host, ctx } = t;
    expect(host.querySelector("[data-testid=host-tabs]")).not.toBeNull();
    expect(host.querySelector(".seqv-tabs")).toBeNull();
    expect(host.querySelector(".seqv-legend")).toBeNull();
    expect(host.querySelector(".seqv-description")).toBeNull();
    expect(host.querySelector(".seqv-player")).not.toBeNull();
    expect(hints).toEqual(["Previous snapshot", "Play", "Next snapshot"]);
    // The menu at a header (Menu "none" keeps the default primitive), its items, a hidden track's chip.
    const h = at(ctx, "chains.main", "t.dom");
    pointer(ctx)("move", "chains.main", h.header, h.y);
    const trigger = host.querySelector(".seqv-track-actions button") as HTMLButtonElement;
    expect(trigger.getAttribute("aria-label")).toBe("Actions for Domains");
    trigger.click();
    const items = [...host.querySelectorAll<HTMLButtonElement>("[role=menuitem]")];
    expect(items.map((i) => i.textContent)).toEqual([
      "Show only what this track covers",
      "Focus",
      "Colors ▸",
      "Hide track",
    ]);
    // Focused, as a pointer press leaves it: the menu closes and still runs the item.
    items[3].focus();
    items[3].click();
    await flush();
    const chip = host.querySelector(".seqv-chips button") as HTMLButtonElement;
    expect(chip.textContent).toBe("Show Domains");
    chip.click();
    await flush();
    expect(host.querySelector(".seqv-chips")).toBeNull();
    expect(host.querySelectorAll("[title]").length).toBe(0);
  });

  it("`UI.TrackActions` orders the menu, with a host's own action; `none` (or nothing that applies) shows no menu", async () => {
    const { ctx } = await opened();
    const ran: string[] = [];
    ctx.config.set(SeqViewerConfig.UI.TrackActions, [
      "hide",
      "rows",
      { id: "log", label: "Log it", run: (_c, section, track) => void ran.push(`${section} ${track}`) },
      { id: "never", label: "Never", when: () => false, run: () => {} },
    ]);
    const host = container();
    let m = D.mountSeqViewer(host, ctx, { surface: stubSurface });
    const h = at(ctx, "chains.main", "t.dom");
    pointer(ctx)("move", "chains.main", h.header, h.y);
    (host.querySelector(".seqv-track-actions button") as HTMLButtonElement).click();
    const items = [...host.querySelectorAll<HTMLButtonElement>("[role=menuitem]")];
    expect(items.map((i) => i.textContent)).toEqual(["Hide track", "Log it"]);
    items[1].click();
    expect(ran).toEqual(["chains.main t.dom"]);
    m.dispose();
    for (const list of ["none", [{ id: "never", label: "Never", when: () => false, run: () => {} }]] as const) {
      ctx.config.set(SeqViewerConfig.UI.TrackActions, list);
      m = D.mountSeqViewer(host, ctx, { surface: stubSurface });
      pointer(ctx)("move", "chains.main", h.header, h.y + 1);
      pointer(ctx)("move", "chains.main", h.header, h.y);
      expect(host.querySelector(".seqv-track-actions button")).toBeNull();
      m.dispose();
    }
  });

  it("Shift+F10 on one item, then on another: the menu's item action runs on the cursor's item each time", async () => {
    const ran: ItemTarget[] = [];
    const show: ActionEntry<ItemTarget> = {
      name: "test:show",
      title: "Show it",
      when: () => true,
      run: (_c, t) => void ran.push(t),
    };
    const spec = headlessSpec();
    const ctx = createSeqViewer(
      { ...spec, registry: [...(spec.registry ?? []), { actions: [show] }] },
      { scheduler: manualScheduler(), text: fixtureMeasurer() },
    );
    cleanups.push(() => ctx.dispose());
    expect((await ctx.setDocument(uiDoc())).ok).toBe(true);
    ctx.config.set(SeqViewerConfig.UI.ItemActions, ["test:show"]);
    const host = container();
    const m = D.mountSeqViewer(host, ctx, { surface: stubSurface });
    cleanups.unshift(() => m.dispose());
    await flush();
    const S = "chains.main";
    const key = keys(ctx);
    ctx.input.focus(ctx.layout.sectionIndex(S), true);
    for (let i = 0; i < 8 && ctx.ui.keyCursor.value?.track !== "t.dom"; i++) key(S, "ArrowDown");
    const h = Hit.create();
    const items: number[] = [];
    for (const pos of [10, 60]) {
      // The keyboard's cursor on the item at `pos`; the track's menu by Shift+F10, its first entry chosen.
      ctx.layout.hitTest(S, xOf(ctx, S, pos), at(ctx, S, "t.dom").y, h);
      ctx.hover.setFromHit(h, "keyboard");
      items.push(ctx.hover.live.current.item.item);
      key(S, "F10", true);
      await flush();
      const entries = [...host.querySelectorAll<HTMLButtonElement>("[role^=menuitem]")];
      expect(entries[0]?.textContent).toBe("Show it");
      entries[0].click();
      await flush();
    }
    expect(items[0]).not.toBe(items[1]);
    expect(ran.map((t) => t.item)).toEqual(items);
  });

  it("the track's menu sits on the hovered line: a rows track's row, on every hover; it stays while held", async () => {
    const { host, ctx } = await mounted();
    const box = () => host.querySelector(".seqv-track-actions") as HTMLElement | null;
    const move = pointer(ctx);
    // An ordinary track: its band.
    const dom = lines(ctx, "chains.main", "t.dom")[0];
    move("move", "chains.main", dom.header, dom.y);
    expect([box()?.style.top, box()?.style.height]).toEqual([`${dom.top}px`, `${dom.height}px`]);
    ctx.input.leave(ctx.layout.sectionIndex("chains.main"));
    ctx.views.setCurrent("msa");
    await flush();
    const section = "msa.main";
    const rows = lines(ctx, section, "t.rows");
    expect(rows.length).toBeGreaterThan(4);
    // Row 3, then row 1, then row 4 of the same track: the same box, on each hovered row.
    move("move", section, rows[3].header, rows[3].y);
    const first = box();
    for (const k of [3, 1, 4]) {
      move("move", section, rows[k].header, rows[k].y);
      expect([k, box()?.style.top, box()?.style.height]).toEqual([k, `${rows[k].top}px`, `${rows[k].height}px`]);
    }
    expect(box()).toBe(first);
    // The pointer leaves the surface for the menu (one dispatch): it stays on row 4, also when the view moves under it.
    ctx.input.leave(ctx.layout.sectionIndex(section));
    box()?.dispatchEvent(new PointerEvent("pointerenter"));
    await flush();
    ctx.viewports.get(section)?.panBy(10, 0, "api");
    expect([box(), box()?.style.top]).toEqual([first, `${rows[4].top}px`]);
    box()?.dispatchEvent(new PointerEvent("pointerleave"));
    await flush();
    expect(box()).toBeNull();
  });

  it("the legend on the canvas's grid: its toggle at the labels' inset, its items from the data left", async () => {
    const { host, ctx } = await mounted();
    const surface = host.querySelector(".seqv-surface") as HTMLElement;
    const left = (sel: string) => (host.querySelector(sel) as HTMLElement).getBoundingClientRect().left;
    const dataLeft = () => ctx.layout.get("chains.main")?.current.dataLeft;
    expect(host.querySelectorAll(".seqv-legend-item").length).toBeGreaterThan(0);
    for (const headerWidth of [150, 220]) {
      ctx.setRenderProps({ headerWidth });
      ctx.layout.get("chains.main"); // As the canvas reads it on its next frame (the layout builds on read).
      await flush();
      const x0 = surface.getBoundingClientRect().left;
      expect([headerWidth, left(".seqv-legend-item") - x0]).toEqual([headerWidth, dataLeft()]);
      expect(left(".seqv-legend summary") - x0).toBe(10);
    }
  });
});
