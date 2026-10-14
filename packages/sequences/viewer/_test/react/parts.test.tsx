/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The React parts' content (UX3; seqspec:interaction#the-default-ui, seqspec:interaction#the-legend): the legend's "+n
// more" and ramps with their captions (UX1's helpers drawn), Export through the browser's action, Busy while regions
// are in flight, the row controls (their layout read outside render; their scroll a custom property written when it
// changes), and the whole UI leaving no subscription and no document listener behind.

import { createSeqViewSpec } from "@molstar/seq-core/builder";
import { act } from "react";
import { describe, expect, it, vi } from "vitest";
import "../../src/skin.css";
import { browserExportAction } from "../../src/dom/index.js";
import { BUSY_HOLD_MS } from "../../src/dom/shared.js";
import { createSeqViewer, KernelRenderers, type LayerRenderer, manualScheduler } from "../../src/index.js";
import { type RowProps, SeqViewerLayout, SeqViewerProvider, type SurfaceFactory } from "../../src/react/index.js";
import { fixtureMeasurer } from "../../src/testing/index.js";
import {
  at,
  harness,
  headlessSpec,
  lines,
  mount,
  openedContext,
  pointer,
  settle,
  stubSurface,
  subscriptions,
  uiDoc,
  WIDTH,
  xOf,
} from "./fixtures.js";

const { opened, cleanups, parent } = harness();
const A = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQMKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ";

/** A context showing `doc`, the default spec with some kernel renderers replaced by name, mounted. */
async function mounted(
  doc: ReturnType<ReturnType<typeof createSeqViewSpec>["getDocument"]>,
  replace: Record<string, LayerRenderer> = {},
) {
  const spec = headlessSpec();
  const [first, ...rest] = spec.registry ?? [];
  const registry = [
    { ...first, renderers: (KernelRenderers as LayerRenderer[]).map((r) => replace[r.name] ?? r) },
    ...rest,
  ];
  const ctx = createSeqViewer({ ...spec, registry }, { scheduler: manualScheduler(), text: fixtureMeasurer() });
  const done = ctx.setDocument(doc);
  for (const ref of ctx.viewports.sections.value) ctx.viewports.get(ref)?.setSize(WIDTH, 220, 1);
  await done;
  const m = await mount(ctx, { surface: stubSurface });
  cleanups.push(() => ctx.dispose(), m.unmount);
  return { ctx, container: m.container };
}

/** `n` ranges, a category each. */
function manyDoc(n: number) {
  const b = createSeqViewSpec({ title: "Many" });
  const a = b.data.sequence({ ref: "a", label: "A", value: A });
  const start = Array.from({ length: n }, (_, k) => k);
  const d = b.data.annotation({
    ref: "d",
    geometry: "ranges",
    target: a,
    value: { start, stop: start.map((x) => x + 1), label: start.map((k) => `cat${k}`) },
  });
  const v = b.view({ ref: "v", label: "V" });
  v.block({ sequence: a });
  v.track({ ref: "v.t", label: "Many" }).ranges({ data: d });
  return b.getDocument();
}
function heatDoc() {
  const b = createSeqViewSpec({ title: "Heat" });
  const a = b.data.sequence({ ref: "a", label: "A", value: A });
  const d = b.data.annotation({
    ref: "d",
    geometry: "values",
    target: a,
    value: { index: [0, 5, 9], value: [0, 5, 20] },
  });
  b.view({ ref: "v", label: "V" }).block({ sequence: a }).track({ label: "Heat" }).color({ data: d });
  return b.getDocument();
}

const inputs = (ramp: HTMLElement) =>
  [...ramp.querySelectorAll<HTMLInputElement>(".seqv-range-input")].map((i) => i.value);

describe("the legend", () => {
  it('a group shows `UI.LegendItems` items and "+n more"; the control expands and collapses it and keeps focus', async () => {
    const { container, ctx } = await mounted(manyDoc(30));
    const items = () => container.querySelectorAll(".seqv-legend-item").length;
    const more = () => container.querySelector(".seqv-legend-more") as HTMLButtonElement;
    expect([items(), more().textContent, more().getAttribute("aria-expanded")]).toEqual([24, "+6 more", "false"]);
    more().focus();
    await act(async () => more().click());
    expect([items(), more().textContent, more().getAttribute("aria-expanded"), [...ctx.ui.legendOpen.value]]).toEqual([
      30,
      "fewer",
      "true",
      ["v.main/v.t"],
    ]);
    expect(document.activeElement).toBe(more());
    await act(async () => more().click());
    expect([items(), more().textContent, document.activeElement === more()]).toEqual([24, "+6 more", true]);
    // The open state is the context's: another UI on it agrees, and the details element's own state is kept.
    await act(async () => ctx.ui.toggleLegend("v.main/v.t"));
    expect(items()).toBe(30);
    const legend = container.querySelector(".seqv-legend") as HTMLDetailsElement;
    legend.open = false;
    await settle();
    await act(async () => ctx.ui.toggleLegend("v.main/v.t"));
    expect((container.querySelector(".seqv-legend") as HTMLDetailsElement).open).toBe(false);
  });

  it("a group at the cap has no control", async () => {
    const { container } = await mounted(manyDoc(24));
    expect(container.querySelectorAll(".seqv-legend-item").length).toBe(24);
    expect(container.querySelector(".seqv-legend-more")).toBeNull();
  });

  it("a continuous palette is a gradient bar with its two ends, not a list; a caption stands beside it", async () => {
    const { container } = await mounted(heatDoc());
    const ramp = container.querySelector(".seqv-ramp") as HTMLElement;
    expect(inputs(ramp)).toEqual(["0", "20"]);
    expect(container.querySelectorAll(".seqv-legend-item").length).toBe(0);
    expect((ramp.querySelector(".seqv-ramp-bar") as HTMLElement).style.background).toContain("linear-gradient(");
    const captioned: LayerRenderer = {
      ...(KernelRenderers.find((r) => r.name === "color") as LayerRenderer),
      legend: (_p, _t, sink) =>
        sink.ramp(Float64Array.of(0, 0x112233, 1, 0xffeeaa), 0, 12.3456, "contacts (number of kinds)"),
    };
    const t = await mounted(heatDoc(), { color: captioned });
    const cap = t.container.querySelector(".seqv-ramp") as HTMLElement;
    expect([inputs(cap), cap.textContent]).toEqual([["0", "12.3"], "contacts (number of kinds)"]);
    expect([...cap.children].map((c) => c.className)).toEqual([
      "seqv-range-input",
      "seqv-ramp-track",
      "seqv-range-input",
      "seqv-quiet",
    ]);
  });

  it("a sample size is a dot of that size", async () => {
    const sized: LayerRenderer = {
      ...(KernelRenderers.find((r) => r.name === "ranges") as LayerRenderer),
      legend: (_p, _t, sink) => sink.size("a pair", 6, 0),
    };
    const { container } = await mounted(manyDoc(3), { ranges: sized });
    const dot = container.querySelector(".seqv-legend-item .seqv-size") as HTMLElement;
    expect([dot.style.width, dot.style.height, dot.style.borderRadius]).toEqual(["6px", "6px", "50%"]);
    expect(dot.parentElement?.textContent).toBe("a pair");
  });
});

describe("the color range editor", () => {
  it("sets the ends from the fields on blur, marks the data's span, refuses a bad number, resets", async () => {
    const { container, ctx } = await mounted(heatDoc());
    const field = (end: string) =>
      container.querySelector(`.seqv-range-input[data-seqv-key^="range-${end}"]`) as HTMLInputElement;
    const type = async (end: string, text: string) => {
      const input = field(end);
      await act(async () => input.focus());
      input.value = text;
      await act(async () => input.blur());
    };
    expect(container.querySelector('[data-seqv-key^="range-reset"]')).toBeNull();
    await type("hi", "40");
    expect([...ctx.ui.colorRanges.value.values()]).toEqual([[0, 40]]);
    const mark = container.querySelector(".seqv-ramp-data") as HTMLElement;
    expect([mark.style.left, mark.style.width]).toEqual(["0%", "50%"]);
    await type("hi", "-3");
    expect(field("hi").value).toBe("40");
    await act(async () => (container.querySelector('[data-seqv-key^="range-reset"]') as HTMLButtonElement).click());
    expect([ctx.ui.colorRanges.value.size, field("hi").value]).toEqual([0, "20"]);
  });
});

describe("Export and Busy", () => {
  it("the Export part is a menu in the bar when the context registers `view.exportSvg`; its item runs the action", async () => {
    const without = await opened();
    expect(without.container.querySelector(".seqv-export")).toBeNull();
    const offered: string[] = [];
    const spec = headlessSpec();
    const ctx = createSeqViewer(
      {
        ...spec,
        registry: [
          ...(spec.registry ?? []),
          { actions: [browserExportAction({ offer: (_svg, name) => void offered.push(name) })] },
        ],
      },
      { scheduler: manualScheduler(), text: fixtureMeasurer() },
    );
    await ctx.setDocument(uiDoc());
    const m = await mount(ctx, { surface: stubSurface });
    cleanups.push(() => ctx.dispose(), m.unmount);
    const trigger = m.container.querySelector(".seqv-export button") as HTMLButtonElement;
    expect([trigger.textContent, trigger.getAttribute("aria-haspopup")]).toEqual(["Export", "menu"]);
    await act(async () => trigger.click());
    const item = m.container.querySelector("[role=menuitem]") as HTMLButtonElement;
    expect(item.textContent).toBe("Export as SVG");
    // A host syncing its URL from `committed` sees nothing of the export (seqspec:viewer#the-render-backend).
    let commits = 0;
    const sub = ctx.viewports.get("chains.main")?.committed.subscribe(() => commits++);
    const seen = commits;
    await act(async () => item.click());
    for (let i = 0; i < 40 && offered.length === 0; i++) await act(async () => new Promise((r) => setTimeout(r, 50)));
    expect(offered).toEqual(["chains.svg"]);
    sub?.unsubscribe();
    expect(commits).toBe(seen);
    expect([m.container.querySelector("[role=menu]"), document.activeElement]).toEqual([null, trigger]);
  });

  for (const [label, ui] of [
    ["a bar with its parts", {}],
    ["a bar with nothing else in it (4 px)", { chrome: "clean", components: { Search: "none" } }],
  ] as const)
    it(`Busy takes no room: nothing moves when it shows or hides, in ${label}`, async () => {
      const { container, ctx } = await opened(ui as never);
      const busy = ctx.tasks.busy as unknown as { set(v: boolean): void };
      // Every box of the viewer but the spinner and the bar's spacer (which gives up what the spinner's gap needs).
      const probe = () =>
        Array.from(container.querySelectorAll(".seqv-root, .seqv-root *"))
          .filter((e) => !e.matches(".seqv-busy, .seqv-spacer"))
          .map((e) => {
            const r = e.getBoundingClientRect();
            return [r.left, r.top, r.width, r.height].map((x) => Math.round(x * 100) / 100);
          });
      const before = probe();
      const barHeight = (container.querySelector(".seqv-bar") as HTMLElement).getBoundingClientRect().height;
      if (label.includes("nothing else")) expect(barHeight).toBeLessThan(10);
      await act(async () => busy.set(true));
      const spinner = container.querySelector(".seqv-busy") as HTMLElement;
      expect(spinner).not.toBeNull();
      expect(probe()).toEqual(before);
      // The spinner is drawn at its full size.
      const r = spinner.getBoundingClientRect();
      expect([r.width, r.height]).toEqual([14, 14]);
      await act(async () => busy.set(false));
      expect(probe()).toEqual(before);
    });

  it("Busy shows while tasks run and while regions in view are in flight, held a moment after the last region", async () => {
    const { container, ctx } = await opened();
    const busy = () => container.querySelector(".seqv-busy");
    expect(busy()).toBeNull();
    await act(async () => {
      ctx.streams.stats.current.inFlight = 2;
      (ctx.streams.stats as unknown as { bump(): void }).bump();
    });
    expect(busy()).not.toBeNull();
    const stats = (n: number) =>
      act(async () => {
        ctx.streams.stats.current.inFlight = n;
        (ctx.streams.stats as unknown as { bump(): void }).bump();
      });
    // The last region landed: held a moment (a stream's next region keeps it on, no flicker), then gone.
    await stats(0);
    expect(busy()).not.toBeNull();
    await act(async () => new Promise((r) => setTimeout(r, BUSY_HOLD_MS / 2)));
    await stats(1);
    await stats(0);
    await act(async () => new Promise((r) => setTimeout(r, BUSY_HOLD_MS / 2 + 20)));
    expect(busy()).not.toBeNull();
    await act(async () => new Promise((r) => setTimeout(r, BUSY_HOLD_MS)));
    expect(busy()).toBeNull();
    const tasks = ctx.tasks.busy as unknown as { set(v: boolean): void };
    await act(async () => tasks.set(true));
    expect(busy()).not.toBeNull();
    await act(async () => tasks.set(false));
    expect(busy()).toBeNull();
  });
});

function RowToggle({ displayRow }: RowProps) {
  return (
    <button type="button" aria-label={`Row ${displayRow}`}>
      3D
    </button>
  );
}

describe("the row controls", () => {
  it("one beside each visible row, a view's own; a vertical scroll moves them by a custom property, written on change", async () => {
    // A surface shorter than the alignment's rows, so that it scrolls.
    const short: SurfaceFactory = (c, section, host) => {
      const s = stubSurface(c, section, host);
      c.viewports.get(section)?.setSize(WIDTH, 60, 1);
      return s;
    };
    const { container, ctx } = await opened({ surface: short, components: { RowControl: RowToggle } });
    await act(async () => ctx.views.setCurrent("msa"));
    await settle();
    expect(container.querySelectorAll(".seqv-row-control").length).toBeGreaterThan(0);
    expect(lines(ctx, "msa.main", "t.rows").length).toBeGreaterThan(0);
    const inner = container.querySelector(".seqv-row-controls > div") as HTMLElement;
    // Every write of the property, even of the value it has (a mutation observer would not see that).
    let owner = Object.getPrototypeOf(inner.style);
    while (owner && !Object.hasOwn(owner, "setProperty")) owner = Object.getPrototypeOf(owner);
    const writes = vi.spyOn(owner as CSSStyleDeclaration, "setProperty");
    const scrolls = () => writes.mock.calls.filter(([name]) => name === "--seqv-scroll-y").length;
    const vp = ctx.viewports.get("msa.main");
    // A vertical scroll writes the offset once, a horizontal pan nothing, a vertical scroll again its new offset.
    await act(async () => vp?.scrollTo(0, 12, "api"));
    await settle();
    expect([scrolls(), inner.style.getPropertyValue("--seqv-scroll-y")]).toEqual([1, "-12px"]);
    await act(async () => vp?.panBy(30, 0, "api"));
    // A frame of another kind (the version moved, the scroll did not): nothing is written.
    await act(async () => (vp?.live as unknown as { bump(): void } | undefined)?.bump());
    await settle();
    expect(scrolls()).toBe(1);
    await act(async () => vp?.scrollTo(0, 24, "api"));
    await settle();
    const total = scrolls();
    writes.mockRestore();
    expect([total, inner.style.getPropertyValue("--seqv-scroll-y"), vp?.live.current.scrollY]).toEqual([
      2,
      "-24px",
      24,
    ]);
  });

  it("reads the layout in a subscription, never in a render: a stale layout rebuilt on read notifies no component mid-render", async () => {
    const errors = vi.spyOn(console, "error").mockImplementation(() => {});
    try {
      const { ctx } = await openedContext();
      cleanups.push(() => ctx.dispose());
      await ctx.actions.run("view.switch", "msa");
      // A parent that renders the UI again with another slot component each time: the section frames and the row controls
      // render again with it.
      const Tabs = [() => <div>a</div>, () => <div>b</div>];
      const p = parent((i) => (
        <SeqViewerProvider
          context={ctx}
          ui={{ surface: stubSurface, components: { RowControl: RowToggle, ViewTabs: Tabs[i % 2] } }}
        >
          <SeqViewerLayout />
        </SeqViewerProvider>
      ));
      await p.render();
      expect(p.container.querySelectorAll(".seqv-row-control").length).toBeGreaterThan(0);
      errors.mockClear();
      // The layout is stale until something reads it.
      const reads: string[] = [];
      const get = ctx.layout.get.bind(ctx.layout);
      const props = ctx.layout.props.bind(ctx.layout);
      const rows = ctx.layout.visibleRows.bind(ctx.layout);
      const stack = () => new Error().stack?.split("\n").find((l) => /RowControls/.test(l)) ?? "";
      vi.spyOn(ctx.layout, "get").mockImplementation((s) => {
        reads.push(stack());
        return get(s);
      });
      vi.spyOn(ctx.layout, "props").mockImplementation((s) => {
        reads.push(stack());
        return props(s);
      });
      vi.spyOn(ctx.layout, "visibleRows").mockImplementation((s) => {
        reads.push(stack());
        return rows(s);
      });
      await act(async () => ctx.setRenderProps({ headerWidth: 200 }));
      await p.render();
      await p.render();
      expect(reads.filter((r) => r !== "")).toEqual([]);
      expect(errors.mock.calls.filter((c) => String(c[0]).includes("while rendering"))).toEqual([]);
    } finally {
      vi.restoreAllMocks();
    }
  });

  it("follows the rows in view and a header width, and a remount of the dataRow", async () => {
    const { container, ctx } = await opened({ components: { RowControl: RowToggle } });
    await act(async () => ctx.views.setCurrent("msa"));
    await settle();
    const box = () => container.querySelector(".seqv-row-controls") as HTMLElement;
    const before = box().style.left;
    await act(async () => {
      ctx.setRenderProps({ headerWidth: 260 });
      ctx.layout.get("msa.main");
    });
    await settle();
    // The box's left is the header's end less the control column: it follows the header.
    expect(Number.parseFloat(box().style.left)).toBeGreaterThanOrEqual(Number.parseFloat(before));
    const labels = [...container.querySelectorAll(".seqv-row-control button")].map((b) => b.getAttribute("aria-label"));
    expect(new Set(labels).size).toBe(labels.length);
  });
});

describe("the tooltip's place", () => {
  it("is right of and below the pointer and within the viewer, by two properties the skin makes a transform of", async () => {
    const { container, ctx } = await opened();
    const x = xOf(ctx, "chains.main", 10);
    const y = at(ctx, "chains.main", "t.res").y;
    await act(async () => pointer(ctx)("move", "chains.main", x, y));
    const tip = container.querySelector(".seqv-tooltip") as HTMLElement;
    const surface = container.querySelector(".seqv-surface") as HTMLElement;
    const r = tip.getBoundingClientRect();
    const s = surface.getBoundingClientRect();
    expect(r.left).toBeGreaterThanOrEqual(s.left + x);
    expect(r.top).toBeGreaterThanOrEqual(s.top + y);
    expect(r.right).toBeLessThanOrEqual(
      (container.querySelector(".seqv-root") as HTMLElement).getBoundingClientRect().right + 1,
    );
    expect(tip.style.getPropertyValue("--seqv-x")).toMatch(/^\d+px$/);
    expect(tip.style.transform).toBe("");
    // At the right edge it flips to the left of the pointer.
    await act(async () => pointer(ctx)("move", "chains.main", xOf(ctx, "chains.main", 85), y));
    const flipped = (container.querySelector(".seqv-tooltip") as HTMLElement).getBoundingClientRect();
    const px = s.left + xOf(ctx, "chains.main", 85);
    expect(flipped.right).toBeLessThanOrEqual(px + 1);
  });
});

describe("the whole UI leaves nothing behind", () => {
  it("every subscription and every document listener it added goes with its unmount, with its leaf state in use", async () => {
    const { ctx } = await openedContext();
    await ctx.snapshots.apply("s1");
    await ctx.snapshots.apply("s2");
    ctx.views.setCurrent("msa");
    await settle();
    ctx.views.setCurrent("chains");
    ctx.report("seqv.test", "a reported problem");
    await settle();
    const before = subscriptions(ctx);
    const net: Record<string, number> = {};
    const count = (name: "addEventListener" | "removeEventListener", by: number) => {
      const real = document[name].bind(document) as (...a: unknown[]) => void;
      return vi.spyOn(document, name).mockImplementation(((type: string, ...rest: unknown[]) => {
        net[type] = (net[type] ?? 0) + by;
        real(type, ...rest);
      }) as typeof document.addEventListener);
    };
    const spies = [count("addEventListener", 1), count("removeEventListener", -1)];
    const m = await mount(ctx, { surface: stubSurface, components: { RowControl: RowToggle } });
    const h = at(ctx, "chains.main", "t.dom");
    await act(async () => pointer(ctx)("move", "chains.main", h.header, h.y));
    await act(async () => (m.container.querySelector(".seqv-track-actions button") as HTMLButtonElement).click());
    expect(m.container.querySelectorAll("[role=menuitem]").length).toBe(4);
    await act(async () => (m.container.querySelector(".seqv-badge") as HTMLButtonElement).click());
    expect(m.container.querySelector(".seqv-issue-list")).not.toBeNull();
    await act(async () =>
      pointer(ctx)("move", "chains.main", xOf(ctx, "chains.main", 12), at(ctx, "chains.main", "t.res").y),
    );
    expect(m.container.querySelector(".seqv-tooltip")).not.toBeNull();
    await act(async () => ctx.snapshots.apply("s1"));
    await settle();
    await m.unmount();
    expect(m.container.childNodes.length).toBe(0);
    expect(subscriptions(ctx)).toBe(before);
    for (const spy of spies) spy.mockRestore();
    expect(Object.entries(net).filter(([, n]) => n !== 0)).toEqual([]);
    ctx.dispose();
  });
});
