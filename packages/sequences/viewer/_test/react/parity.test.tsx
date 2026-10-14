/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// React/DOM parity (UX3; seqspec:interaction#the-dom-element): for every case of the format's examples (and the UI
// fixture), the React and
// DOM UIs mounted on the same context produce the same part structure and the same text: in each view, at each
// snapshot, with a track's menu open, with the hover's label, with the issue list open. Compared: elements, their
// classes, roles, ARIA and data attributes (ids, `aria-controls` and styles aside), and their text.

import { createSeqViewSpec } from "@molstar/seq-core/builder";
import { act, type ReactNode } from "react";
import { createRoot } from "react-dom/client";
import { afterAll, afterEach, describe, expect, it } from "vitest";
import { page, userEvent } from "vitest/browser";
import { CASES } from "../../../core/_test/format/cases.ts";
import { fixtureLoaders, MolstarFixture } from "../../../core/_test/format/load-fixtures.ts";
import * as D from "../../src/dom/index.ts";
import { mountSeqViewer } from "../../src/dom/index.ts";
import {
  createSeqViewer,
  DefaultSeqViewerSpec,
  EntryKind,
  manualScheduler,
  SeqViewerConfig,
  type SeqViewerContext,
} from "../../src/index.ts";
// The skin, for the legend's grid (its layout is CSS).
import "../../src/skin.css";

import { AxisGeometry } from "../../src/axis.ts";
import * as R from "../../src/react/index.ts";
import { fixtureMeasurer } from "../../src/testing/index.ts";
import { keys } from "../dom/fixtures.ts";
import { distinctDoc, kindsDoc } from "../legend-fixtures.ts";
import { headlessSpec, lines, mount, pointer, settle, stubSurface, subscriptions, uiDoc, WIDTH } from "./fixtures.tsx";

const ATTRS = /^(class|role|type|disabled|tabindex|href|target|rel|open|aria-(?!controls$).+|data-(?!testid$).+)$/;
/** An element tree as lines: tag, the compared attributes (sorted), then adjacent text merged into one line. */
function shape(n: Element, depth = 0, out: string[] = []): string[] {
  const pad = "  ".repeat(depth);
  const attrs = [...n.attributes].filter((a) => ATTRS.test(a.name)).map((a) => `${a.name}="${a.value}"`);
  out.push(`${pad}<${n.tagName.toLowerCase()}${attrs.length ? ` ${attrs.sort().join(" ")}` : ""}>`);
  let text = "";
  const flush = () => {
    if (text) out.push(`${pad}  ${JSON.stringify(text)}`);
    text = "";
  };
  for (const c of n.childNodes) {
    if (c.nodeType === Node.TEXT_NODE) text += c.textContent ?? "";
    else if (c instanceof Element) {
      flush();
      shape(c, depth + 1, out);
    }
  }
  flush();
  return out;
}

/** The default spec with `full` chrome (the legend's tests; the default is clean). */
function fullChrome() {
  const spec = DefaultSeqViewerSpec();
  return { ...spec, config: [...(spec.config ?? []), [SeqViewerConfig.UI.Chrome, "full"] as const] };
}

async function caseContext(doc: Parameters<SeqViewerContext["setDocument"]>[0]) {
  const spec = headlessSpec();
  const ctx = createSeqViewer(
    { ...spec, registry: [...(spec.registry ?? []), { svs: MolstarFixture }] },
    { scheduler: manualScheduler(), loaders: fixtureLoaders(), text: fixtureMeasurer() },
  );
  await ctx.setDocument(doc);
  return ctx;
}

const LONG = "SRC_HUMAN/151-248 with a name longer than any label column";
const LONG_TRACK = "Residues, as the sequence gives them, letter by letter";
/** An alignment's rows with one name too long for the label column, and a track's label too long for it. */
function cutDoc() {
  const b = createSeqViewSpec({ title: "Cut" });
  const aln = b.data.alignment({ ref: "aln", value: ["MKTAYIAKQR", "MKTAHIAKQR"], labels: [LONG, "B"] });
  const v = b.view({ ref: "v", label: "V" }).block({ alignment: aln });
  v.track({ ref: "v.letters", label: LONG_TRACK }).letters();
  v.track({ ref: "v.rows", label: "Rows" }).rows();
  return b.getDocument();
}

const cleanups: (() => Promise<void> | void)[] = [];
afterEach(async () => {
  for (const c of cleanups.splice(0)) await c();
});

/** Both UIs on one context; `same(step)` compares them and returns the shape (for checks that it is not empty). */
async function both(ctx: SeqViewerContext) {
  const react = await mount(ctx, { surface: stubSurface });
  const dom = document.createElement("div");
  dom.style.width = `${WIDTH}px`;
  document.body.appendChild(dom);
  const m = mountSeqViewer(dom, ctx, { surface: stubSurface });
  cleanups.push(
    react.unmount,
    () => {
      m.dispose();
      dom.remove();
    },
    () => ctx.dispose(),
  );
  const roots = () => [react.container, dom].map((c) => c.firstElementChild as Element);
  return {
    roots,
    same: async (step: string) => {
      await settle();
      const [r, d] = roots().map((e) => shape(e).join("\n"));
      expect(`${step}\n${d}`).toBe(`${step}\n${r}`);
      return r;
    },
    /** A click on the element at `selector` in both UIs. */
    click: async (selector: string) => {
      for (const root of roots()) await act(async () => (root.querySelector(selector) as HTMLElement | null)?.click());
    },
    /**
     * Each UI's menu opened by `selector`, one after the other, and compared: an open menu takes focus, so the one that
     * opens second closes the first (focus leaving a menu closes it, seqspec:interaction#slots-and-chrome), and both
     * are not open at once.
     */
    sameOpen: async (step: string, selector: string) => {
      const shapes: string[] = [];
      for (const root of roots()) {
        const trigger = () => root.querySelector(selector) as HTMLElement | null;
        await act(async () => trigger()?.click());
        await settle();
        shapes.push(shape(root).join("\n"));
        await act(async () => trigger()?.click());
        // Closed, focus is on its trigger: not left in one UI only, where it would keep that UI's box.
        await act(async () => (document.activeElement as HTMLElement | null)?.blur());
      }
      expect(`${step}\n${shapes[1]}`).toBe(`${step}\n${shapes[0]}`);
      return shapes[0];
    },
  };
}

/** The visible track entries' middles; `header` the header column's middle (null: no header column). */
function entries(ctx: SeqViewerContext, section: string): { header: number | null; x: number; y: number }[] {
  const l = ctx.layout.get(section)?.current;
  const vp = ctx.viewports.get(section)?.live.current;
  const out: { header: number | null; x: number; y: number }[] = [];
  for (let e = 0; l && vp && e < l.entryCount; e++) {
    if (l.entryKind[e] !== EntryKind.track || l.entryHeight[e] < 4) continue;
    const top = e < l.pinnedCount ? l.entryTop[e] : l.scrolledTop + l.entryTop[e] - vp.scrollY;
    const x = l.headerWidth + (vp.width - l.headerWidth) * 0.3;
    const header = l.headerWidth < 8 ? null : l.headerWidth / 2;
    if (top >= 0 && top < vp.height) out.push({ header, x, y: top + l.entryHeight[e] / 2 });
  }
  return out;
}

async function parity(name: string, doc: Parameters<SeqViewerContext["setDocument"]>[0]) {
  const ctx = await caseContext(doc);
  const ui = await both(ctx);
  const seen = new Set<string>();
  const step = async (s: string) => {
    for (const line of (await ui.same(`${name}: ${s}`)).split("\n")) seen.add(line.trim());
  };
  await step("initial");
  for (const v of ctx.views.list.value) {
    await act(async () => ctx.views.setCurrent(v.ref));
    await step(`view ${v.ref}`);
    const section = ctx.viewports.sections.value[0];
    const list = section ? entries(ctx, section) : [];
    const move = pointer(ctx);
    const h = list.find((e) => e.header !== null);
    if (h?.header) {
      await act(async () => move("move", section, h.header as number, h.y));
      await step(`view ${v.ref}, a header hovered`);
      for (const line of (
        await ui.sameOpen(`${name}: view ${v.ref}, its menu open`, ".seqv-track-actions .seqv-icon")
      ).split("\n"))
        seen.add(line.trim());
    }
    // The first entry with a label at 30% of its width.
    for (const e of list) {
      await act(async () => move("move", section, e.x, e.y));
      if (ctx.labels.current.value) break;
    }
    await step(`view ${v.ref}, a position hovered`);
    if (section) await act(async () => ctx.input.leave(ctx.layout.sectionIndex(section)));
  }
  for (const s of ctx.snapshots.list.value) {
    await act(async () => ctx.snapshots.apply(s.ref));
    await step(`snapshot ${s.ref}`);
  }
  await act(async () => ctx.report("seqv.render.labels", "a quiet note"));
  await step("a note alone, its quiet badge");
  await act(async () => ctx.report("seqv.test", "a reported problem"));
  await ui.click(".seqv-badge");
  await step("the issue list open, the note last under its heading");
  return seen;
}

const PARTS = [
  "seqv-tabs",
  "seqv-section-header",
  "seqv-chips",
  "seqv-track-actions",
  "seqv-menu-list",
  "seqv-tooltip",
].concat(["seqv-legend", "seqv-player", "seqv-description", "seqv-issue-list", "seqv-state", "seqv-switch"]);
const covered: string[] = [];
afterAll(() => console.info(`parity covered:\n${covered.join("\n")}`));

describe("React/DOM parity (U2)", () => {
  it("the hidden tracks' chips sit in the header's row, left of the switches, in its padding and the filter chip's look, in both UIs", async () => {
    const ctx = await caseContext(uiDoc());
    const ui = await both(ctx);
    await act(async () => {
      await ctx.actions.run("track.filterView", "t.dom");
      await ctx.actions.run("track.hide", "t.res");
    });
    await ui.same("a track hidden under a filter");
    for (const root of ui.roots()) {
      const header = root.querySelector(".seqv-section-header") as HTMLElement;
      const chips = header.querySelector(":scope > .seqv-chips") as HTMLElement;
      expect(chips?.textContent).toBe("Show Residues");
      // The only chips of the section sit in the header's row, no row of their own.
      expect(root.querySelectorAll(".seqv-chips").length).toBe(1);
      const [h, c] = [header, chips].map((e) => e.getBoundingClientRect());
      const padding = Number.parseFloat(getComputedStyle(header).paddingLeft);
      const filter = header.querySelector("[data-seqv-key=filter]") as HTMLElement;
      const first = header.querySelector(".seqv-switch, .seqv-switch-trigger") as HTMLElement;
      const f = filter.getBoundingClientRect();
      // Right of the filter chip, left of the first switch, one row, inside the header's padding.
      expect(c.left).toBeGreaterThanOrEqual(f.right);
      expect(c.right).toBeLessThanOrEqual(first.getBoundingClientRect().left + 0.5);
      expect(c.left).toBeGreaterThanOrEqual(h.left + padding);
      expect(c.top).toBeGreaterThanOrEqual(h.top - 0.5);
      expect(c.bottom).toBeLessThanOrEqual(h.bottom + 0.5);
      expect(Math.abs(c.top + c.height / 2 - (f.top + f.height / 2))).toBeLessThan(1);
      // The chip is the filter chip's: the same class, so the same border, radius and size.
      const chip = chips.querySelector("button") as HTMLElement;
      const look = (e: HTMLElement) => {
        const cs = getComputedStyle(e);
        return [e.className, cs.borderRadius, cs.borderTopColor, cs.fontSize];
      };
      expect(look(chip)).toEqual(look(filter));
    }
  });

  it("the track's menu: both UIs place it on the hovered line, a rows track's row on each hover", async () => {
    const ctx = await caseContext(uiDoc());
    const ui = await both(ctx);
    const move = pointer(ctx);
    const placed = () =>
      ui.roots().map((r) => {
        const b = r.querySelector(".seqv-track-actions") as HTMLElement | null;
        return b && [b.style.left, b.style.top, b.style.height, b.style.visibility];
      });
    const dom = lines(ctx, "chains.main", "t.dom")[0];
    await act(async () => move("move", "chains.main", dom.header, dom.y));
    const [r, d] = placed();
    expect(r?.[1]).toBe(`${dom.top}px`);
    expect(d).toEqual(r);
    await act(async () => ctx.input.leave(ctx.layout.sectionIndex("chains.main")));
    await act(async () => ctx.views.setCurrent("msa"));
    await settle();
    const rows = lines(ctx, "msa.main", "t.rows");
    for (const k of [2, 5, 0, 3]) {
      await act(async () => move("move", "msa.main", rows[k].header, rows[k].y));
      const [r, d] = placed();
      expect([k, r?.[1], r?.[2]]).toEqual([k, `${rows[k].top}px`, `${rows[k].height}px`]);
      expect(d).toEqual(r);
    }
  });
  it("a split track's rows: each UI's track menu toggles them, by a click or the keyboard, both menus following", async () => {
    const ctx = await caseContext(CASES.case13());
    const ui = await both(ctx);
    const move = pointer(ctx);
    const S = "features.main";
    const NV = "features.natural-variants";
    const t = ctx.graph.assembled.value?.view?.resolved.tracks.ref.indexOf(
      ctx.graph.assembled.value.set.strings.find(NV),
    );
    /** The track's entries: its label's, then a row per value when expanded. */
    const count = () => {
      const l = ctx.layout.get(S)?.current;
      return l ? l.entryTrack.subarray(0, l.entryCount).filter((x) => x === t).length : 0;
    };
    const item = (r: Element) => r.querySelector(".seqv-menu-item") as HTMLElement | null;
    /** Hovers the track's first line. */
    const hover = async () => {
      const l = ctx.layout.get(S)?.current;
      const e = l ? l.entryTrack.indexOf(t ?? -1) : -1;
      await act(async () => ctx.viewports.get(S)?.scrollTo(0, l ? l.entryTop[e] : 0, "api"));
      await settle();
      const line = lines(ctx, S, NV)[0];
      await act(async () => move("move", S, line.header, line.y));
    };
    /** One UI's menu open (an open menu has focus, and focus leaving a menu closes it: the other's closes). */
    const open = async (r: Element) => {
      await hover();
      if (!item(r)) await act(async () => (r.querySelector(".seqv-track-actions .seqv-icon") as HTMLElement).click());
      await settle();
      return item(r)?.textContent;
    };
    await settle();
    const expanded = count();
    expect(expanded).toBeGreaterThan(10);
    // Hovering lane k (the label's line is first): both UIs put the "…" on lane k's line, the hover names lane k.
    await hover();
    const lanes = lines(ctx, S, NV);
    for (const k of [0, 7, 15, lanes.length - 2]) {
      const line = lanes[k + 1];
      await act(async () => move("move", S, line.header, line.y));
      await settle();
      expect(ctx.hover.live.current.split).toBe(k);
      const tops = ui.roots().map((r) => (r.querySelector(".seqv-track-actions") as HTMLElement | null)?.style.top);
      expect([k, ...tops]).toEqual([k, `${line.top}px`, `${line.top}px`]);
    }
    const [react, dom] = ui.roots();
    expect([await open(react), await open(dom)]).toEqual(["Collapse rows", "Collapse rows"]);
    await ui.sameOpen("case 13: a split track's menu open", ".seqv-track-actions .seqv-icon");
    // The DOM's item by a click: one compact row; the React's by the keyboard: the rows again.
    await open(dom);
    await act(async () => item(dom)?.click());
    await settle();
    expect(count()).toBe(1);
    expect([await open(react), await open(dom)]).toEqual(["Expand rows", "Expand rows"]);
    await open(react);
    item(react)?.focus();
    await act(async () => userEvent.keyboard("{Enter}"));
    await settle();
    expect(count()).toBe(expanded);
  });
  it("the legend: both UIs start its items at the section's data left", async () => {
    const ctx = await caseContext(uiDoc());
    const ui = await both(ctx);
    for (const headerWidth of [150, 200]) {
      await act(async () => {
        ctx.setRenderProps({ headerWidth });
        ctx.layout.get("chains.main"); // As the canvas reads it on its next frame (the layout builds on read).
      });
      await settle();
      const lefts = ui.roots().map((r) => {
        const x0 = (r.querySelector(".seqv-surface") as HTMLElement).getBoundingClientRect().left;
        return [...r.querySelectorAll(".seqv-legend-group")].map((g) => [
          g.getAttribute("data-section"),
          // The first item, or a group's color scheme switch (seqspec:interaction#color-schemes-at-runtime): from the
          // data left either way.
          (
            (g.querySelector(".seqv-legend-items") as HTMLElement).firstElementChild as HTMLElement
          ).getBoundingClientRect().left - x0,
        ]);
      });
      expect(lefts[0].length).toBeGreaterThan(0);
      expect(lefts[1]).toEqual(lefts[0]);
      for (const [section, x] of lefts[0]) expect(x).toBe(ctx.layout.get(String(section))?.current.dataLeft);
    }
  });
  it("a compact legend: 25 tracks of three kinds are three groups, the same in both UIs, closed and open", async () => {
    const ctx = await caseContext(kindsDoc());
    const ui = await both(ctx);
    await ui.same("closed");
    for (const r of ui.roots()) {
      const names = [...r.querySelectorAll(".seqv-legend-group .seqv-legend-name")].map((x) => x.textContent);
      expect(names).toEqual(["Genes (5)", "Transcripts (10)", "Domains (10)"]);
    }
    const sel = '[data-seqv-key="colors:v.main/tx1.1"]';
    for (const r of ui.roots()) {
      const t = r.querySelector(sel) as HTMLElement;
      const icon = t.querySelector(".seqv-legend-icon");
      // The title is the button: the group's name whole, then the palette icon labeled with the scheme shown.
      expect([t.tagName, t.getAttribute("aria-haspopup"), t.textContent]).toEqual([
        "BUTTON",
        "menu",
        "Transcripts (10)",
      ]);
      expect([icon?.getAttribute("role"), icon?.getAttribute("aria-label")]).toEqual(["img", "Colors: document"]);
      // A hint (the full name) only when the name is cut: here it is whole.
      expect(t.closest(".seqv-hint")).toBeNull();
    }
    expect(await ui.sameOpen("open", sel)).toContain("Okabe–Ito (colorblind-safe)");
  });
  it("a group's title is never cut while the legend has room: a long name gets its own line; cut past it, its hint is the full name", async () => {
    const ctx = await caseContext(kindsDoc());
    const ui = await both(ctx);
    const sel = '[data-seqv-key="colors:v.main/tx1.1"]';
    const state = (r: Element) => {
      const g = (r.querySelector(sel) as HTMLElement).closest(".seqv-legend-group") as HTMLElement;
      const name = g.querySelector(".seqv-legend-name") as HTMLElement;
      return [
        g.hasAttribute("data-seqv-long"),
        g.hasAttribute("data-seqv-cut"),
        name.scrollWidth <= name.clientWidth,
        g.querySelector(".seqv-hint-bubble")?.textContent ?? null,
      ];
    };
    // A label column narrower than "Transcripts (10)": the title on its own line, whole.
    await act(async () => {
      ctx.setRenderProps({ headerWidth: 60 });
      ctx.layout.get("v.main");
    });
    await settle();
    for (const r of ui.roots()) expect(state(r)).toEqual([true, false, true, null]);
    // A legend narrower than the name: cut, and its hint shows the full name.
    for (const r of ui.roots()) (r.parentElement as HTMLElement).style.width = "70px";
    await settle();
    await new Promise((ok) => requestAnimationFrame(() => requestAnimationFrame(ok)));
    for (const r of ui.roots()) expect(state(r)).toEqual([true, true, false, "Transcripts (10)"]);
    await ui.same("cut: the title in its hint");
  });
  it("a group's title by the keyboard, in each UI: Enter or Space opens, the arrows move, Escape closes, a pick goes to every track", async () => {
    const ctx = await caseContext(kindsDoc());
    const ui = await both(ctx);
    const sel = '[data-seqv-key="colors:v.main/tx1.1"]';
    const layers = Array.from({ length: 10 }, (_, k) => `l.tx${1 + (k >> 1)}.${1 + (k & 1)}`);
    for (const [k, r] of ui.roots().entries()) {
      const trigger = () => r.querySelector(sel) as HTMLElement;
      const focused = () => document.activeElement as HTMLElement | null;
      trigger().focus();
      await act(async () => userEvent.keyboard(k === 0 ? "{Enter}" : " "));
      await settle();
      expect([trigger().getAttribute("aria-expanded"), focused()?.textContent]).toEqual(["true", "document"]);
      expect(focused()?.getAttribute("aria-checked")).toBe("true");
      await act(async () => userEvent.keyboard("{ArrowDown}"));
      expect(focused()?.textContent).toBe("Category 10");
      await act(async () => userEvent.keyboard("{Escape}"));
      await settle();
      expect([r.querySelector("[role=menu]"), focused()]).toEqual([null, trigger()]);
      await act(async () => userEvent.keyboard("{Enter}"));
      await settle();
      await act(async () => userEvent.keyboard("{ArrowDown}{ArrowDown}{Enter}"));
      await act(async () => ctx.graph.update());
      await settle();
      // Every transcript's layer; the icon's label names it; focus stays on the same group's title.
      expect(layers.map((l) => ctx.ui.colorSchemes.value.get(l))).toEqual(layers.map(() => "okabe_ito"));
      expect([trigger().querySelector(".seqv-legend-icon")?.getAttribute("aria-label"), focused()]).toEqual([
        "Colors: Okabe–Ito",
        trigger(),
      ]);
      await act(async () => {
        for (const l of layers) ctx.ui.setColorScheme(l, null);
      });
      await settle();
    }
    await ui.same("reset");
  });
  it("a menu is in the top layer: nothing clips it, at the legend's bottom edge and in a narrow viewport, in both UIs", async () => {
    const ctx = await caseContext(kindsDoc());
    const ui = await both(ctx);
    const sel = '[data-seqv-key="colors:v.main/dom1.1"]';
    // Each root in a box that clips (as a card or a scroll area would), the legend's last group at its bottom edge.
    for (const r of ui.roots()) {
      const box = r.parentElement as HTMLElement;
      box.style.overflow = "hidden";
      box.style.height = `${(r.querySelector(sel) as HTMLElement).getBoundingClientRect().bottom - box.getBoundingClientRect().top + 2}px`;
    }
    const check = async (r: Element) => {
      const trigger = r.querySelector(sel) as HTMLElement;
      trigger.scrollIntoView({ block: "end" });
      await act(async () => trigger.click());
      await settle();
      const list = r.querySelector("[role=menu]") as HTMLElement;
      const box = list.getBoundingClientRect();
      const t = trigger.getBoundingClientRect();
      const [vw, vh] = [document.documentElement.clientWidth, document.documentElement.clientHeight];
      // Shown whole: in the top layer, inside the viewport, and what is at its corners is the menu's.
      expect(list.matches(":popover-open")).toBe(true);
      expect([box.left >= 0, box.top >= 0, box.right <= vw, box.bottom <= vh]).toEqual([true, true, true, true]);
      for (const [x, y] of [
        [box.left + 4, box.top + 4],
        [box.right - 4, box.bottom - 4],
      ])
        expect(list.contains(document.elementFromPoint(x, y))).toBe(true);
      // Below the trigger, or flipped above it where the viewport ends.
      expect(box.top >= t.bottom || box.bottom <= t.top).toBe(true);
      await act(async () => trigger.click());
      await act(async () => (document.activeElement as HTMLElement | null)?.blur());
      return box.top >= t.bottom ? "below" : "above";
    };
    for (const r of ui.roots()) await check(r);
    // A narrow, short viewport: the menu flips above and keeps within the right edge.
    await page.viewport(320, 260);
    try {
      for (const r of ui.roots()) expect(await check(r)).toBe("above");
    } finally {
      await page.viewport(1200, 900);
    }
  });
  it("a track header's tooltip: its own legend under its line, never over its ⋯, which hides it; by the keyboard too", async () => {
    const ctx = await caseContext(kindsDoc());
    const ui = await both(ctx);
    const move = pointer(ctx);
    const line = lines(ctx, "v.main", "tx1.2")[0];
    await act(async () => move("move", "v.main", line.header, line.y));
    await settle();
    const legendOf = (r: Element) => {
      const tip = r.querySelector(".seqv-tooltip") as HTMLElement;
      return [
        tip.querySelector(".seqv-tooltip-title")?.textContent,
        [...tip.querySelectorAll(".seqv-tooltip-legend .seqv-legend-item")].map((x) => x.textContent),
      ];
    };
    for (const r of ui.roots()) {
      expect(legendOf(r)).toEqual(["Transcript G1.2", ["UTR", "exon"]]);
      // The ⋯ is on the hovered line, the tooltip below it: what is at the button's middle is the button.
      const b = r.querySelector(".seqv-track-actions .seqv-icon") as HTMLElement;
      const at = b.getBoundingClientRect();
      const tip = (r.querySelector(".seqv-tooltip") as HTMLElement).getBoundingClientRect();
      expect(tip.top >= at.bottom || tip.bottom <= at.top).toBe(true);
      expect(document.elementFromPoint(at.left + at.width / 2, at.top + at.height / 2)).toBe(b);
    }
    await ui.same("a header hovered: its legend");
    // The pointer onto the ⋯: the tooltip gives way, and the menu opens.
    for (const r of ui.roots()) {
      const b = r.querySelector(".seqv-track-actions .seqv-icon") as HTMLElement;
      await userEvent.hover(b);
      expect(getComputedStyle(r.querySelector(".seqv-tooltip") as HTMLElement).visibility).toBe("hidden");
      await act(async () => b.click());
      await settle();
      expect(r.querySelector("[role=menu]")).not.toBeNull();
      await act(async () => b.click());
      await act(async () => (document.activeElement as HTMLElement | null)?.blur());
    }
    await userEvent.unhover(ui.roots()[1].querySelector(".seqv-track-actions .seqv-icon") as HTMLElement);
    await act(async () => ctx.input.leave(ctx.layout.sectionIndex("v.main")));
    await settle();
    for (const r of ui.roots()) expect(r.querySelector(".seqv-tooltip")).toBeNull();
    // The keyboard cursor on a header: Home goes to it, Down moves along the headers.
    const key = keys(ctx);
    await act(async () => {
      ctx.input.focus(ctx.layout.sectionIndex("v.main"), true);
      key("v.main", "Home");
      for (let i = 0; i < 30 && ctx.ui.keyCursor.value?.track !== "dom2.1"; i++) key("v.main", "ArrowDown");
    });
    await settle();
    expect(ctx.ui.keyCursor.value?.header).toBe(true);
    for (const r of ui.roots()) expect(legendOf(r)).toEqual(["G2.1 domain", ["domain"]]);
    await ui.same("the keyboard cursor on a header: its legend");
  });
  it("a header the canvas cuts: its tooltip's first line is the full text; whole and with no legend, no tooltip", async () => {
    const ctx = await caseContext(cutDoc());
    const ui = await both(ctx);
    const move = pointer(ctx);
    const title = async (track: string, k: number) => {
      const line = lines(ctx, "v.main", track)[k];
      await act(async () => move("move", "v.main", 12, line.y));
      await settle();
      return ui.roots().map((r) => r.querySelector(".seqv-tooltip .seqv-tooltip-title")?.textContent ?? null);
    };
    // A row's name and a track's label cut by the label column: each whole in its tooltip.
    expect(await title("v.rows", 0)).toEqual([LONG, LONG]);
    await ui.same("a cut row's name in its tooltip");
    expect(await title("v.letters", 0)).toEqual([LONG_TRACK, LONG_TRACK]);
    // Whole, and the rows list nothing: no tooltip.
    expect(await title("v.rows", 1)).toEqual([null, null]);
  });
  it("the legend: collapsed by default, opened and closed by `seqv:legend`, the toggle kept across a patch", async () => {
    const ctx = createSeqViewer(fullChrome(), { scheduler: manualScheduler(), text: fixtureMeasurer() });
    await ctx.setDocument(kindsDoc());
    const ui = await both(ctx);
    const open = () => ui.roots().map((r) => (r.querySelector(".seqv-legend") as HTMLDetailsElement).open);
    expect(open()).toEqual([false, false]);
    await ui.same("collapsed");
    const run = (o: unknown) =>
      act(async () => void (await ctx.commands.run({ command: "seqv:legend", params: { open: o } }, "host:links")));
    await run(true);
    expect([open(), ctx.ui.legendExpanded.value]).toEqual([[true, true], true]);
    await ui.same("opened by the command");
    await run("toggle");
    expect(open()).toEqual([false, false]);
    await run("sideways"); // refused: nothing changes
    expect([open(), ctx.ui.legendExpanded.value]).toEqual([[false, false], false]);
    // The person's toggle in one UI: both follow, and a patch, then a new document, keep it.
    await act(async () => (ui.roots()[1].querySelector(".seqv-legend > summary") as HTMLElement).click());
    await settle();
    expect(open()).toEqual([true, true]);
    const ops = [{ op: "set" as const, ref: "tx1.1", params: { label: "T" } }];
    expect((await ctx.applyPatch({ format: "seq_view_spec_patch", ops })).ok).toBe(true);
    await settle();
    expect(open()).toEqual([true, true]);
    await act(async () => void (await ctx.setDocument(kindsDoc())));
    await settle();
    expect(open()).toEqual([true, true]);
  });
  it("the legend's `hidden` option leaves the part out of either UI; `UI.Legend` sets every UI's", async () => {
    const ctx = createSeqViewer(fullChrome(), { scheduler: manualScheduler(), text: fixtureMeasurer() });
    await ctx.setDocument(kindsDoc());
    const react = await mount(ctx, { surface: stubSurface, legend: "hidden" });
    const host = document.createElement("div");
    document.body.appendChild(host);
    const m = mountSeqViewer(host, ctx, { surface: stubSurface, legend: "expanded" });
    cleanups.push(react.unmount, () => {
      m.dispose();
      host.remove();
    });
    await settle();
    const dom = host.querySelector(".seqv-legend") as HTMLDetailsElement;
    expect([react.container.querySelector(".seqv-legend"), dom.open]).toEqual([null, true]);
    ctx.config.set(SeqViewerConfig.UI.Legend, "hidden");
    const again = document.createElement("div");
    document.body.appendChild(again);
    const n = mountSeqViewer(again, ctx, { surface: stubSurface });
    cleanups.push(() => {
      n.dispose();
      again.remove();
      ctx.dispose();
    });
    await settle();
    expect(again.querySelector(".seqv-legend")).toBeNull();
  });
  it("a track with no legend row: no group, and its menu's Colors ▸ lists its schemes in place, in both UIs", async () => {
    const ctx = await caseContext(kindsDoc(true));
    const ui = await both(ctx);
    for (const r of ui.roots())
      expect([...r.querySelectorAll(".seqv-legend-name")].map((x) => x.textContent)).toEqual([
        "Transcripts (10)",
        "Domains (10)",
      ]);
    const move = pointer(ctx);
    const items = (r: Element) => [...r.querySelectorAll<HTMLElement>("[role^=menuitem]")];
    const open = async (r: Element) => {
      const line = lines(ctx, "v.main", "gene2")[0];
      await act(async () => move("move", "v.main", line.header, line.y));
      await act(async () => (r.querySelector(".seqv-track-actions .seqv-icon") as HTMLElement).click());
      await settle();
    };
    const shapes: string[] = [];
    for (const [k, r] of ui.roots().entries()) {
      await open(r);
      const colors = items(r).find((i) => i.textContent === "Colors ▸") as HTMLElement;
      expect(colors).toBeDefined();
      colors.focus();
      await act(async () => (k === 0 ? colors.click() : userEvent.keyboard("{Enter}")));
      await settle();
      // The same menu, now the schemes after a way back, focus on its first item.
      expect(items(r).map((i) => i.textContent)).toEqual([
        "‹ Colors",
        "document",
        "Category 10",
        "Okabe–Ito (colorblind-safe)",
        "Tableau 10",
      ]);
      expect(document.activeElement).toBe(items(r)[0]);
      shapes.push(shape(r.querySelector("[role=menu]") as Element).join("\n"));
      // Back, and the menu again; closed, it opens at its top.
      await act(async () => items(r)[0].click());
      await settle();
      expect(items(r)[0].textContent).toBe("Show only what this track covers");
      await act(async () => userEvent.keyboard("{Escape}"));
      await settle();
      await open(r);
      expect(items(r).map((i) => i.textContent)).toContain("Colors ▸");
      await act(async () => (items(r).find((i) => i.textContent === "Colors ▸") as HTMLElement).click());
      await settle();
      await act(async () => (items(r).find((i) => i.textContent === "Tableau 10") as HTMLElement).click());
      await settle();
      expect([r.querySelector("[role=menu]"), ctx.ui.colorSchemes.value.get("l.gene2")]).toEqual([null, "tableau10"]);
      await act(async () => ctx.ui.setColorScheme("l.gene2", null));
      await act(async () => (document.activeElement as HTMLElement | null)?.blur());
    }
    expect(shapes[1]).toBe(shapes[0]);
  });
  it('the legend past `UI.LegendGroups`: "+n more" in both UIs, and either opens both', async () => {
    const ctx = await caseContext(distinctDoc(12));
    const ui = await both(ctx);
    const count = () => ui.roots().map((r) => r.querySelectorAll(".seqv-legend-track:not(:empty)").length);
    await ui.same("closed");
    expect(count()).toEqual([8, 8]);
    const more = '[data-seqv-key="more:v.main/*"]';
    for (const r of ui.roots()) expect(r.querySelector(more)?.textContent).toBe("+4 more");
    await act(async () => (ui.roots()[1].querySelector(more) as HTMLElement).click());
    await ui.same("open");
    expect(count()).toEqual([12, 12]);
    for (const r of ui.roots())
      expect([r.querySelector(more)?.textContent, r.querySelector(more)?.getAttribute("aria-expanded")]).toEqual([
        "fewer",
        "true",
      ]);
  });
  it("the legend in a wrapped section: a group is one row, the name then its items, the next rows under the first item", async () => {
    // Case 1's first view is wrapped text: no label column (`dataLeft` the text's inset), the Regions track's legend.
    const ctx = await caseContext(CASES.case1());
    const ui = await both(ctx);
    await act(async () => {
      ctx.layout.get("text.main");
    });
    await settle();
    const l = ctx.layout.get("text.main")?.current;
    expect(l?.geometry).toBe(AxisGeometry.wrapped);
    for (const r of ui.roots()) {
      const g = [...r.querySelectorAll('.seqv-legend-group[data-section="text.main"]')].find(
        (x) => x.querySelector(".seqv-legend-name")?.textContent === "Regions",
      ) as HTMLElement;
      const name = (g.querySelector(".seqv-legend-name") as HTMLElement).getBoundingClientRect();
      const cell = (g.querySelector(".seqv-legend-track") as HTMLElement).getBoundingClientRect();
      const swatch = (g.querySelector(".seqv-swatch") as HTMLElement).getBoundingClientRect();
      const x0 = (r.querySelector(".seqv-surface") as HTMLElement).getBoundingClientRect().left;
      // One row: the name's text at the data area's left (or the label inset, the larger), the first swatch after it,
      // both centered on the row.
      expect(cell.left - x0).toBe(0);
      expect(name.left - x0).toBeGreaterThanOrEqual(l?.dataLeft ?? 0);
      expect(name.left - x0).toBeLessThan((l?.dataLeft ?? 0) + 12);
      expect(name.width).toBeGreaterThan(30);
      expect(swatch.left).toBeGreaterThanOrEqual(name.right);
      expect(Math.abs(name.top + name.height / 2 - (swatch.top + swatch.height / 2))).toBeLessThan(3);
      expect(g.getBoundingClientRect().height).toBeLessThan(cell.height * 1.5);
    }
  });
  it("a wrapped legend group's many items wrap onto further rows under the first item, not under the name", async () => {
    const b = createSeqViewSpec({ title: "Many" });
    const seq = b.data.sequence({ ref: "s", value: "ACDEFGHIKLMNPQRSTVWY".repeat(6) });
    const regions = b.data.annotation({
      ref: "regions",
      geometry: "ranges",
      target: seq,
      value: Array.from({ length: 12 }, (_, i) => ({ start: i * 8, stop: i * 8 + 6, label: `Region number ${i + 1}` })),
    });
    const text = b.view({ ref: "text", label: "Text", layout: "text" });
    text.track({ label: "Sequence" }).letters();
    text.track({ label: "Regions" }).ranges({
      data: regions,
      display: "underline",
      color: { theme: "field_palette", params: { field: "label", palette: { kind: "categorical", colors: "Set2" } } },
    });
    const ctx = await caseContext(b.getDocument());
    const ui = await both(ctx);
    await act(async () => {
      ctx.layout.get("text.main");
    });
    await settle();
    expect(ctx.layout.get("text.main")?.current.geometry).toBe(AxisGeometry.wrapped);
    for (const r of ui.roots()) {
      const g = [...r.querySelectorAll('.seqv-legend-group[data-section="text.main"]')].find(
        (x) => x.querySelector(".seqv-legend-name")?.textContent === "Regions",
      ) as HTMLElement;
      const items = [...g.querySelectorAll(".seqv-legend-item")].map((x) => x.getBoundingClientRect());
      const name = (g.querySelector(".seqv-legend-name") as HTMLElement).getBoundingClientRect();
      const tops = [...new Set(items.map((i) => Math.round(i.top)))];
      expect(tops.length).toBeGreaterThan(1);
      // The first of each row starts at the first item's left, which is after the name.
      expect(items[0].left).toBeGreaterThanOrEqual(name.right);
      for (const t of tops.slice(1)) {
        const first = items.find((i) => Math.round(i.top) === t) as DOMRect;
        expect(first.left).toBe(items[0].left);
      }
      // The name sits on the first row only.
      expect(Math.abs(name.top + name.height / 2 - (items[0].top + items[0].height / 2))).toBeLessThan(3);
    }
  });
  it("each chrome part placed by a host outside the layout: the same in both UIs, with the UI's scope and without", async () => {
    const ctx = await caseContext(uiDoc());
    cleanups.push(() => ctx.dispose());
    // The context offers the SVG export, as `SeqViewer` and the element register it.
    ctx.register({ actions: [D.browserExportAction({ offer: () => {} })] });
    const S = { section: "chains.main" };
    const search = async () => ({ results: [] });
    const react = (c: typeof ctx): ReactNode => (
      <>
        <R.ViewTabs context={c} />
        <R.SnapshotPlayer context={c} />
        <R.Description context={c} />
        <R.SectionHeader context={c} {...S} />
        <R.SectionSwitches context={c} {...S} />
        <R.HiddenChips context={c} {...S} />
        <R.TrackActions context={c} {...S} track="t.dom" />
        <R.Legend context={c} />
        <R.Search context={c} search={search} />
        <R.Busy context={c} />
        <R.Export context={c} />
        <R.IssueBadge context={c} />
        <R.EmptyState context={c} />
        <R.ErrorState context={c} />
        <R.Tooltip context={c} />
        <R.HoverCard context={c} />
      </>
    );
    const dom = (host: HTMLElement, ui?: D.DomUI) => [
      D.ViewTabs(host, ctx, { ui }),
      D.SnapshotPlayer(host, ctx, { ui }),
      D.Description(host, ctx, { ui }),
      D.SectionHeader(host, ctx, { ui, ...S }),
      D.SectionSwitches(host, ctx, { ui, ...S }),
      D.HiddenChips(host, ctx, { ui, ...S }),
      D.TrackActions(host, ctx, { ui, ...S, track: "t.dom" }),
      D.Legend(host, ctx, { ui }),
      D.Search(host, ctx, { ui, search }),
      D.Busy(host, ctx, { ui }),
      D.Export(host, ctx, { ui }),
      D.IssueBadge(host, ctx, { ui }),
      D.EmptyState(host, ctx, { ui }),
      D.ErrorState(host, ctx, { ui }),
      D.Tooltip(host, ctx, { ui }),
      D.HoverCardPart(host, ctx, { ui }),
    ];
    await act(async () => {
      await ctx.actions.run("track.hide", "t.res");
      ctx.report("seqv.test", "a reported problem");
    });
    await settle();
    const before = subscriptions(ctx);
    for (const scoped of [true, false]) {
      const [r, d] = [document.createElement("div"), document.createElement("div")];
      document.body.append(r, d);
      const root = createRoot(r);
      const parts = react(ctx);
      await act(async () =>
        root.render(
          scoped ? (
            <R.SeqViewerProvider context={ctx} ui={{ surface: stubSurface }}>
              {parts}
            </R.SeqViewerProvider>
          ) : (
            parts
          ),
        ),
      );
      const handles = dom(d, scoped ? D.createDomUI(ctx, { surface: stubSurface }) : undefined);
      await settle();
      const [a, b] = [r, d].map((e) => shape(e).join("\n"));
      expect(`${scoped}\n${b}`).toBe(`${scoped}\n${a}`);
      for (const c of [
        "seqv-tabs",
        "seqv-player",
        "seqv-section-header",
        "seqv-switches",
        "seqv-chips",
        "seqv-legend",
        "seqv-search",
        "seqv-export",
        "seqv-issues",
        '"track-actions"',
      ])
        expect([scoped, c, a.includes(c)]).toEqual([scoped, c, true]);
      await act(async () => root.unmount());
      for (const h of handles) h.dispose();
      expect(d.childNodes.length).toBe(0);
      // Every part, scoped or not, leaves no subscription behind.
      expect([scoped, subscriptions(ctx)]).toEqual([scoped, before]);
      r.remove();
      d.remove();
    }
  });

  it("the UI fixture", async () => {
    const seen = await parity("ui", uiDoc());
    for (const c of ["seqv-tabs", "seqv-player", "seqv-description", "seqv-legend", "seqv-tooltip", "seqv-menu-list"])
      expect([...seen].some((s) => s.includes(c))).toBe(true);
  });
  for (const name of Object.keys(CASES))
    it(`case ${name.slice(4)}`, async () => {
      const seen = await parity(name, CASES[name]());
      expect([...seen].some((s) => s.includes("seqv-section"))).toBe(true);
      const parts = PARTS.filter((p) => [...seen].some((s) => s.includes(`"${p}`) || s.includes(` ${p}`)));
      covered.push(`${name}: ${parts.join(" ")}`);
    });
});
