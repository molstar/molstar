/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The `Search` part (LS2; seqspec:interaction#the-default-ui, seqspec:interaction#slots-and-chrome) in both UIs on one
// context: absent without the UI spec's `search`; with it, a box that calls `search(text, signal)` once the text
// settles (the previous call aborted), lights the answer through `ctx.results` ("search"), steps with Enter and
// Shift+Enter, and saves the results as a location set. The two UIs give the same structure and text. The status
// (`searching…`, `no results`, `2 / 14`, the error) sits in a slot of a fixed width inside the field, so no state moves
// the field or the buttons; and both parts take a `search` of their own, so the box stands anywhere with a context.

import { act } from "react";
import { createRoot } from "react-dom/client";
import { afterEach, describe, expect, it } from "vitest";
import { userEvent } from "vitest/browser";
import { Search as DomSearch, mountSeqViewer } from "../../src/dom/index.ts";
import type { Result, SearchFn, SeqViewerContext } from "../../src/index.ts";
import { Search as ReactSearch } from "../../src/react/index.ts";
import "../../src/skin.css";

import { SEARCH_SETTLE_MS } from "../../src/search-box.ts";
import { openedContext } from "../dom/fixtures.ts";
import { mount, settle, stubSurface, WIDTH } from "./fixtures.tsx";

const PROTEIN = "MSEEKLKELAEKAIRLAKEGNPEEALKLLEEALRLAKESGDPELIARVLELLAEVTKQSVHAYRLDDGKWFLQPGGRVSLEMAPTDLQAAH";
const cleanups: (() => Promise<void> | void)[] = [];
afterEach(async () => {
  for (const c of cleanups.splice(0)) await c();
});

/** The application's search: every occurrence of the text in the kinase, as its residues. */
function searchOf(ctx: SeqViewerContext, asked: string[]): SearchFn {
  return async (text, signal) => {
    asked.push(text);
    const results: Result[] = [];
    for (let i = PROTEIN.indexOf(text); i >= 0; i = PROTEIN.indexOf(text, i + 1))
      results.push(ctx.refs.toLoci("a", { start: i, stop: i + text.length }) as Result);
    if (signal.aborted) throw new Error("aborted");
    if (text === "boom") throw new Error("The search service is not reachable at the moment");
    return { results, labels: results.map((_, k) => `hit ${k + 1}`) };
  };
}

async function both(search?: (ctx: SeqViewerContext) => SearchFn) {
  const { ctx, advance } = await openedContext();
  const react = await mount(ctx, { surface: stubSurface, search: search?.(ctx) });
  const dom = document.createElement("div");
  dom.style.width = `${WIDTH}px`;
  document.body.appendChild(dom);
  const m = mountSeqViewer(dom, ctx, { surface: stubSurface, search: search?.(ctx) });
  cleanups.push(react.unmount, () => {
    m.dispose();
    dom.remove();
    ctx.dispose();
  });
  const roots = [react.container, dom].map((c) => c.firstElementChild as HTMLElement);
  const text = () => roots.map((r) => r.querySelector(".seqv-search")?.textContent ?? null);
  const settled = async () => {
    await act(async () => advance(SEARCH_SETTLE_MS));
    await settle();
  };
  return { ctx, roots, text, settled };
}

describe("the Search part (LS2)", () => {
  it("is left out of both UIs when the UI spec gives no `search`", async () => {
    const t = await both();
    expect(t.text()).toEqual([null, null]);
  });

  it("searches as the text settles, steps with Enter, and keeps the results, the same in both UIs", async () => {
    const asked: string[] = [];
    const t = await both((ctx) => searchOf(ctx, asked));
    const inputs = t.roots.map((r) => r.querySelector(".seqv-search-input") as HTMLInputElement);
    for (const input of inputs) await userEvent.type(input, "EAL");
    await t.settled();
    // Each UI's box called once, with the settled text only.
    expect(asked).toEqual(["EAL", "EAL"]);
    expect(t.text()).toEqual(
      Array(2).fill(
        "– / 2‹Previous result›Next resultKeepKeep these results highlightedRemove every kept search highlight",
      ),
    );
    expect(t.ctx.highlight.channels.value.get("results:search")?.style).toBe("highlight");

    await userEvent.type(inputs[0], "{Enter}");
    await settle();
    expect(t.text().map((s) => s?.slice(0, 5))).toEqual(["1 / 2", "1 / 2"]);
    expect(t.ctx.focus.current.value.get("chains.main")).toMatchObject({ start: 23, stop: 26 });
    await userEvent.type(inputs[1], "{Shift>}{Enter}{/Shift}");
    await settle();
    expect(t.text().map((s) => s?.slice(0, 5))).toEqual(["2 / 2", "2 / 2"]);
    expect(t.ctx.highlight.channels.value.get("results:search:current")?.style).toBe("accent");

    await act(async () => (t.roots[1].querySelector('[data-seqv-key="save"]') as HTMLElement).click());
    await settle();
    const data = t.ctx.document.shown.value?.root.children?.find((c) => c.kind === "data");
    const saved = data?.children?.find((c) => c.ref === "chains.main.keep.eal");
    expect(saved?.params?.value).toEqual([
      { target: "a/unit", start: 23, stop: 26, id: "r1", label: "hit 1" },
      { target: "a/unit", start: 30, stop: 33, id: "r2", label: "hit 2" },
    ]);

    // Cleared: the results go, and the box says nothing.
    await userEvent.clear(inputs[0]);
    await settle();
    expect(t.ctx.highlight.channels.value.has("results:search")).toBe(false);
  });
  it("lists a kept set under Highlights in the legend, replaces it when the same query is kept again, and removes it (Undo returns it)", async () => {
    const t = await both((ctx) => searchOf(ctx, []));
    const inputs = t.roots.map((r) => r.querySelector(".seqv-search-input") as HTMLInputElement);
    const keep = (i: number) =>
      act(async () => (t.roots[i].querySelector('[data-seqv-key="save"]') as HTMLElement).click());
    const kept = () => {
      const data = t.ctx.document.shown.value?.root.children?.find((c) => c.kind === "data");
      return (data?.children ?? []).filter((c) => c.ref?.startsWith("chains.main.keep.")).map((c) => c.ref);
    };
    const highlights = () =>
      t.ctx.document.shown.value?.root.children?.flatMap((c) => c.children ?? []).find((c) => c.ref === "chains.main")
        ?.params?.highlights;
    const rows = () =>
      t.roots.map((r) =>
        [...r.querySelectorAll(".seqv-legend-group")].find(
          (g) => g.querySelector(".seqv-legend-name")?.textContent === "Highlights",
        ),
      );
    for (const input of inputs) await userEvent.type(input, "EAL");
    await t.settled();
    expect(kept()).toEqual([]);
    expect(rows()).toEqual([undefined, undefined]);
    await keep(0);
    await settle();
    expect(kept()).toEqual(["chains.main.keep.eal"]);
    expect(highlights()).toEqual([{ data: "chains.main.keep.eal" }]);
    // Both legends list it by the query, with a ✕; the same structure.
    expect(rows().map((g) => g?.textContent)).toEqual(Array(2).fill("HighlightsEAL✕Remove EAL from the highlights"));
    expect(rows().map((g) => g?.querySelector(".seqv-swatch")?.getAttribute("data-style"))).toEqual(["band", "band"]);
    // The same query kept again replaces the set; another query is another set.
    await keep(1);
    await settle();
    expect(kept()).toEqual(["chains.main.keep.eal"]);
    expect((highlights() as unknown[]).length).toBe(1);
    for (const input of inputs) {
      await userEvent.clear(input);
      await userEvent.type(input, "KL");
    }
    await t.settled();
    await keep(0);
    await settle();
    expect(kept()).toEqual(["chains.main.keep.eal", "chains.main.keep.kl"]);
    // ✕ in either UI removes the entry and its set in one patch; Undo returns both.
    await act(async () => (t.roots[1].querySelector('[data-seqv-key="hl:chains.main/0"]') as HTMLElement).click());
    await settle();
    expect(kept()).toEqual(["chains.main.keep.kl"]);
    expect(highlights()).toEqual([{ data: "chains.main.keep.kl" }]);
    expect(rows().map((g) => g?.querySelectorAll(".seqv-legend-item").length)).toEqual([1, 1]);
    await act(async () => {
      await t.ctx.history.undo();
    });
    await settle();
    expect(kept()).toEqual(["chains.main.keep.eal", "chains.main.keep.kl"]);
    expect(rows().map((g) => g?.querySelectorAll(".seqv-legend-item").length)).toEqual([2, 2]);
  });

  it('clears what it kept where it was kept: "Clear kept" removes every kept set, by the keyboard too, in both UIs', async () => {
    const t = await both((ctx) => searchOf(ctx, []));
    const inputs = t.roots.map((r) => r.querySelector(".seqv-search-input") as HTMLInputElement);
    const clears = () => t.roots.map((r) => r.querySelector('[data-seqv-key="clear-kept"]') as HTMLElement);
    const shown = () => clears().map((b) => getComputedStyle(b).visibility);
    const highlights = () =>
      t.ctx.document.shown.value?.root.children?.flatMap((c) => c.children ?? []).find((c) => c.ref === "chains.main")
        ?.params?.highlights ?? [];
    expect(shown()).toEqual(["hidden", "hidden"]);
    for (const q of ["EAL", "KL"]) {
      for (const input of inputs) {
        await userEvent.clear(input);
        await userEvent.type(input, q);
      }
      await t.settled();
      await act(async () => (t.roots[0].querySelector('[data-seqv-key="save"]') as HTMLElement).click());
      await settle();
    }
    expect((highlights() as unknown[]).length).toBe(2);
    expect(shown()).toEqual(["visible", "visible"]);
    // By the keyboard, in the DOM UI: focused and Enter.
    clears()[1].focus();
    await userEvent.keyboard("{Enter}");
    await settle();
    expect(highlights()).toEqual([]);
    expect(shown()).toEqual(["hidden", "hidden"]);
    // The button hid while it had focus: focus went to its field, never to the page.
    expect(document.activeElement).toBe(inputs[1]);
    // One patch: one Undo brings both back; a click in the React UI clears them.
    await act(async () => {
      await t.ctx.history.undo();
    });
    await settle();
    expect((highlights() as unknown[]).length).toBe(2);
    clears()[0].focus();
    await act(async () => clears()[0].click());
    await settle();
    expect(highlights()).toEqual([]);
    expect(document.activeElement).toBe(inputs[0]);
  });

  it("says `n+` when the search stopped at its cap", async () => {
    const t = await both((ctx) => async (text, signal) => ({ ...(await searchOf(ctx, [])(text, signal)), more: true }));
    const inputs = t.roots.map((r) => r.querySelector(".seqv-search-input") as HTMLInputElement);
    for (const input of inputs) await userEvent.type(input, "EAL");
    await t.settled();
    expect(t.text().map((s) => s?.slice(0, 6))).toEqual(["– / 2+", "– / 2+"]);
  });

  it("keeps the field and the buttons where they are in every state (a fixed status slot)", async () => {
    const t = await both((ctx) => searchOf(ctx, []));
    const inputs = t.roots.map((r) => r.querySelector(".seqv-search-input") as HTMLInputElement);
    const geometry = (r: HTMLElement) => {
      const box = (sel: string) => r.querySelector(sel)?.getBoundingClientRect();
      const status = box(".seqv-search-status");
      const field = box(".seqv-search-field");
      const prev = box('[data-seqv-key="prev"]');
      const save = box('[data-seqv-key="save"]');
      return {
        status: status?.width,
        field: [field?.left, field?.width, field?.height],
        prev: prev?.left,
        save: save?.left,
        shownStatus: r.querySelector(".seqv-search-status")?.textContent,
      };
    };
    const seen: string[] = [];
    const states = async (type: string | null) => {
      for (const input of inputs) {
        await userEvent.clear(input);
        if (type) await userEvent.type(input, type);
      }
    };
    const record = () => {
      const [react, dom] = t.roots.map(geometry);
      expect(react).toEqual(dom);
      seen.push(react.shownStatus ?? "");
      return react;
    };
    const idle = record();
    await states("EAL");
    const searching = record(); // before the text settles
    await t.settled();
    const found = record();
    await states("ZZZ");
    await t.settled();
    const none = record();
    await states("boom");
    await t.settled();
    const failed = record();
    expect(seen).toEqual([
      "",
      "searching…",
      "– / 2",
      "no results",
      "The search service is not reachable at the moment",
    ]);
    // The slot has a width, and the same one with every text; nothing else moves.
    expect(idle.status).toBeGreaterThan(40);
    for (const g of [searching, found, none, failed]) expect(g).toEqual({ ...idle, shownStatus: g.shownStatus });
  });

  it("shrinks with a narrow bar, the status slot kept, and takes a row of its own below 480 px", async () => {
    const t = await both((ctx) => searchOf(ctx, []));
    const field = (r: HTMLElement) => r.querySelector(".seqv-search-field")?.getBoundingClientRect() as DOMRect;
    const slot = (r: HTMLElement) => r.querySelector(".seqv-search-status")?.getBoundingClientRect().width;
    const em = parseFloat(getComputedStyle(t.roots[0]).fontSize);
    const widths: Record<number, number[]> = {};
    const content: Record<number, number> = {};
    for (const w of [900, 600, 400]) {
      for (const r of t.roots) (r.parentElement as HTMLElement).style.width = `${w}px`;
      await settle();
      widths[w] = t.roots.map((r) => field(r).width);
      const bar = t.roots[0].querySelector(".seqv-bar") as HTMLElement;
      content[w] = bar.clientWidth - 2 * Number.parseFloat(getComputedStyle(bar).paddingLeft);
      // Both UIs the same, and the status slot as wide as ever.
      expect(widths[w][0]).toBeCloseTo(widths[w][1], 0);
      expect(slot(t.roots[0])).toBeCloseTo(slot(t.roots[1]) ?? 0, 0);
      expect(slot(t.roots[0])).toBeGreaterThan(40);
    }
    expect(widths[900][0]).toBeCloseTo(20 * em, 0);
    // 40% of the bar's content box (the canvas's edge and its padding inside 600 px), less than 20em.
    expect(widths[600][0]).toBeCloseTo(0.4 * content[600], 0);
    expect(widths[600][0]).toBeLessThan(widths[900][0]);
    // Below 480: its own row, the field taking what the buttons leave.
    for (const r of t.roots) {
      const bar = r.querySelector(".seqv-bar") as HTMLElement;
      const search = r.querySelector(".seqv-search")?.getBoundingClientRect() as DOMRect;
      expect(search.width).toBeGreaterThan(bar.getBoundingClientRect().width * 0.8);
    }
    // Not held to 40% of the bar there: it takes what its buttons leave of the row.
    expect(widths[400][0]).toBeGreaterThan(0.5 * content[400]);
  });

  it("shows Keep at its place only when the results can be kept, in both UIs", async () => {
    const t = await both((ctx) => searchOf(ctx, []));
    const saves = () => t.roots.map((r) => r.querySelector('[data-seqv-key="save"]') as HTMLElement);
    const hidden = () => saves().map((b) => getComputedStyle(b).visibility);
    expect(hidden()).toEqual(["hidden", "hidden"]);
    const inputs = t.roots.map((r) => r.querySelector(".seqv-search-input") as HTMLInputElement);
    for (const input of inputs) await userEvent.type(input, "EAL");
    await t.settled();
    expect(hidden()).toEqual(["visible", "visible"]);
  });
});

describe("the Search part outside the bar (LS2)", () => {
  it("takes its own `search`, in React and in the DOM, with no viewer around it", async () => {
    const { ctx, advance } = await openedContext();
    const asked: string[] = [];
    const search = searchOf(ctx, asked);
    const react = document.createElement("div");
    const dom = document.createElement("div");
    document.body.append(react, dom);
    const root = createRoot(react);
    await act(async () => root.render(<ReactSearch context={ctx} search={search} />));
    const part = DomSearch(dom, ctx, { search });
    cleanups.push(async () => {
      await act(async () => root.unmount());
      part.dispose();
      react.remove();
      dom.remove();
      ctx.dispose();
    });
    const texts = () => [react, dom].map((c) => c.querySelector(".seqv-search")?.textContent ?? null);
    // Both stand alone, the same in structure and text.
    expect(texts()).toEqual(
      Array(2).fill("‹Previous result›Next resultKeepKeep these results highlightedRemove every kept search highlight"),
    );
    for (const c of [react, dom]) await userEvent.type(c.querySelector(".seqv-search-input") as HTMLElement, "EAL");
    await act(async () => advance(SEARCH_SETTLE_MS));
    await settle();
    expect(asked).toEqual(["EAL", "EAL"]);
    expect(texts()).toEqual(
      Array(2).fill(
        "– / 2‹Previous result›Next resultKeepKeep these results highlightedRemove every kept search highlight",
      ),
    );
    await userEvent.type(dom.querySelector(".seqv-search-input") as HTMLElement, "{Enter}");
    await settle();
    expect(ctx.focus.current.value.get("chains.main")).toMatchObject({ start: 23, stop: 26 });
  });

  it("draws nothing without a `search`, from either source", async () => {
    const { ctx } = await openedContext();
    const react = document.createElement("div");
    const dom = document.createElement("div");
    document.body.append(react, dom);
    const root = createRoot(react);
    await act(async () => root.render(<ReactSearch context={ctx} />));
    const part = DomSearch(dom, ctx, {});
    expect([react.firstChild, dom.firstChild]).toEqual([null, null]);
    part.dispose();
    await act(async () => root.unmount());
    react.remove();
    dom.remove();
    ctx.dispose();
  });
});
