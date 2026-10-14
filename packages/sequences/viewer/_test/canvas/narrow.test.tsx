/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The viewer's own parts on a phone: in a frame 308 px wide (a story on a 390 px screen), the view tabs, the legend's
// value ramp and the Search box never overflow the viewer's root, in both default UIs on real surfaces
// (seqspec:interaction#the-default-ui, seqspec:interaction#the-legend, seqspec:interaction#slots-and-chrome).

import { createSeqViewSpec } from "@molstar/seq-core/builder";
import { act } from "react";
import { createRoot } from "react-dom/client";
import { afterEach, describe, expect, it } from "vitest";
import { windowScheduler } from "../../src/canvas/index.ts";
import * as D from "../../src/dom/index.ts";
import "../../src/skin.css";
import { createSeqViewer, DefaultSeqViewerSpec, type SearchFn, SeqViewerConfig } from "../../src/index.ts";
import { Search as ReactSearch, SeqViewerLayout, SeqViewerProvider } from "../../src/react/index.ts";

(globalThis as { IS_REACT_ACT_ENVIRONMENT?: boolean }).IS_REACT_ACT_ENVIRONMENT = true;

const cleanups: (() => void | Promise<void>)[] = [];
afterEach(async () => {
  for (const c of cleanups.splice(0).reverse()) await c();
});

const A = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ";
/** Several views with one-word names (a tab does not wrap its text further), and a heat track (a value ramp, captioned). */
function doc() {
  const b = createSeqViewSpec({ title: "Narrow" });
  const a = b.data.sequence({ ref: "a", label: "A", value: A });
  const d = b.data.annotation({
    ref: "d",
    geometry: "values",
    target: a,
    value: { index: [0, 5, 9], value: [0, 5.5, 20.25] },
  });
  for (const name of ["Overviewchain", "Domainslist", "Variantspos", "Contactsall"])
    b.view({ ref: name.toLowerCase(), label: name }).block({ sequence: a }).track({ label: "Heat" }).color({ data: d });
  return b.getDocument();
}

/** A viewer 308 px wide (a story's frame on a 390 px screen; the UI named) with the legend expanded and a search; resolved once the surfaces are sized. */
async function open(ui: "react" | "dom") {
  const spec = DefaultSeqViewerSpec();
  const ctx = createSeqViewer(
    {
      ...spec,
      config: [...(spec.config ?? []), [SeqViewerConfig.UI.Legend, "expanded"], [SeqViewerConfig.UI.Chrome, "full"]],
    },
    { scheduler: windowScheduler() },
  );
  const r = await ctx.setDocument(doc());
  if (!r.ok) throw new Error(JSON.stringify(r.issues));
  const container = document.createElement("div");
  container.style.cssText = "width:308px;--seq-font-ui:Arial;--seq-font-mono:'Courier New'";
  document.body.appendChild(container);
  cleanups.push(
    () => ctx.dispose(),
    () => container.remove(),
  );
  const search: SearchFn = async () => ({ results: [] });
  if (ui === "react") {
    const root = createRoot(container);
    await act(async () =>
      root.render(
        <SeqViewerProvider context={ctx} ui={{ search }}>
          <SeqViewerLayout />
        </SeqViewerProvider>,
      ),
    );
    cleanups.push(() => act(async () => root.unmount()));
  } else {
    const m = D.mountSeqViewer(container, ctx, { search });
    cleanups.push(() => m.dispose());
  }
  for (let k = 0; k < 30; k++) await new Promise((res) => requestAnimationFrame(res));
  return { ctx, container, root: container.querySelector(".seqv-root") as HTMLElement };
}

/** What sticks out of `root` on the right: the visible elements' classes. */
function overflow(root: HTMLElement): string[] {
  const edge = root.getBoundingClientRect().right + 0.5;
  const out: string[] = [];
  for (const el of root.querySelectorAll<HTMLElement>("*")) {
    // Not what is hidden until a hover (the hint's bubble, the tooltip), nor the screen-reader text.
    if (el.closest(".seqv-sr-only, [popover]") || getComputedStyle(el).visibility === "hidden") continue;
    const r = el.getBoundingClientRect();
    if (r.width > 0 && r.right > edge) out.push(`${el.className || el.tagName} right ${r.right.toFixed(1)} > ${edge}`);
  }
  return out;
}

/** The fixed status slot's left and right edges lie inside the field's. */
function statusInside(search: HTMLElement): boolean[] {
  const field = (search.querySelector(".seqv-search-field") as HTMLElement).getBoundingClientRect();
  const status = (search.querySelector(".seqv-search-status") as HTMLElement).getBoundingClientRect();
  return [status.left >= field.left, status.right <= field.right + 0.5];
}

describe.each(["react", "dom"] as const)("at 390 px (%s)", (ui) => {
  it("the view tabs stay inside the bar", async () => {
    const { root } = await open(ui);
    expect(root.querySelectorAll(".seqv-tab").length).toBe(4);
    expect(overflow(root.querySelector(".seqv-bar") as HTMLElement)).toEqual([]);
  });

  it("the legend's value ramp stays inside the legend, its bar shrinking before its labels do", async () => {
    const { root } = await open(ui);
    const ramp = root.querySelector(".seqv-ramp") as HTMLElement;
    expect(ramp).toBeTruthy();
    expect(overflow(root.querySelector(".seqv-legend") as HTMLElement)).toEqual([]);
    expect((ramp.querySelector(".seqv-ramp-bar") as HTMLElement).getBoundingClientRect().width).toBeGreaterThanOrEqual(
      40,
    );
  });

  it("the Search box in the bar takes the bar's width, its status slot still inside the field", async () => {
    const { root } = await open(ui);
    const search = root.querySelector(".seqv-search") as HTMLElement;
    expect(search).toBeTruthy();
    expect(overflow(root.querySelector(".seqv-bar") as HTMLElement)).toEqual([]);
    expect(statusInside(search)).toEqual([true, true]);
  });

  it("the Search part placed on its own shrinks to the room it has", async () => {
    const { ctx } = await open(ui);
    const host = document.createElement("div");
    // A flex row, as a host's toolbar is (a flex item does not shrink below its content unless it may).
    host.style.cssText = "display:flex;width:308px;--seq-font-ui:Arial";
    document.body.appendChild(host);
    cleanups.push(() => host.remove());
    const search: SearchFn = async () => ({ results: [] });
    if (ui === "react") {
      const r = createRoot(host);
      await act(async () => r.render(<ReactSearch context={ctx} search={search} />));
      cleanups.push(() => act(async () => r.unmount()));
    } else {
      const m = D.Search(host, ctx, { search });
      cleanups.push(() => m.dispose());
    }
    const box = host.querySelector(".seqv-search") as HTMLElement;
    expect(box).toBeTruthy();
    expect(overflow(host)).toEqual([]);
    expect(statusInside(box)).toEqual([true, true]);
  });

  it("the whole viewer does not overflow", async () => {
    const { root } = await open(ui);
    expect(overflow(root)).toEqual([]);
  });
});
