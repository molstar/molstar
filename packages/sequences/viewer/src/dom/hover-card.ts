/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The hover card's DOM part and placement (seqspec:interaction#the-hover-card), the React part's too: the card at the
// point it opened, right
// of and below the pointer by a few pixels (over the item's edge, so the pointer crosses into it), flipped where it
// would leave the viewer. Escape and a press outside close it; its buttons are ordinary buttons.

import { createHoverCard, type HoverCard, type HoverCardState } from "../hover-card.ts";
import type { SeqViewerContext } from "../index.ts";
import { type DomPart, h, off, type PartProps } from "./core.ts";

/** Places the card at its state's point on its section's surface, inside the viewer (its offset parent). */
export function placeCard(surfaces: ReadonlyMap<string, HTMLElement>, el: HTMLElement, s: HoverCardState): void {
  const surface = surfaces.get(s.section);
  const parent = el.offsetParent as HTMLElement | null;
  if (!surface || !parent) return;
  let ox = 0;
  let oy = 0;
  for (let e: Element | null = surface; e instanceof HTMLElement && e !== parent; e = e.offsetParent) {
    ox += e.offsetLeft;
    oy += e.offsetTop;
  }
  let x = ox + s.x + 4;
  let y = oy + s.y + 4;
  if (x + el.offsetWidth > parent.clientWidth) x = ox + s.x - 4 - el.offsetWidth;
  if (y + el.offsetHeight > parent.clientHeight) y = oy + s.y - 4 - el.offsetHeight;
  el.style.transform = `translate(${Math.max(0, Math.round(x))}px, ${Math.max(0, Math.round(y))}px)`;
}

/** Escape anywhere and a press outside the card close it while it is open; returns the unlisten. */
export function dismissOn(card: HoverCard, el: () => HTMLElement | null): () => void {
  const doc = globalThis.document;
  const key = (e: KeyboardEvent) => e.key === "Escape" && card.current.value && card.close();
  const press = (e: PointerEvent) => {
    const box = el();
    if (card.current.value && box && !box.contains(e.target as Node)) card.close();
  };
  doc.addEventListener("keydown", key, true);
  doc.addEventListener("pointerdown", press, true);
  return () => {
    doc.removeEventListener("keydown", key, true);
    doc.removeEventListener("pointerdown", press, true);
  };
}

/** The card's content: the hover's label, then a button per action. */
function content(card: HoverCard, s: HoverCardState): HTMLElement[] {
  const label = s.label;
  return [
    ...(label?.title ? [h("div", { class: "seqv-tooltip-title" }, label.title)] : []),
    ...(label?.lines ?? []).map((l) =>
      h("div", { class: "seqv-tooltip-line" }, h("span", { class: "seqv-quiet" }, l.key), ` ${l.value}`),
    ),
    h(
      "div",
      { class: "seqv-hover-card-actions" },
      ...s.actions.map((a) => {
        const b = h("button", { type: "button", class: "seqv-hover-card-action" }, a.title);
        b.addEventListener("click", () => card.run(a.name));
        return b;
      }),
    ),
  ];
}

export const HoverCardPart: DomPart<PartProps> = (host: HTMLElement, ctx: SeqViewerContext, { ui }) => {
  // A canvas overlay, placed at its UI's surfaces: nothing without a UI (as React's without a provider).
  if (!ui) return { dispose() {} };
  const card = createHoverCard(ctx);
  let el: HTMLElement | null = null;
  const sub = card.current.subscribe((s) => {
    el?.remove();
    el = null;
    if (!s) return;
    const box = h("div", { class: "seqv-hover-card", role: "dialog", "aria-label": s.label?.title ?? "Item" });
    box.addEventListener("pointerenter", () => card.enter());
    box.addEventListener("pointerleave", () => card.leave());
    box.append(...content(card, s));
    host.append(box);
    placeCard(ui.surfaces, box, s);
    el = box;
  });
  const unlisten = dismissOn(card, () => el);
  return {
    dispose: () => {
      unlisten();
      off([sub, card]);
      el?.remove();
    },
  };
};
