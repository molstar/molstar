/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// `ViewSurface` (seqspec:interaction#the-default-ui) as a `DomPart`: the current view's sections in a column, each a
// frame with its header (its row holds the hidden tracks' chips) and a mounted surface, plus the chrome placed over the
// canvas: the row controls and the track header's menu. Those follow the viewport by styles set in subscriptions; their
// nodes change only when the rows in view or the hovered header change. As src/react/section.tsx.

import { Height } from "@molstar/seq-core/format";
import { type SeqViewerContext, shownSections, type Unsubscribable } from "../index.ts";
import {
  type DomHandle,
  type DomPart,
  type DomSlots,
  type DomUI,
  derived,
  h,
  later,
  lose,
  off,
  type RowProps,
  setAttr,
  watch,
} from "./core.ts";
import { focusIndicator, placeScroll, setStyle, trackBox } from "./shared.ts";

interface Frame {
  readonly root: HTMLElement;
  dispose(): void;
}

/** The canvas: it needs its UI (`createDomUI`), whose surfaces map and composition it fills and reads. */
export const ViewSurface: DomPart<{ readonly ui: DomUI }> = (host, ctx, { ui }) => {
  const root = h("div", { class: "seqv-view" });
  host.append(root);
  const frames = new Map<string, Frame>();
  const subs = derived(
    [ctx.graph.assembled],
    () => shownSections(ctx),
    (refs) => {
      for (const [ref, f] of frames) {
        if (refs.includes(ref)) continue;
        // Focus inside a section that leaves lands on the viewer's root (seqspec:interaction#accessibility).
        lose(f.root, root, () => f.dispose());
        frames.delete(ref);
      }
      refs.forEach((ref, i) => {
        const at = () => root.children[i] ?? null;
        const f = frames.get(ref) ?? sectionFrame(ctx, ui, ref, (el) => root.insertBefore(el, at()));
        frames.set(ref, f);
        if (at() !== f.root) root.insertBefore(f.root, at());
      });
    },
  );
  return {
    dispose() {
      off([...subs, ...frames.values()]);
      root.remove();
    },
  };
};

/** A section's frame; `place` puts it in the document before its surface mounts (the engine measures its host). */
function sectionFrame(ctx: SeqViewerContext, ui: DomUI, section: string, place: (el: HTMLElement) => void): Frame {
  // The section's chrome, as the layout's composition places it (seqspec:interaction#slots-and-chrome).
  const placed = <K extends "SectionHeader" | "RowControl" | "TrackActions">(k: K) =>
    ui.placed.has(k) ? ui.slots[k] : null;
  const [SectionHeader, RowControl, TrackActions] = [
    placed("SectionHeader"),
    placed("RowControl"),
    placed("TrackActions"),
  ];
  const root = h("section", { class: "seqv-section" });
  const parts: { dispose(): void }[] = [];
  if (SectionHeader) parts.push(SectionHeader(root, ctx, { ui, section }));
  const el = h("div", { class: "seqv-surface" });
  const stage = h("div", { class: "seqv-stage" }, el);
  root.append(stage);
  place(root);
  const s = ui.surface(ctx, section, el);
  ui.surfaces.set(section, s.element);
  const subs = [
    ctx.layout.contentHeight(section)?.subscribe((height) => {
      setStyle(el, "height", `${Math.max(40, Math.ceil(height))}px`);
      setAttr(el, "data-height", ctx.layout.get(section)?.current.heightMode === Height.box ? "box" : "content");
    }),
  ];
  const indicator = focusIndicator(el, s.element);
  if (RowControl) parts.push(rowControls(stage, ctx, ui, section, RowControl));
  if (TrackActions) parts.push(trackActionsAt(stage, ctx, ui, section, TrackActions));
  return {
    root,
    dispose() {
      off([...parts, ...subs]);
      indicator();
      if (ui.surfaces.get(section) === s.element) ui.surfaces.delete(section);
      s.dispose();
      root.remove();
    },
  };
}

/**
 * The `TrackActions` slot at the hovered or cursor's header's right edge, for the track `trackBox` decides: kept while
 * the pointer is over it, while its menu is open (an `aria-expanded` control in it, so a host's menu counts) and, with
 * nothing hovered, while focus is inside it; the keyboard cursor's blur clears the hover when focus moves into the menu.
 */
function trackActionsAt(
  stage: HTMLElement,
  ctx: SeqViewerContext,
  ui: DomUI,
  section: string,
  Actions: DomSlots["TrackActions"],
): { dispose(): void } {
  let box: HTMLElement | null = null;
  let part: { dispose(): void } | null = null;
  let follow: (Unsubscribable | undefined)[] = [];
  const clear = () => {
    const [b, p, f] = [box, part, follow];
    [box, part, follow] = [null, null, []];
    const drop = () => {
      off([p, ...f]);
      b?.remove();
    };
    // Focus inside the box (its trigger, its menu) lands on the section's surface.
    if (b) lose(b, stage, drop);
    else drop();
  };
  const tracker = trackBox(ctx, section, {
    box: () => box,
    show(track) {
      clear();
      if (!track) return;
      const b = h("div", { class: "seqv-track-actions" });
      b.addEventListener("pointerenter", () => tracker.hold(track));
      b.addEventListener("pointerleave", () => tracker.hold(null));
      b.addEventListener("focusout", tracker.later);
      stage.append(b);
      box = b;
      part = Actions(b, ctx, { ui, section, track, row: () => tracker.row() });
      // The box is new: placed once now, then as the viewport or the layout moves.
      tracker.place();
      const [vp, layout] = [ctx.viewports.get(section), ctx.layout.get(section)];
      follow = [
        vp ? later(vp.live.version, tracker.place) : undefined,
        layout ? later(layout.version, () => tracker.update()) : undefined,
      ];
    },
  });
  const subs = watch([ctx.hover.live.version, ctx.ui.keyCursor], () => tracker.update());
  return {
    dispose() {
      off(subs);
      tracker.dispose();
      clear();
    },
  };
}

interface Control {
  readonly el: HTMLElement;
  readonly dataRow: number;
  readonly part: DomHandle<RowProps>;
}
/** The `RowControl` slot beside each visible row of the section's rows layer, in the header's last pixels. */
function rowControls(
  stage: HTMLElement,
  ctx: SeqViewerContext,
  ui: DomUI,
  section: string,
  Control: DomSlots["RowControl"],
): { dispose(): void } {
  const rows = ctx.layout.visibleRows(section);
  const layout = ctx.layout.get(section);
  const vp = ctx.viewports.get(section);
  if (!rows || !layout) return { dispose() {} };
  const inner = h("div");
  const box = h("div", { class: "seqv-row-controls" }, inner);
  stage.append(box);
  const shown = new Map<number, Control>();
  const drop = (d: number, c: Control) => [c.part.dispose(), c.el.remove(), shown.delete(d)];
  // Edge-rate: when the rows in view change, not per scroll. Layout is read here, in a subscription, never in a render.
  const update = () => {
    const [r, l, w] = [rows.current, layout.current, ctx.layout.props(section).controlWidth];
    setStyle(box, "top", `${l.scrolledTop}px`);
    setStyle(box, "left", `${l.headerWidth - w}px`);
    setStyle(box, "width", `${w}px`);
    const want = new Map<number, number>();
    for (let i = 0; i < r.count; i++) want.set(r.displayRow[i], r.dataRow[i]);
    for (const [d, c] of shown) if (want.get(d) !== c.dataRow) drop(d, c);
    for (let i = 0; i < r.count; i++) {
      const [displayRow, dataRow] = [r.displayRow[i], r.dataRow[i]];
      let c = shown.get(displayRow);
      if (!c) {
        const el = h("div", { class: "seqv-row-control" });
        inner.insertBefore(el, inner.children[i] ?? null);
        c = { el, dataRow, part: Control(el, ctx, { ui, section, displayRow, dataRow }) };
        shown.set(displayRow, c);
      } else if (inner.children[i] !== c.el) inner.insertBefore(c.el, inner.children[i] ?? null);
      setStyle(c.el, "top", `${r.top[i]}px`);
      setStyle(c.el, "height", `${r.height[i]}px`);
    }
  };
  // The scroll follows by a custom property (skin.css), written when the quarter pixel changes: no string per frame.
  const at = { quarter: Number.NaN, page: Number.NaN };
  const subs = [
    rows.version.subscribe(update),
    vp?.live.version.subscribe(() => placeScroll(inner, vp.live.current.scrollY, at)),
  ];
  return {
    dispose() {
      off(subs);
      for (const [d, c] of shown) drop(d, c);
      box.remove();
    },
  };
}
