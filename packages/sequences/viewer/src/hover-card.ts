/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The hover card (seqspec:interaction#the-hover-card): an interactive card for an item that has actions, beside the
// tooltip, which never takes the pointer. A host names item actions (`UI.ItemActions`: registered `ActionEntry`s whose
// target
// is an `ItemTarget`); resting on such an item for `UI.HoverCardDelayMs` opens the card, the hover's label and a
// button per action that applies. It stays while the pointer goes from the item into it (it opens where the pointer
// is, over the item's edge, and a grace of `HOVER_CARD_GRACE_MS` covers the gap) and while the pointer is in it; it
// closes on leaving both, on Escape or a press outside (the UI parts), on a pan or zoom, and on a new document. Both
// UIs draw it from this one controller; a host may replace the part, restyle it, or leave it out.

import type { ActionEntry } from "./actions.ts";
import { SeqViewerConfig } from "./config.ts";
import type { SeqViewerContext } from "./context.ts";
import { HoverKind } from "./hit.ts";
import type { HoverLabel } from "./labels.ts";
import { type SeqState, state } from "./state.ts";

/** An annotation's item, what an item action acts on. */
export interface ItemTarget {
  /** The section and track it was hovered in; null: not known. */
  readonly section: string | null;
  readonly track: string | null;
  /** The annotation's ref and the item's index in it; `id` its `id`, null: none. */
  readonly annotation: string;
  readonly item: number;
  readonly id: string | null;
  /** The split row it was hovered on (a value's row of an expanded split track); -1: none (the track collapsed). */
  readonly split: number;
}
export interface HoverCardState {
  readonly target: ItemTarget;
  readonly actions: readonly { readonly name: string; readonly title: string }[];
  /** The hover's label when it opened. */
  readonly label: HoverLabel | null;
  /** Where it opened: the section's ref and the pointer there, CSS px on its surface. */
  readonly section: string;
  readonly x: number;
  readonly y: number;
}
export interface HoverCard {
  readonly current: SeqState<HoverCardState | null>;
  /** The pointer came into the card, or left it (the UI parts call these). */
  enter(): void;
  leave(): void;
  /** Runs one of its actions on its target (origin `pointer`), then closes. */
  run(name: string): void;
  close(): void;
  dispose(): void;
}

/**
 * The item the hover holds (an item or a position of one, `owner` holding it) as an action's target, with the item
 * actions `UI.ItemActions` names that apply to it; null: no item, or none applies. Edge rate: the card's opening and a
 * track menu's (seqspec:interaction#the-hover-card).
 */
export function itemActionsOf(
  ctx: SeqViewerContext,
  owner: string,
): { readonly target: ItemTarget; readonly actions: { name: string; title: string }[] } | null {
  const h = ctx.hover.live.current;
  const view = ctx.graph.assembled.value?.view;
  const set = view?.set;
  const isItem = h.owner === owner && (h.kind === HoverKind.item || h.kind === HoverKind.seq);
  const a = isItem ? set?.annotations[h.item.annotation] : undefined;
  const section = ctx.layout.sectionRef(h.section);
  if (!view || !set || !a || section === null) return null;
  const item = h.item.item;
  const key = a.itemKey[item];
  const target: ItemTarget = {
    section,
    track: h.track >= 0 ? set.strings.get(view.resolved.tracks.ref[h.track]) : null,
    annotation: set.strings.get(a.key),
    item,
    id: key >= 0 ? set.strings.get(key) : null,
    split: h.split >= 0 ? h.split : -1,
  };
  const actions: { name: string; title: string }[] = [];
  for (const name of ctx.config.get(SeqViewerConfig.UI.ItemActions)) {
    const entry = ctx.registry.get("actions", name) as ActionEntry<ItemTarget> | undefined;
    if (ctx.actions.can(name, target))
      actions.push({ name, title: entry?.titleFor?.(ctx, target) ?? entry?.title ?? name });
  }
  return actions.length > 0 ? { target, actions } : null;
}

/** The grace that lets the pointer cross from the item into the card, and the delay before it closes. */
export const HOVER_CARD_GRACE_MS = 300;

/**
 * The controller, one per mounted UI. Inert when `UI.ItemActions` is empty (read at mount). The hover handler runs at
 * pointer rate and compares numbers only: moving within one item allocates nothing; a new item sets one timer.
 */
export function createHoverCard(ctx: SeqViewerContext): HoverCard {
  const current = state<HoverCardState | null>(null);
  const names = ctx.config.get(SeqViewerConfig.UI.ItemActions);
  const delay = ctx.config.get(SeqViewerConfig.UI.HoverCardDelayMs);
  const actions = ctx.actions;
  // The item the pointer rests on (annotation id, item; -1: none) and the open card's.
  let annotation = -1;
  let item = -1;
  let openAnnotation = -1;
  let openItem = -1;
  let inCard = false;
  let opening = 0;
  let closing = 0;
  const cancel = (id: number) => id && ctx.scheduler.clearTimer(id);
  /** The card goes; a pending opening (the next item's) stays. */
  const hide = () => {
    cancel(closing);
    closing = 0;
    openAnnotation = openItem = -1;
    inCard = false;
    if (current.value) current.set(null);
  };
  const close = () => {
    cancel(opening);
    opening = 0;
    hide();
  };
  const closeLater = () => {
    if (closing || !current.value) return;
    closing = ctx.scheduler.setTimer(() => {
      closing = 0;
      if (!inCard && !(annotation === openAnnotation && item === openItem)) hide();
    }, HOVER_CARD_GRACE_MS);
  };
  const open = () => {
    opening = 0;
    // The pointer went into the card it crossed to: that card stays.
    if (inCard) return;
    const found = itemActionsOf(ctx, "pointer");
    if (!found || found.target.item !== item) return;
    const p = ctx.pointer.current;
    [openAnnotation, openItem] = [annotation, item];
    const { target, actions: list } = found;
    current.set({
      target,
      actions: list,
      label: ctx.labels.current.value,
      section: target.section as string,
      x: p.x,
      y: p.y,
    });
  };
  const onHover = () => {
    const h = ctx.hover.live.current;
    const isItem =
      h.owner === "pointer" && (h.kind === HoverKind.item || h.kind === HoverKind.seq) && h.item.annotation >= 0;
    const a = isItem ? h.item.annotation : -1;
    const i = isItem ? h.item.item : -1;
    if (a === annotation && i === item) return;
    annotation = a;
    item = i;
    cancel(opening);
    opening = 0;
    if (a === openAnnotation && i === openItem) {
      cancel(closing);
      closing = 0;
      return;
    }
    if (current.value) {
      if (inCard) return;
      closeLater();
    }
    // A new item gets its own delay even while another's card is closing, so moving from one item to the next opens
    // the next one's card (the old card closes after its grace, the new one opens after the delay).
    if (isItem) opening = ctx.scheduler.setTimer(open, delay);
  };
  const moves: { unsubscribe(): void }[] = [];
  // A pan or zoom (any change of a viewport, a resize included), or a new document, closes it; a reassembly of the same
  // document (a skin, an arrival) does not (seqspec:interaction#replays-and-carries). The viewports' replay guard only
  // matters when the subscription is made: `sections` emits on a new document, which closes the card anyway.
  let doc = ctx.document.shown.value;
  const subs =
    names.length === 0
      ? []
      : [
          ctx.hover.live.version.subscribe(onHover),
          ctx.document.shown.subscribe((d) => {
            if (d === doc) return;
            doc = d;
            close();
          }),
          ctx.viewports.sections.subscribe((refs) => {
            for (const s of moves) s.unsubscribe();
            moves.length = 0;
            for (const ref of refs) {
              const vp = ctx.viewports.get(ref);
              if (!vp) continue;
              let replay = true;
              moves.push(vp.live.version.subscribe(() => !replay && current.value && close()));
              replay = false;
            }
          }),
        ];
  return {
    current,
    enter() {
      inCard = true;
      cancel(closing);
      cancel(opening);
      closing = opening = 0;
    },
    leave() {
      inCard = false;
      closeLater();
    },
    run(name) {
      const c = current.value;
      if (!c) return;
      close();
      actions.run(name, c.target, "pointer");
    },
    close,
    dispose() {
      close();
      for (const s of [...subs, ...moves]) s.unsubscribe();
    },
  };
}
