/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The hover card's port (_spec/README.md "Ports", Hover card): the viewer computes the card as data (a title, rows of a
// label and a value, the announcer's text) and its place (the anchor it stays clear of, the pointer, the view, a side
// from `placeCard`), headless; the port draws it. The default adapter is `@molstar/ui`'s surface over the canvas
// (`mount-hover-card.ts`); a host passes its own, drawn with its components and styles; a test records the calls.

/** One row of the card: what it names and its value (`kabat`, `82B`). */
export interface HoverCardRow {
  readonly label: string;
  readonly value: string;
  /** The row names an item at the place (a feature under it, a group's member): what `hover.cardItems` counts. */
  readonly item?: boolean;
}

/** What the card says: its title, its rows, and the sentence the live region says of the same place. */
export interface HoverCard {
  readonly title: string;
  /** Its rows: past `hover.cardItems` item rows, the first of them and a row `+k more` (`more` holds the rest). */
  readonly rows: readonly HoverCardRow[];
  /** The announcer's words for the place: the card and the live region read the same (every item). */
  readonly text: string;
  /** The item rows the card leaves to its `+k more` row: a host's card may list them. */
  readonly more?: readonly HoverCardRow[];
  /**
   * Matches past the cap not made into rows (a tooltip table's rows over a bin, counted by its index): the `+k more` row
   * says `k = more.length + moreCount`.
   */
  readonly moreCount?: number;
}

/**
 * The card a port draws: at most `maxItems` item rows, the rest counted in a row `+k more` and kept in `more` (with the
 * matches never made, `moreCount`); the card itself when it has no more.
 */
export function cutItems(card: HoverCard, maxItems: number): HoverCard {
  let items = 0;
  for (const r of card.rows) if (r.item) items++;
  const uncounted = card.moreCount ?? 0;
  if (items <= maxItems && !uncounted) return card;
  const rows: HoverCardRow[] = [];
  const more: HoverCardRow[] = [];
  let n = 0;
  for (const r of card.rows) {
    if (!r.item) rows.push(r);
    else if (n++ < maxItems) rows.push(r);
    else more.push(r);
  }
  rows.push({ label: `+${(more.length + uncounted).toLocaleString('en-US')} more`, value: '' });
  return { ...card, rows, more };
}

/** A box in the canvas's CSS px. */
export interface CardBox {
  readonly x: number;
  readonly y: number;
  readonly width: number;
  readonly height: number;
}

/** Where a card goes: above or below the anchor, else beside it (left or right) when neither has room. */
export type CardSide = 'above' | 'below' | 'right' | 'left';

/** The card's place as boxes: what it stays clear of, the pointer, the view, and the side `placeCard` suggests. */
export interface HoverCardPlace {
  /** The hovered cell or item (the keyboard cursor's cell): the card is never over it. */
  readonly anchor: CardBox;
  /** The pointer (the anchor's middle for the keyboard cursor). */
  readonly pointer: { readonly x: number; readonly y: number };
  /** The view the card stays inside: the canvas. */
  readonly view: CardBox;
  /** The side and position `placeCard` gives a card of `CardEstimate`'s size; an adapter places its own by its size. */
  readonly side: CardSide;
  readonly x: number;
  readonly y: number;
}

/** What the viewer asks of the port. */
export interface HoverCardPort {
  /** Shows the card at its place (one card: a show replaces the last). */
  show(card: HoverCard, place: HoverCardPlace): void;
  hide(): void;
  /** False: the host owns it (one adapter shared by several viewers), and the viewer leaves it. */
  readonly disposeWithViewer?: boolean;
  dispose?(): void;
}

/**
 * The port's factory: given the element the canvas is mounted in (none headless). Called once a viewer: when it mounts,
 * or headless at the first card (a mount after a headless card disposes that one first, when the viewer owns it).
 */
export type HoverCardPortFactory = (container: HTMLElement | undefined) => HoverCardPort;

/** The size `HoverCardPlace` suggests a side for: a title and a few rows. */
export const CardEstimate = { width: 200, height: 96 } as const;

/** The gap between the card and the anchor, and the pointer, CSS px. */
const CardGap = 8;

/**
 * The card's position (its top-left) and side: beside the pointer (its left at the pointer, or its right where the
 * view ends), above the anchor, else below it, else to its right or left; kept inside the view. Never over the anchor
 * while the view has room on some side; with none, above it, clamped.
 */
export function placeCard(
  anchor: CardBox,
  pointer: { readonly x: number; readonly y: number },
  view: CardBox,
  size: { readonly width: number; readonly height: number },
  gap = CardGap,
): { side: CardSide; x: number; y: number } {
  const clamp = (v: number, lo: number, hi: number) => Math.max(lo, Math.min(hi, v));
  const right = view.x + view.width;
  const bottom = view.y + view.height;
  // across: from the pointer on, else ending at it; inside the view
  const across = () => {
    const x = pointer.x + gap + size.width <= right ? pointer.x + gap : pointer.x - gap - size.width;
    return clamp(x, view.x, right - size.width);
  };
  const above = anchor.y - gap - size.height;
  if (above >= view.y) return { side: 'above', x: across(), y: above };
  const below = anchor.y + anchor.height + gap;
  if (below + size.height <= bottom) return { side: 'below', x: across(), y: below };
  // neither above nor below: beside the anchor, level with the pointer
  const y = clamp(pointer.y - size.height / 2, view.y, bottom - size.height);
  const toRight = anchor.x + anchor.width + gap;
  if (toRight + size.width <= right) return { side: 'right', x: toRight, y };
  const toLeft = anchor.x - gap - size.width;
  if (toLeft >= view.x) return { side: 'left', x: toLeft, y };
  return { side: 'above', x: across(), y: clamp(above, view.y, bottom - size.height) };
}
