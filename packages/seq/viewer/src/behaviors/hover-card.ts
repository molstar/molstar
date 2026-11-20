/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type ReadonlyState, State } from '@molstar/core/reactive';
import type { GridLay, RowsLay, SectionLayout, ViewLayout } from '@molstar/seqx-core/render';
import { DomHost } from '@molstar/ui/dom/host';
import { DefaultAnnouncerOptions } from '@molstar/seqx-viewer/announcer';
import { describerFor, type Place } from '@molstar/seqx-viewer/describe';
import type { Cursor } from '@molstar/seqx-viewer/models';
import { twoAxesLayer } from './two-axes.js';
import {
  type CardBox,
  CardEstimate,
  type HoverCard,
  type HoverCardPlace,
  type HoverCardPort,
  type HoverCardPortFactory,
  placeCard,
} from '@molstar/seqx-viewer/hover-card';
import type { BehaviorContext, BehaviorProvider } from '@molstar/seqx-viewer/viewer-spec';

// The hover card (_spec/README.md "Behaviors", hoverCard; "Ports", Hover card): what the pointer rests on, in words,
// after the hover rests (the announcer's rest), one card. Its content is the announcer's description (`describe.ts`)
// and its place `placeCard`'s, both headless; the port draws it. Hidden while a drag pans or zooms (the input tells it of
// the pointer and its drag), on leave, on Esc, a press and a wheel; shown at the keyboard cursor when a key moves it;
// a dragged rectangle on a grid shows its summary (`gridSummary`). Per pointer move it stores two numbers: the card is
// made on rest, never per move, and one timer waits for the rest however many moves come (a move pushes the rest on).
// No card outlives what it describes: a frame laid out again (a zoom, a pan, an edit, a load, data arriving) places the
// cursor's card again and sends a hover's away until the hover rests again. The port is a param, given by the default
// spec (`hoverCardOverlay()`): this module draws nothing itself.

export interface HoverCardParams {
  /** The port's adapter (the default spec gives `@molstar/ui`'s card over the canvas; none: nothing is drawn). */
  readonly port?: HoverCardPortFactory;
  /** How long the hover rests before the card shows, ms (the announcer's, 500). */
  readonly restMs?: number;
}

/** What is shown: the card and its place. */
export interface ShownHoverCard {
  readonly card: HoverCard;
  readonly place: HoverCardPlace;
}

export interface HoverCardBehavior {
  /** What the card shows now; undefined while it is hidden. */
  readonly shown: ReadonlyState<ShownHoverCard | undefined>;
  /**
   * The pointer is at `at` (x, y in the canvas), a drag running or not (the input, each move): copied, no allocation
   * (in an array: fractional numbers passed as arguments are boxed).
   */
  pointer(at: Float64Array, dragging: boolean): void;
  /** The pointer left the canvas. */
  leave(): void;
  /** A press, a wheel or Esc: the card goes until the hover changes and rests again. */
  interrupt(): void;
  /** A drag ended. */
  released(): void;
  /** A key moved the keyboard cursor: the card shows at its cell after the next frame (the one that revealed it). */
  cursorMoved(): void;
  /** A frame was drawn (the viewer, each frame), laid out again or not: no allocation. */
  framed(relaid: boolean): void;
  /** The mounted canvas's element: the port's adapter is made over it. */
  attach(container: HTMLElement): () => void;
  dispose(): void;
}

/** How long a rectangle's summary waits after a change while it is dragged, ms: the summary once a frame or so. */
const SummaryMs = 50;

/** A box `r` px around a point. */
const around = (x: number, y: number, rx: number, ry: number): CardBox => ({
  x: x - rx,
  y: y - ry,
  width: 2 * rx,
  height: 2 * ry,
});

/**
 * The box of the hovered mark in the canvas: in a tracks layout, the item's or the cell's columns over its band (a rows
 * band's row); elsewhere (a grid, a ring, lines of text) a box around the pointer as large as a cell on each side.
 * Read on rest from the layout; never a frame's work.
 */
function anchorOf(layout: ViewLayout, place: Place, pointer: { readonly x: number; readonly y: number }): CardBox {
  const section = layout.sections.find((s) => s.section.address === place.section);
  const px = pointer.x;
  const py = pointer.y;
  if (!section) return around(px, py, 6, 6);
  const scale = section.axis ? Math.abs(section.axis.scale) : 0;
  const cell = Math.max(6, Math.min(64, scale));
  const track = section.tracks.find((t) => t.layers.some((l) => l.key === place.layer));
  if (section.layout.name !== 'tracks' || !track || !section.axis) return around(px, py, cell, cell);
  const layer = track.layers.find((l) => l.key === place.layer)!;
  let x0 = Number.NaN;
  let x1 = Number.NaN;
  let y = track.box.y;
  let height = track.box.height;
  const item = place.element && layer.type.items?.(layer.lay).find((d) => d.key === place.element!.key);
  if (item) {
    x0 = item.x0;
    x1 = item.x1;
    if (item.y !== undefined && item.height !== undefined) {
      y = item.y;
      height = item.height;
    }
  } else if (!Number.isNaN(place.position)) {
    const p = Math.floor(place.position);
    const a = section.axis.x(p);
    const b = section.axis.x(p + 1);
    x0 = Math.min(a, b);
    x1 = Math.max(a, b);
  }
  // a rows band is tall: the row under the pointer, not the band
  if (layer.type.name === 'rows' && !item) {
    const r = Math.max(1, (layer.lay as RowsLay).rowScale);
    y = Math.max(track.box.y, py - r);
    height = Math.min(track.box.y + track.box.height, py + r) - y;
  }
  const box = { x: x0, y, width: x1 - x0, height };
  // what was found holds the pointer, else the layout is not the one hit (a strip's frame): the pointer's box
  const holds =
    Number.isFinite(x0) &&
    Number.isFinite(x1) &&
    px >= box.x - 1 &&
    px <= box.x + box.width + 1 &&
    py >= box.y - 1 &&
    py <= box.y + box.height + 1;
  return holds ? box : around(px, py, cell, cell);
}

/**
 * The box of the keyboard cursor's cell: the cell's columns over the section's tracks; on a section of two axes its
 * row's cell (a rows cell, a grid's), or its row's header (the label column to the data).
 */
function cellBox(section: SectionLayout, c: Cursor): CardBox | undefined {
  const l = c.ySystem !== undefined && c.yPosition !== undefined ? twoAxesLayer(section, c.ySystem) : undefined;
  if (l?.type.name === 'rows' && c.yPosition !== undefined) {
    const lay = l.lay as RowsLay;
    // where the cursor is drawn: its row (a hidden one's nearest shown), its cell or its header (the label column)
    const i = lay.cursorRow(c.yPosition);
    if (i < 0) return undefined;
    const y = lay.rowY(i);
    const x0 = c.part === 'header' ? lay.labelX : lay.axis.x(c.position);
    const x1 = c.part === 'header' ? lay.numberX : lay.axis.x(c.position + 1);
    if (!Number.isFinite(x0) || !Number.isFinite(x1)) return undefined;
    return { x: Math.min(x0, x1), y, width: Math.abs(x1 - x0), height: lay.rowScale };
  }
  if (l?.type.name === 'grid' && c.yPosition !== undefined) {
    const lay = l.lay as GridLay;
    const ax = lay.axis;
    const ay = lay.axisY;
    if (!ax || !ay || !ax.cell(c.position)) return undefined;
    const x = Math.min(ax.sx0, ax.sx1);
    const width = Math.abs(ax.sx1 - ax.sx0);
    if (!ay.cell(c.yPosition)) return undefined;
    return { x, y: Math.min(ay.sx0, ay.sx1), width, height: Math.abs(ay.sx1 - ay.sx0) };
  }
  const t = section.tracks[0];
  if (!section.axis || !t) return undefined;
  const a = section.axis.x(c.position);
  const b = section.axis.x(c.position + 1);
  if (!Number.isFinite(a) || !Number.isFinite(b)) return undefined;
  const top = section.tracks[0].box.y;
  const last = section.tracks[section.tracks.length - 1].box;
  return { x: Math.min(a, b), y: top, width: Math.abs(b - a), height: last.y + last.height - top };
}

/**
 * What a hover hits, as a card names it: its section, layer, target, element, block and state, and its cell (the
 * position's, unless an element), on one axis or two. `update` keeps the place's and says whether it was the same.
 * Strings by reference and numbers in an array: no allocation a move.
 */
class HitIdentity {
  private section: string | undefined;
  private layer: string | undefined;
  private target: string | undefined;
  private ySystem: string | undefined;
  private state: string | undefined;
  private owner: string | undefined;
  private key: string | number | undefined;
  private members: unknown;
  /** The model the hit was read from: a load or an edit keeping the section, layer and key is another hit. */
  private data: unknown;
  /** The cell along x and y (−1: none), the block and the y block (−1: none). */
  private readonly cells = new Float64Array([-1, -1, -1, -1]);

  update(p: Place, data: unknown): boolean {
    const e = p.element;
    const x = e || Number.isNaN(p.position) ? -1 : Math.floor(p.position);
    const y = p.yPosition === undefined || Number.isNaN(p.yPosition) ? -1 : Math.floor(p.yPosition);
    const c = this.cells;
    const same =
      this.section === p.section &&
      this.layer === p.layer &&
      this.target === p.target &&
      this.ySystem === p.ySystem &&
      this.state === p.state &&
      this.members === p.members &&
      this.data === data &&
      this.owner === e?.owner &&
      this.key === e?.key &&
      c[0] === x &&
      c[1] === y &&
      c[2] === (p.block ?? -1) &&
      c[3] === (p.yBlock ?? -1);
    this.section = p.section;
    this.layer = p.layer;
    this.target = p.target;
    this.ySystem = p.ySystem;
    this.state = p.state;
    this.members = p.members;
    this.data = data;
    this.owner = e?.owner;
    this.key = e?.key;
    c[0] = x;
    c[1] = y;
    c[2] = p.block ?? -1;
    c[3] = p.yBlock ?? -1;
    // nothing hovered is never the same hover (a leave and a return show the card again)
    return same && p.section !== undefined;
  }
}

/** What the card shows, by what asked for it. */
type Showing = 'none' | 'hover' | 'cursor' | 'rect';

function hoverCardBehavior(context: BehaviorContext, params: HoverCardParams = {}): HoverCardBehavior {
  const { hover, selection, config } = context;
  const restMs = params.restMs ?? DefaultAnnouncerOptions.restMs;
  // the viewer's one describer: the live region says the same words from it
  const words = describerFor(context);
  const shown = new State<ShownHoverCard | undefined>(undefined);
  // the port, made once (Ports, Hover card): at the mount, or headless at the first card; `made` says it was asked
  let port: HoverCardPort | undefined;
  let made = false;
  // the pointer: NaN while it is off the canvas (kept in an array: a number written to a closure's variable is boxed)
  const at = new Float64Array([Number.NaN, Number.NaN]);
  let dragging = false;
  let showing: Showing = 'none';
  /** A key moved the keyboard cursor: its card is made after the next frame (the one that revealed it). */
  let cursorDue = false;
  let changedAt = 0;
  let timer: unknown;
  let summaryTimer: unknown;
  let disposed = false;
  const sameHit = new HitIdentity();

  const portOf = (): HoverCardPort | undefined => {
    if (!made) {
      made = true;
      port = params.port?.(undefined);
    }
    return port;
  };

  /** Lets the port go: disposed when the viewer owns it (`disposeWithViewer` not false). */
  const release = () => {
    const p = port;
    port = undefined;
    made = false;
    if (p && p.disposeWithViewer !== false) p.dispose?.();
  };

  const show = (card: HoverCard, anchor: CardBox, pointer: { x: number; y: number }, by: Showing) => {
    const layout = context.layout();
    if (!layout || !config.hover.card) return;
    const view = { x: 0, y: 0, width: layout.width, height: layout.height };
    const placed = placeCard(anchor, pointer, view, CardEstimate);
    const place: HoverCardPlace = { anchor, pointer, view, ...placed };
    showing = by;
    shown.set({ card, place });
    portOf()?.show(card, place);
  };

  const hide = () => {
    showing = 'none';
    if (shown.value === undefined) return;
    shown.set(undefined);
    portOf()?.hide();
  };

  const stopTimer = () => {
    DomHost.clearTimeout(timer);
    timer = undefined;
  };

  /** The hover's rest begins (again): the card shows `restMs` after it, unless something moves first. */
  const arm = () => {
    changedAt = DomHost.now();
    if (timer === undefined && hover.active) timer = DomHost.setTimeout(rest, restMs);
  };

  const rest = () => {
    timer = undefined;
    if (disposed) return;
    const wait = changedAt + restMs - DomHost.now();
    // a move since pushed the rest on: wait the rest of it
    if (wait > 0) {
      timer = DomHost.setTimeout(rest, wait);
      return;
    }
    if (dragging || Number.isNaN(at[0]) || !hover.active || !config.hover.card) return;
    const layout = context.layout();
    const card = words.describe(hover.current);
    if (!layout || !card) return;
    const pointer = { x: at[0], y: at[1] };
    show(card, anchorOf(layout, hover.current, pointer), pointer, 'hover');
  };

  /** The card at the keyboard cursor's cell, from the last frame. */
  const showCursor = () => {
    const c = selection.cursor.value;
    const layout = context.layout();
    const section = c && layout?.sections.find((s) => s.section.address === c.section);
    if (!c || !section) return hide();
    const card = words.describe({ cursor: c });
    const box = card && cellBox(section, c);
    if (!card || !box) return hide();
    show(card, box, { x: box.x + box.width / 2, y: box.y + box.height / 2 }, 'cursor');
  };

  const stopHover = hover.version.listen(() => {
    // a re-hit of what was hovered (a zoom or a pan under a resting pointer: the same element, cell or target) is no
    // change: the card stays (a cursor's is placed again with the frame); the same place of another model is a change
    if (sameHit.update(hover.current, context.spec.model.value)) return;
    // a rectangle's summary stays while its drag moves the hover
    if (!dragging) hide();
    cursorDue = false;
    arm();
  });

  // a rectangle dragged on a grid: its summary, a little after each change while the drag runs
  const summary = () => {
    summaryTimer = undefined;
    if (disposed || !dragging || !config.hover.card) return;
    const card = words.describe({ rect: selection.value.value.location });
    if (!card) return;
    const pointer = { x: at[0], y: at[1] };
    show(card, around(at[0], at[1], 6, 6), pointer, 'rect');
  };
  const stopSelection = selection.value.listen((v) => {
    if (!dragging || !v.location.rects.length || summaryTimer !== undefined) return;
    summaryTimer = DomHost.setTimeout(summary, SummaryMs);
  });

  return {
    shown: shown.asReadonly,
    pointer(p, drag) {
      at[0] = p[0];
      at[1] = p[1];
      if (drag && !dragging) {
        dragging = true;
        stopTimer();
        hide();
      } else if (!drag && dragging) dragging = false;
    },
    leave() {
      at[0] = at[1] = Number.NaN;
      stopTimer();
      hide();
    },
    interrupt() {
      stopTimer();
      cursorDue = false;
      hide();
    },
    released() {
      if (!dragging) return;
      dragging = false;
      hide();
    },
    cursorMoved() {
      cursorDue = true;
    },
    framed(relaid) {
      if (cursorDue) {
        cursorDue = false;
        showCursor();
        return;
      }
      if (!relaid || showing === 'none' || showing === 'rect') return;
      // the view or the data moved under the card: the cursor's is placed again; a hover's goes until it rests again
      if (showing === 'cursor') showCursor();
      else {
        hide();
        arm();
      }
    },
    attach(container) {
      if (!config.hover.card) return () => {};
      // made once: a card shown before the mount (headless) gives its port up first
      if (made) release();
      made = true;
      port = params.port?.(container);
      return () => {
        hide();
        release();
      };
    },
    dispose() {
      disposed = true;
      stopHover();
      stopSelection();
      stopTimer();
      DomHost.clearTimeout(summaryTimer);
      release();
    },
  };
}

/**
 * The hover card (`hoverCard`, in the default spec with `hoverCardOverlay()` as its port; `hover.card` turns it off).
 * Without a `port` nothing is drawn (the card's state still moves): a host gives its adapter, and bundles no default.
 */
export const HoverCardProvider: BehaviorProvider<HoverCardParams | undefined, HoverCardBehavior> = {
  name: 'hoverCard',
  create: (context, params) => hoverCardBehavior(context, params),
};
