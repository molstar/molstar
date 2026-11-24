/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type ReadonlyState, State } from '@molstar/core/reactive';
import { DomHost } from '@molstar/ui/dom/host';
import { locationText } from '@molstar/seqx-core/render';
import { type Describer, describerFor } from '@molstar/seqx-viewer/describe';
import type { SeqViewer } from './viewer.js';

// What the viewer says to a reader who does not see the canvas (_spec/README.md "Input", Announcements): one polite
// live region's text, made here DOM-free and shown by the mount. Three things are said: the keyboard cursor's place
// (its section, system, position, letter and the items there) as it moves, the hover's (when it rests: a pointer
// crossing the picture says nothing) and a selection's change (its size and extent, once). Everything goes through one
// rate limit: a message never follows another within `minGapMs`, the latest one waiting its turn (a held arrow key
// says where it is each `minGapMs`, not each key). Why one text a time: a live region queues what it is given, and a
// reader hears the queue after the pointer has moved on. A message waiting its turn is replaced by a later one of the same
// or a higher priority only (a selection over the cursor over the hover: a hover never replaces a selection's final
// message), and the pointer leaving drops a waiting hover. While a gesture is dragged the cursor and the hover say
// nothing, and the selection says its end once the drag ends. Each message is a token (a counter and its text), so
// equal words said again are a new message: the mount clears the region and sets it again.

export interface AnnouncerOptions {
  /** How long the hover rests before it is said, ms. */
  readonly restMs?: number;
  /** The least time between two messages, ms. */
  readonly minGapMs?: number;
}

/** A message said: its counter makes equal texts apart. */
export interface Said {
  readonly id: number;
  readonly text: string;
}

export type AnnouncePriority = 'hover' | 'notice' | 'cursor' | 'selection';
/** A part's notice (`viewer.notify`) over the hover: a hover never replaces it, nor does the pointer leaving drop it. */
const Rank: Readonly<Record<AnnouncePriority, number>> = { hover: 1, notice: 2, cursor: 3, selection: 4 };

const CursorKeys = /^cursor\.(prev|next|up|down|page-up|page-down|first-row|last-row|extend-.+|add-.+)$/;

export const DefaultAnnouncerOptions = { restMs: 500, minGapMs: 700 } as const;

export class Announcer {
  private _count = 0;
  private readonly _message = new State<Said>({ id: 0, text: '' });
  /** What the live region says now: each message a new token, so equal words are said again. */
  readonly message: ReadonlyState<Said> = this._message.asReadonly;
  private readonly _viewer: SeqViewer;
  private readonly _rest: number;
  private readonly _gap: number;
  private _last = Number.NEGATIVE_INFINITY;
  private _pending: { text: string; priority: AnnouncePriority } | undefined;
  private _dragTimer: unknown;
  private _selectionDue = false;
  private _gapTimer: unknown;
  private _restTimer: unknown;
  private readonly _stop: (() => void)[] = [];
  private _selected = '';
  /** The keyboard cursor when the last command ran: a command that moved it says its message. */
  private _cursor: unknown;
  private _disposed = false;

  constructor(viewer: SeqViewer, options: AnnouncerOptions = {}) {
    this._viewer = viewer;
    this._rest = options.restMs ?? DefaultAnnouncerOptions.restMs;
    this._gap = options.minGapMs ?? DefaultAnnouncerOptions.minGapMs;
    // the viewer's one describer: the hover card says the same words from it
    this._words = describerFor({
      spec: viewer.spec,
      store: viewer.store,
      addresses: viewer.addresses,
      view: viewer.view,
      layout: () => viewer.frame?.layout,
    });
    this._cursor = viewer.selection.cursor.value;
    this._stop.push(
      // a cursor key is said whether or not the cursor moved (at the end of its axis it stays: the key is still heard);
      // the pointer's own cursor is shown, not said
      viewer.commandRan.listen(({ id, message }) => {
        const cursor = viewer.selection.cursor.value;
        const moved = cursor !== this._cursor;
        this._cursor = cursor;
        if (this.dragging) return;
        if (CursorKeys.test(id)) this.say(this.cursorText(), 'cursor');
        // another command that moved the keyboard cursor says its own words (`rondo:partner`: where, and why there)
        else if (moved && cursor && message) this.say(message, 'cursor');
      }),
      // a part's state (a ring's links loading), once: above the hover, below the cursor and the selection
      viewer.notice.listen(({ text }) => this.say(text, 'notice')),
      viewer.hover.version.listen(() => this.hovered()),
      viewer.selection.value.listen(() => this.selectionChanged()),
    );
  }

  /** The text of a say, at once when the rate limit allows, else when it does (the latest of those waiting). */
  say(text: string, priority: AnnouncePriority = 'cursor'): void {
    if (this._disposed || !text) return;
    const wait = this._last + this._gap - DomHost.now();
    if (wait <= 0) {
      this._pending = undefined;
      this._emit(text);
      return;
    }
    // a waiting message of a higher priority is kept
    if (this._pending && Rank[this._pending.priority] > Rank[priority]) return;
    this._pending = { text, priority };
    if (this._gapTimer === undefined)
      this._gapTimer = DomHost.setTimeout(() => {
        this._gapTimer = undefined;
        const t = this._pending;
        this._pending = undefined;
        if (t && !this._disposed) this._emit(t.text);
      }, wait);
  }

  private _emit(text: string): void {
    this._last = DomHost.now();
    this._message.set({ id: ++this._count, text });
  }

  // ---- the hover, when it rests ----

  private hovered(): void {
    DomHost.clearTimeout(this._restTimer);
    this._restTimer = undefined;
    if (!this._viewer.hover.active) {
      // the pointer left: what waited to say its hover is dropped
      if (this._pending?.priority === 'hover') {
        this._pending = undefined;
        DomHost.clearTimeout(this._gapTimer);
        this._gapTimer = undefined;
      }
      return;
    }
    if (this.dragging) return;
    this._restTimer = DomHost.setTimeout(() => {
      this._restTimer = undefined;
      if (!this.dragging) this.say(this.hoverText(), 'hover');
    }, this._rest);
  }

  /** Whether a gesture is being dragged (a drag-select, a press on the canvas). */
  private get dragging(): boolean {
    const sel = this._viewer.selection;
    return sel.dragging || sel.rectDragging || this._viewer.input.dragging;
  }

  // ---- the selection, once a change ----

  private selectionChanged(): void {
    // during a drag the selection changes with each move: it says its end once, when the drag ends
    if (this.dragging) {
      this._selectionDue = true;
      if (this._dragTimer === undefined) this.waitForDragEnd();
      return;
    }
    this.sayOfSelection();
  }

  private waitForDragEnd(): void {
    this._dragTimer = DomHost.setTimeout(() => {
      this._dragTimer = undefined;
      if (this._disposed) return;
      if (this.dragging) return this.waitForDragEnd();
      if (this._selectionDue) this.sayOfSelection();
      this._selectionDue = false;
    }, 100);
  }

  private sayOfSelection(): void {
    const shown = this._viewer.selection.value.value.shown;
    const text = shown.parts.length
      ? `Selected ${locationText(this._viewer.spec.model.value, shown)}`
      : this._selected
        ? 'Selection cleared'
        : '';
    if (text === this._selected && shown.parts.length) return;
    this._selected = shown.parts.length ? text : '';
    this.say(text, 'selection');
  }

  // ---- the words: the hover card's (`describe.ts`), one description for both ----

  private readonly _words: Describer;

  /** `protein, p53 12 (kabat 10A), letter L, items: Domain A, p53 3–8; Domain B, p53 5–12`. */
  placeText(section: string, system: string, position: number): string {
    return this._words.text({ section, system, position });
  }

  /** The cursor's place: its cell, with its row on a section of two axes, or its row's header. */
  cursorText(): string {
    const c = this._viewer.selection.cursor.value;
    return c ? this._words.text({ cursor: c }) : '';
  }

  hoverText(): string {
    return this._viewer.hover.active ? this._words.text(this._viewer.hover.current) : '';
  }

  dispose(): void {
    this._disposed = true;
    for (const stop of this._stop) stop();
    this._stop.length = 0;
    DomHost.clearTimeout(this._gapTimer);
    DomHost.clearTimeout(this._restTimer);
    DomHost.clearTimeout(this._dragTimer);
  }
}
