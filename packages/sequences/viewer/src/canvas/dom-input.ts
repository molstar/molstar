/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// DOM input → `ctx.input` (seqspec:viewer#the-render-backend; seqspec:interaction#input-events): pointer, wheel, key,
// leave and focus events filled into reused, DOM-free records with x and y in the surface's CSS px. The primary pointer
// only (a second finger is none);
// a press subscribes the document's moves and release (in the capture phase) until its release, so the drag follows the
// pointer anywhere on the page, read in the pressing surface's px; the surface's own moves are its hover while nothing
// is pressed. A press never starts the browser's drag and drop where the surface is the target. Nothing here decides
// what a gesture does: behaviors do, from the base events. Handlers are bound once; an event allocates nothing.

import { type InputSink, Inputs, Modifier, PointerType } from '../input.js';

const LINE_PX = 16;

const modifiers = (e: MouseEvent | KeyboardEvent) =>
  (e.shiftKey ? Modifier.shift : 0) |
  (e.ctrlKey ? Modifier.ctrl : 0) |
  (e.altKey ? Modifier.alt : 0) |
  (e.metaKey ? Modifier.meta : 0);

/** The target is where the page types (an input, a textarea, a select, an editable element): its keys are not ours. */
function typing(t: EventTarget | null): boolean {
  if (!(t instanceof HTMLElement)) return false;
  if (t.isContentEditable) return true;
  const n = t.tagName;
  return n === 'INPUT' || n === 'TEXTAREA' || n === 'SELECT';
}

/**
 * Listens on `el` (the surface) and the document; `section` is read per event (a view change moves the section's
 * index). `stop` removes every listener; with `notify`, it also ends a press this surface started (a cancel) and tells
 * the dispatcher the surface lost focus if it had it (removing a focused element fires no blur).
 */
export class DomInput {
  private readonly el: HTMLElement;
  private readonly sink: InputSink;
  private readonly section: () => number;
  private readonly p = Inputs.pointer();
  private readonly w = Inputs.wheel();
  private readonly k = Inputs.key();
  private pressed = false;
  private focused = false;
  /**
   * The element's rect during a press, so a move anywhere reads in its px: read at the press and again when anything
   * scrolls (a page scroll mid-drag), never per move.
   */
  private rect: DOMRect | null = null;

  constructor(el: HTMLElement, sink: InputSink, section: () => number) {
    this.el = el;
    this.sink = sink;
    this.section = section;
    el.addEventListener('pointerdown', this.down);
    el.addEventListener('pointermove', this.move);
    el.addEventListener('pointerleave', this.leave);
    el.addEventListener('wheel', this.wheel, { passive: false });
    el.addEventListener('keydown', this.key);
    el.addEventListener('focus', this.focus);
    el.addEventListener('blur', this.blur);
    el.addEventListener('dragstart', this.noDrag);
  }

  stop(notify: boolean): void {
    const el = this.el;
    el.removeEventListener('pointerdown', this.down);
    el.removeEventListener('pointermove', this.move);
    el.removeEventListener('pointerleave', this.leave);
    el.removeEventListener('wheel', this.wheel);
    el.removeEventListener('keydown', this.key);
    el.removeEventListener('focus', this.focus);
    el.removeEventListener('blur', this.blur);
    el.removeEventListener('dragstart', this.noDrag);
    if (notify && this.pressed) this.sink.cancel(this.section());
    if (notify && this.focused) this.sink.focus(this.section(), false);
    this.release();
    this.focused = false;
  }

  /** The press over: the document's listeners removed. */
  private release() {
    if (this.pressed) {
      document.removeEventListener('pointermove', this.docMove, true);
      document.removeEventListener('pointerup', this.up, true);
      document.removeEventListener('pointercancel', this.up, true);
      document.removeEventListener('scroll', this.scrolled, true);
    }
    this.pressed = false;
    this.rect = null;
  }

  private fill(e: PointerEvent, type: number): void {
    const p = this.p;
    p.type = type;
    p.section = this.section();
    const r = this.rect;
    if (r) {
      // Outside the surface during a press: from the rect read at the press, as the target may be the document.
      p.x = e.clientX - r.left;
      p.y = e.clientY - r.top;
    } else {
      p.x = e.offsetX;
      p.y = e.offsetY;
    }
    p.buttons = e.buttons;
    p.button = type === PointerType.move ? -1 : e.button;
    p.modifiers = modifiers(e);
    // A press's and a release's time; a move's is not read (`timeStamp` is a fractional double: reading it boxes one per
    // move).
    if (type !== PointerType.move) p.time = e.timeStamp;
    p.pointerId = e.pointerId;
    this.sink.pointer(p);
  }
  private readonly down = (e: PointerEvent) => {
    if (!e.isPrimary) return;
    if (!this.pressed) {
      document.addEventListener('pointermove', this.docMove, true);
      document.addEventListener('pointerup', this.up, true);
      document.addEventListener('pointercancel', this.up, true);
      document.addEventListener('scroll', this.scrolled, true);
    }
    this.pressed = true;
    this.rect = this.el.getBoundingClientRect();
    this.el.focus({ preventScroll: true });
    this.fill(e, PointerType.down);
  };
  /**
   * The surface's own moves: its hover while nothing is pressed here. A move with a button down is another surface's
   * press (followed there, in its own px).
   */
  private readonly move = (e: PointerEvent) => {
    if (e.isPrimary && !this.pressed && e.buttons === 0) this.fill(e, PointerType.move);
  };
  /** A press never starts the browser's drag and drop (which cancels the pointer), where the surface is its target. */
  private readonly noDrag = (e: Event) => e.preventDefault();
  /**
   * The document's moves during a press, anywhere on the page. One with no button down follows a release the page never
   * delivered: passed on (the dispatcher ends the press there), and the press is over here too.
   */
  private readonly docMove = (e: PointerEvent) => {
    if (!this.pressed || !e.isPrimary) return;
    this.fill(e, PointerType.move);
    if (e.buttons === 0) this.release();
  };
  private readonly up = (e: PointerEvent) => {
    if (!this.pressed || !e.isPrimary) return;
    this.fill(e, e.type === 'pointercancel' ? PointerType.cancel : PointerType.up);
    this.release();
  };
  /** Something scrolled during a press: the rect again, so the drag's moves stay in the surface's px. */
  private readonly scrolled = () => {
    if (this.pressed) this.rect = this.el.getBoundingClientRect();
  };
  private readonly leave = () => {
    if (!this.pressed) this.sink.leave(this.section());
  };
  private readonly wheel = (e: WheelEvent) => {
    const scale = e.deltaMode === 1 ? LINE_PX : e.deltaMode === 2 ? this.el.clientHeight : 1;
    const w = this.w;
    w.section = this.section();
    w.x = e.offsetX;
    w.y = e.offsetY;
    w.deltaX = e.deltaX * scale;
    w.deltaY = e.deltaY * scale;
    w.modifiers = modifiers(e);
    w.time = e.timeStamp;
    if (this.sink.wheel(w)) e.preventDefault();
  };
  private readonly key = (e: KeyboardEvent) => {
    const k = this.k;
    k.section = this.section();
    k.key = e.key;
    k.code = e.code;
    k.modifiers = modifiers(e);
    k.typing = typing(e.target);
    k.time = e.timeStamp;
    if (this.sink.key(k)) e.preventDefault();
  };
  private readonly focus = () => {
    this.focused = true;
    this.sink.focus(this.section(), true);
  };
  private readonly blur = () => {
    this.focused = false;
    this.sink.focus(this.section(), false);
  };
}
