/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Gestures that lose their release (prism AGENTS.md, gestures attach once; the scroll area's drag bug): a gesture's
// tests press its target, move, lose the release (nothing comes; or only the capture is lost), N times, and check that
// the listeners on every target (the element, the document, the window) are as many as before, and that a whole
// gesture after still works. jsdom has no pointer capture: it is stubbed, and a `PointerEvent` made if it has none.

/** Net listeners on every target, counted through a patched `addEventListener`/`removeEventListener`. */
export function countListeners(): { readonly count: number; stop(): void } {
  const proto = EventTarget.prototype;
  const add = proto.addEventListener;
  const remove = proto.removeEventListener;
  const live = new Map<EventTarget, Set<string>>();
  const ids = new WeakMap<object, number>();
  let next = 1;
  const idOf = (fn: unknown) => {
    if (typeof fn !== 'function' && (typeof fn !== 'object' || fn === null)) return 0;
    let id = ids.get(fn as object);
    if (!id) ids.set(fn as object, (id = next++));
    return id;
  };
  const keyOf = (type: string, fn: unknown, o: unknown) =>
    `${type}|${idOf(fn)}|${typeof o === 'boolean' ? o : !!(o as { capture?: boolean } | undefined)?.capture}`;
  let count = 0;
  proto.addEventListener = function (this: EventTarget, type: string, fn: unknown, o?: unknown) {
    if (fn) {
      let keys = live.get(this);
      if (!keys) live.set(this, (keys = new Set()));
      const k = keyOf(type, fn, o);
      if (!keys.has(k)) {
        keys.add(k);
        count++;
      }
    }
    return add.call(this, type, fn as EventListener, o as AddEventListenerOptions);
  } as typeof proto.addEventListener;
  proto.removeEventListener = function (this: EventTarget, type: string, fn: unknown, o?: unknown) {
    if (live.get(this)?.delete(keyOf(type, fn, o))) count--;
    return remove.call(this, type, fn as EventListener, o as EventListenerOptions);
  } as typeof proto.removeEventListener;
  return {
    get count() {
      return count;
    },
    stop() {
      proto.addEventListener = add;
      proto.removeEventListener = remove;
    },
  };
}

/** jsdom: pointer capture stubbed (kept per element), and `PointerEvent` when it has none. */
export function stubPointerCapture(): void {
  const proto = HTMLElement.prototype as HTMLElement & { __captured?: Set<number> };
  if (!('PointerEvent' in globalThis)) {
    (globalThis as { PointerEvent?: unknown }).PointerEvent = class extends MouseEvent {
      readonly pointerId: number;
      readonly pointerType: string;
      constructor(type: string, init: PointerEventInit = {}) {
        super(type, init);
        this.pointerId = init.pointerId ?? 1;
        this.pointerType = init.pointerType ?? 'mouse';
      }
    };
  }
  const held = new WeakMap<Element, Set<number>>();
  proto.setPointerCapture = function (this: Element, id: number) {
    let s = held.get(this);
    if (!s) held.set(this, (s = new Set()));
    s.add(id);
  };
  proto.releasePointerCapture = function (this: Element, id: number) {
    held.get(this)?.delete(id);
  };
  proto.hasPointerCapture = function (this: Element, id: number) {
    return !!held.get(this)?.has(id);
  };
}

export interface PointerAt {
  x: number;
  y: number;
}

/** A pointer event at (x, y), the main button held unless `buttons` says otherwise. */
export function pointer(type: string, at: PointerAt, buttons = 1): PointerEvent {
  return new PointerEvent(type, {
    bubbles: true,
    cancelable: true,
    composed: true,
    clientX: at.x,
    clientY: at.y,
    button: 0,
    buttons: type === 'pointerup' ? 0 : buttons,
    pointerId: 1,
    pointerType: 'mouse',
  } as PointerEventInit);
}

export interface GestureCase {
  /** The element pressed (found again each time: a gesture may re-render). */
  target(): Element;
  /** Where the press is, and where the move goes. */
  from: PointerAt;
  to: PointerAt;
  /** After a whole gesture (press, move, release): whether it did what it does (a value moved, a size changed). */
  worked(): boolean;
  /** Resets what `worked` reads, before the whole gesture. */
  reset?(): void;
}

/**
 * Presses, moves and loses the release `n` times (each other time nothing comes; between, only the capture is lost),
 * then makes a whole gesture. Returns the listeners added meanwhile (should be 0) and whether the gesture worked.
 */
export function lostReleases(c: GestureCase, n = 20): { grown: number; worked: boolean } {
  const counter = countListeners();
  try {
    const before = counter.count;
    for (let i = 0; i < n; i++) {
      const el = c.target();
      el.dispatchEvent(pointer('pointerdown', c.from));
      el.dispatchEvent(pointer('pointermove', c.to));
      // The release lost: nothing comes (a window left, a native dialog), or only the capture goes.
      if (i % 2 === 1) el.dispatchEvent(pointer('lostpointercapture', c.to, 0));
    }
    const grown = counter.count - before;
    c.reset?.();
    const el = c.target();
    el.dispatchEvent(pointer('pointerdown', c.from));
    el.dispatchEvent(pointer('pointermove', c.to));
    el.dispatchEvent(pointer('pointerup', c.to));
    return { grown, worked: c.worked() };
  } finally {
    counter.stop();
  }
}
