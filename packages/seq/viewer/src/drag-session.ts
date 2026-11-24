/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A drag ends wherever the button is released (_spec/README.md "Input"): pointer capture fails silently where it is
// refused (an iframe boundary, a pointer already up), so a drag a mount starts also listens on the element's own
// window, for as long as it runs.

/** A drag's temporary listeners: on the window of `element` (an iframe's own) and its document, only while a drag runs. */
export interface DragSession {
  /** True between `start` and the end. */
  readonly active: boolean;
  /** Begins listening (a drag still running is ended first, once). Call on a press the viewer took. */
  start(): void;
  /** Ends the drag: the listeners are removed and `onEnd` runs, once per `start`. No-op when none runs. */
  end(): void;
  /** Takes the listeners away without calling `onEnd` (an unmount). */
  dispose(): void;
}

/**
 * The session for drags on `element`: `onEnd` runs once when the button goes up or the pointer is cancelled anywhere
 * in the window, a `pointermove` arrives with no button held, the window loses its focus or the document turns hidden,
 * or `end` is called (the element's own `pointerup`, `pointercancel` and `lostpointercapture`).
 */
export function createDragSession(element: Element, onEnd: () => void): DragSession {
  const doc = element.ownerDocument;
  const win = doc.defaultView;
  let active = false;
  const ender = () => end();
  const onMove = (e: Event) => {
    if ((e as PointerEvent).buttons === 0) end();
  };
  const onVisibility = () => {
    if (doc.visibilityState === 'hidden') end();
  };
  const detach = () => {
    win?.removeEventListener('pointerup', ender);
    win?.removeEventListener('pointercancel', ender);
    win?.removeEventListener('pointermove', onMove);
    win?.removeEventListener('blur', ender);
    doc.removeEventListener('visibilitychange', onVisibility);
  };
  const end = () => {
    if (!active) return;
    active = false;
    detach();
    onEnd();
  };
  return {
    get active() {
      return active;
    },
    start() {
      end();
      active = true;
      win?.addEventListener('pointerup', ender);
      win?.addEventListener('pointercancel', ender);
      win?.addEventListener('pointermove', onMove);
      win?.addEventListener('blur', ender);
      doc.addEventListener('visibilitychange', onVisibility);
    },
    end,
    dispose() {
      active = false;
      detach();
    },
  };
}
