/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// DOM input into IN2's dispatcher (CS2; seqspec:viewer#the-render-backend, seqspec:interaction#input-events), in
// Chromium: pointer, wheel, key, leave and focus events become reused records on `ctx.input` (the same objects every
// call), the primary pointer only, a press followed on the document until its release (its listeners there for the
// press only, the rect read again on a scroll, a move with no button ending it; drag and drop never starts where the
// surface is its target; another surface's press is never this one's), the page's default prevented
// only when a subscriber handled it, and a stopped listener (or a disposed surface) hearing nothing more; then through
// a mounted surface, where a drag pans the viewport, focus shows the keyboard's cursor and `ui.cursor` sets the
// element's cursor.

import { createSeqViewSpec } from "@molstar/seq-core/builder";
import { OrderedSet } from "@molstar/seq-core/model";
import { afterEach, describe, expect, it } from "vitest";
import { userEvent } from "vitest/browser";
import { DomInput } from "../../src/canvas/dom-input.js";
import {
  type InputSink,
  type KeyInput,
  Modifier,
  type PointerInput,
  PointerType,
  type WheelInput,
} from "../../src/index.js";
import { runX } from "../../src/runs.js";
import { syntheticFamily } from "../families.js";
import { close, type Opened, opened, settle } from "./open.js";

/** A sink that records what it was given (copies, since the records are reused) and the records' identities. */
function recorder(handled = false) {
  const log: string[] = [];
  const seen = { pointer: new Set<PointerInput>(), wheel: new Set<WheelInput>(), key: new Set<KeyInput>() };
  const sink: InputSink = {
    pointer(e) {
      seen.pointer.add(e);
      log.push(
        `pointer ${e.type} ${e.section} ${e.x},${e.y} b${e.button}/${e.buttons} m${e.modifiers} #${e.pointerId}`,
      );
    },
    wheel(e) {
      seen.wheel.add(e);
      log.push(`wheel ${e.section} ${e.x},${e.y} ${e.deltaX},${e.deltaY} m${e.modifiers}`);
      return handled;
    },
    key(e) {
      seen.key.add(e);
      log.push(`key ${e.key} ${e.code} m${e.modifiers} typing:${e.typing}`);
      return handled;
    },
    focus: (s, f) => log.push(`focus ${s} ${f}`),
    leave: (s) => log.push(`leave ${s}`),
    cancel: (s) => log.push(`cancel ${s}`),
  };
  return { sink, log, seen };
}

const hosts: HTMLElement[] = [];
afterEach(() => {
  for (const h of hosts.splice(0)) h.remove();
});
function element() {
  const el = document.createElement("div");
  el.tabIndex = 0;
  el.style.cssText = "position:absolute;left:100px;top:50px;width:300px;height:200px";
  document.body.appendChild(el);
  hosts.push(el);
  return el;
}
const ptr = (type: string, x: number, y: number, o: PointerEventInit = {}) =>
  new PointerEvent(type, {
    clientX: x,
    clientY: y,
    isPrimary: true,
    pointerId: 1,
    bubbles: true,
    cancelable: true,
    ...o,
  });

describe("DomInput", () => {
  it("fills the pointer record, in the surface's px, and follows a press onto the document until its release", () => {
    const el = element();
    const r = recorder();
    const input = new DomInput(el, r.sink, () => 3);
    el.dispatchEvent(ptr("pointermove", 130, 90));
    el.dispatchEvent(ptr("pointerdown", 130, 90, { buttons: 1, button: 0, shiftKey: true }));
    // Off the surface during the press: from the rect read at the press.
    document.body.dispatchEvent(ptr("pointermove", 500, 400, { buttons: 1 }));
    document.body.dispatchEvent(ptr("pointerup", 500, 400, { button: 0 }));
    // After the release the document is not followed.
    document.body.dispatchEvent(ptr("pointermove", 520, 410));
    el.dispatchEvent(ptr("pointerleave", 0, 0));
    expect(r.log).toEqual([
      "pointer 1 3 30,40 b-1/0 m0 #1",
      "focus 3 true",
      "pointer 0 3 30,40 b0/1 m1 #1",
      "pointer 1 3 400,350 b-1/1 m0 #1",
      "pointer 2 3 400,350 b0/0 m0 #1",
      "leave 3",
    ]);
    // One record for every pointer call.
    expect(r.seen.pointer.size).toBe(1);
    input.stop(true);
  });

  it("forwards the primary pointer only; a second finger is none", () => {
    const el = element();
    const r = recorder();
    const input = new DomInput(el, r.sink, () => 0);
    el.dispatchEvent(ptr("pointerdown", 120, 80, { isPrimary: false, pointerId: 2, buttons: 1 }));
    el.dispatchEvent(ptr("pointermove", 125, 80, { isPrimary: false, pointerId: 2, buttons: 1 }));
    document.dispatchEvent(ptr("pointerup", 125, 80, { isPrimary: false, pointerId: 2 }));
    expect(r.log).toEqual([]);
    el.dispatchEvent(ptr("pointerdown", 120, 80, { buttons: 1, button: 0 }));
    // The other finger's release does not end the primary's press.
    document.dispatchEvent(ptr("pointerup", 125, 80, { isPrimary: false, pointerId: 2 }));
    el.dispatchEvent(ptr("pointerleave", 0, 0));
    expect(r.log.map((l) => l.split(" ")[0])).toEqual(["focus", "pointer"]);
    input.stop(false);
  });

  it("a pointercancel ends a press as a cancel; a leave during a press waits for the release", () => {
    const el = element();
    const r = recorder();
    const input = new DomInput(el, r.sink, () => 0);
    el.dispatchEvent(ptr("pointerdown", 120, 80, { buttons: 1, button: 0 }));
    el.dispatchEvent(ptr("pointerleave", 0, 0));
    document.dispatchEvent(ptr("pointercancel", 120, 80));
    el.dispatchEvent(ptr("pointerleave", 0, 0));
    expect(r.log.map((l) => l.split(" ").slice(0, 2).join(" "))).toEqual([
      "focus 0",
      "pointer 0",
      "pointer 3",
      "leave 0",
    ]);
    input.stop(false);
  });

  it("a press subscribes the document's moves and release until its release; a page scroll re-reads the rect", () => {
    const el = element();
    const r = recorder();
    const input = new DomInput(el, r.sink, () => 0);
    const added: string[] = [];
    const removed: string[] = [];
    const [add, remove] = [document.addEventListener, document.removeEventListener];
    document.addEventListener = function (this: Document, type: string, ...rest: never[]) {
      added.push(type);
      return add.call(this, type, ...(rest as [never]));
    } as typeof document.addEventListener;
    document.removeEventListener = function (this: Document, type: string, ...rest: never[]) {
      removed.push(type);
      return remove.call(this, type, ...(rest as [never]));
    } as typeof document.removeEventListener;
    try {
      // Not pressed: the document's moves are none of this surface's.
      document.body.dispatchEvent(ptr("pointermove", 500, 400));
      el.dispatchEvent(ptr("pointerdown", 120, 80, { buttons: 1, button: 0 }));
      expect(added.sort()).toEqual(["pointercancel", "pointermove", "pointerup", "scroll"]);
      document.body.dispatchEvent(ptr("pointermove", 500, 400, { buttons: 1 }));
      // The page scrolled the surface 30 px up: the next move reads the new rect.
      el.style.top = "20px";
      document.dispatchEvent(new Event("scroll"));
      document.body.dispatchEvent(ptr("pointermove", 500, 400, { buttons: 1 }));
      document.body.dispatchEvent(ptr("pointerup", 500, 400, { button: 0 }));
      expect(removed.sort()).toEqual(["pointercancel", "pointermove", "pointerup", "scroll"]);
      document.body.dispatchEvent(ptr("pointermove", 520, 410, { buttons: 0 }));
      expect(r.log.filter((l) => l.startsWith("pointer"))).toEqual([
        "pointer 0 0 20,30 b0/1 m0 #1",
        "pointer 1 0 400,350 b-1/1 m0 #1",
        "pointer 1 0 400,380 b-1/1 m0 #1",
        "pointer 2 0 400,380 b0/0 m0 #1",
      ]);
    } finally {
      document.addEventListener = add;
      document.removeEventListener = remove;
    }
    input.stop(false);
  });

  it("a move with no button during a press (a release the page lost) is passed on, and the press is over here", () => {
    const el = element();
    const r = recorder();
    const input = new DomInput(el, r.sink, () => 0);
    el.dispatchEvent(ptr("pointerdown", 120, 80, { buttons: 1, button: 0 }));
    document.body.dispatchEvent(ptr("pointermove", 500, 400, { buttons: 0 }));
    // Over: the document's moves are no longer this surface's; a leave is a leave again.
    document.body.dispatchEvent(ptr("pointermove", 510, 400, { buttons: 0 }));
    el.dispatchEvent(ptr("pointerleave", 0, 0));
    expect(r.log.filter((l) => !l.startsWith("focus"))).toEqual([
      "pointer 0 0 20,30 b0/1 m0 #1",
      "pointer 1 0 400,350 b-1/0 m0 #1",
      "leave 0",
    ]);
    input.stop(false);
  });

  it("a press never starts the browser's drag and drop where the surface is its target", () => {
    const el = element();
    const r = recorder();
    const input = new DomInput(el, r.sink, () => 0);
    const drag = new Event("dragstart", { bubbles: true, cancelable: true });
    el.dispatchEvent(drag);
    expect(drag.defaultPrevented).toBe(true);
    input.stop(false);
  });

  it("a move with a button down over a surface that holds no press is the press's surface's, in its own px", () => {
    const a = element();
    const b = element();
    b.style.top = "300px";
    const ra = recorder();
    const rb = recorder();
    const ia = new DomInput(a, ra.sink, () => 0);
    const ib = new DomInput(b, rb.sink, () => 1);
    a.dispatchEvent(ptr("pointerdown", 120, 80, { buttons: 1, button: 0 }));
    b.dispatchEvent(ptr("pointermove", 150, 320, { buttons: 1 }));
    document.dispatchEvent(ptr("pointerup", 150, 320, { button: 0 }));
    expect(rb.log.filter((l) => l.startsWith("pointer"))).toEqual([]);
    expect(ra.log.filter((l) => l.startsWith("pointer"))).toEqual([
      "pointer 0 0 20,30 b0/1 m0 #1",
      "pointer 1 0 50,270 b-1/1 m0 #1",
      "pointer 2 0 50,270 b0/0 m0 #1",
    ]);
    // Without a button, a move over it is its hover.
    b.dispatchEvent(ptr("pointermove", 150, 320, { buttons: 0 }));
    expect(rb.log.filter((l) => l.startsWith("pointer"))).toEqual(["pointer 1 1 50,20 b-1/0 m0 #1"]);
    ia.stop(false);
    ib.stop(false);
  });

  it("wheel deltas in CSS px (lines and pages converted), the page's default prevented only when handled", () => {
    const el = element();
    for (const handled of [false, true]) {
      const r = recorder(handled);
      const input = new DomInput(el, r.sink, () => 1);
      const px = new WheelEvent("wheel", { deltaX: 3, deltaY: 40, ctrlKey: true, cancelable: true, bubbles: true });
      const lines = new WheelEvent("wheel", { deltaY: 2, deltaMode: 1, cancelable: true, bubbles: true });
      const page = new WheelEvent("wheel", { deltaY: 1, deltaMode: 2, cancelable: true, bubbles: true });
      for (const e of [px, lines, page]) el.dispatchEvent(e);
      expect(r.log.map((l) => l.replace(/ -?\d+,-?\d+ /, " "))).toEqual([
        `wheel 1 3,40 m${Modifier.ctrl}`,
        "wheel 1 0,32 m0",
        "wheel 1 0,200 m0",
      ]);
      expect([px, lines, page].map((e) => e.defaultPrevented)).toEqual(Array(3).fill(handled));
      expect(r.seen.wheel.size).toBe(1);
      input.stop(false);
    }
  });

  it("keys: the key and the physical code, modifiers, whether the page is typing, the default prevented when handled", () => {
    const el = element();
    const box = document.createElement("input");
    el.appendChild(box);
    for (const handled of [false, true]) {
      const r = recorder(handled);
      const input = new DomInput(el, r.sink, () => 2);
      const key = new KeyboardEvent("keydown", {
        key: "=",
        code: "Equal",
        altKey: true,
        cancelable: true,
        bubbles: true,
      });
      el.dispatchEvent(key);
      const typed = new KeyboardEvent("keydown", { key: "a", code: "KeyA", cancelable: true, bubbles: true });
      box.dispatchEvent(typed);
      expect(r.log).toEqual([`key = Equal m${Modifier.alt} typing:false`, "key a KeyA m0 typing:true"]);
      expect([key.defaultPrevented, typed.defaultPrevented]).toEqual([handled, handled]);
      expect(r.seen.key.size).toBe(1);
      input.stop(false);
    }
  });

  it("focus and blur reach the dispatcher; a press focuses the surface", () => {
    const el = element();
    const r = recorder();
    const input = new DomInput(el, r.sink, () => 4);
    el.dispatchEvent(ptr("pointerdown", 120, 80, { buttons: 1, button: 0 }));
    document.dispatchEvent(ptr("pointerup", 120, 80));
    expect(document.activeElement).toBe(el);
    el.blur();
    expect(r.log.filter((l) => l.startsWith("focus"))).toEqual(["focus 4 true", "focus 4 false"]);
    input.stop(false);
  });

  it("a stopped listener hears nothing; stopping cancels a press and ends a focus it held", () => {
    const el = element();
    const r = recorder();
    const input = new DomInput(el, r.sink, () => 5);
    el.dispatchEvent(ptr("pointerdown", 120, 80, { buttons: 1, button: 0 }));
    r.log.length = 0;
    input.stop(true);
    expect(r.log).toEqual(["cancel 5", "focus 5 false"]);
    el.dispatchEvent(ptr("pointermove", 125, 85));
    document.dispatchEvent(ptr("pointerup", 125, 85));
    el.dispatchEvent(new KeyboardEvent("keydown", { key: "x" }));
    el.dispatchEvent(new WheelEvent("wheel", { deltaY: 1 }));
    expect(r.log.length).toBe(2);
  });
});

describe("through a mounted surface", () => {
  let t: Opened | null = null;
  afterEach(() => {
    if (t) close(t);
    t = null;
  });

  it("a drag pans the viewport, a press focuses the surface, and the cursor behavior sets the element's cursor", async () => {
    t = await opened(syntheticFamily(40, 8), { width: 700, height: 300 });
    const o = t;
    const vp = o.ctx.viewports.get(o.section);
    if (!vp) throw new Error("no viewport");
    vp.setView(14, 0, "api");
    settle(o);
    const l = o.ctx.layout.get(o.section)?.current;
    if (!l) throw new Error("no layout");
    const box = o.surface.element.getBoundingClientRect();
    const x = box.left + l.dataLeft + 120;
    const y = box.top + l.scrolledTop + 40;
    const before = vp.live.current.scrollX;
    o.surface.element.dispatchEvent(ptr("pointermove", x, y));
    o.surface.element.dispatchEvent(ptr("pointerdown", x, y, { buttons: 1, button: 0 }));
    for (let k = 1; k <= 6; k++) document.dispatchEvent(ptr("pointermove", x - 20 * k, y, { buttons: 1 }));
    document.dispatchEvent(ptr("pointerup", x - 120, y, { button: 0 }));
    settle(o);
    expect(vp.live.current.scrollX).toBeGreaterThan(before + 60);
    expect(document.activeElement).toBe(o.surface.element);
    // The cursor behavior (set with the hover) reaches the element.
    o.ctx.ui.setCursor(o.section, "grab");
    expect(o.surface.element.style.cursor).toBe("grab");
  });

  it("focus shows the keyboard's cursor and blur clears it; disposing a focused surface clears it too", async () => {
    t = await opened(syntheticFamily(40, 8), { width: 700, height: 300 });
    const o = t;
    o.surface.element.focus();
    const key = (k: string) => o.surface.element.dispatchEvent(new KeyboardEvent("keydown", { key: k, bubbles: true }));
    key("ArrowRight");
    expect(o.ctx.ui.keyCursor.value).not.toBeNull();
    o.surface.element.blur();
    expect(o.ctx.ui.keyCursor.value).toBeNull();
    o.surface.element.focus();
    key("ArrowRight");
    expect(o.ctx.ui.keyCursor.value).not.toBeNull();
    o.surface.dispose();
    expect(o.ctx.ui.keyCursor.value).toBeNull();
  });

  it("a text selection drag goes on outside the surface (the document, another element) and commits its release there", async () => {
    const b = createSeqViewSpec();
    const a = b.data.sequence({ ref: "a", value: "ACDEFGHIKL".repeat(30) });
    b.view({ ref: "v", label: "V", layout: "text", height: "box" }).block({ sequence: a }).track().letters();
    const other = document.createElement("div");
    other.style.cssText = "position:absolute;left:0;top:400px;width:900px;height:200px";
    document.body.appendChild(other);
    hosts.push(other);
    t = await opened(b.getDocument(), { width: 500, height: 300 });
    const o = t;
    const l = o.ctx.layout.get(o.section)?.current;
    const ax = o.ctx.layout.axis(o.section);
    const run = ax?.runs[0];
    if (!l || !ax || !run) throw new Error("no layout");
    const box = o.surface.element.getBoundingClientRect();
    // Position 5 on line 0 (its letters' row, below the group labels).
    const x = box.left + l.dataLeft + runX(run, 0, 5, ax.gap) + run.advance / 2;
    const y = box.top + l.scrolledTop + l.entryTop[1] + l.entryHeight[1] / 2;
    const el = o.surface.element;
    el.dispatchEvent(ptr("pointerdown", x, y, { buttons: 1, button: 0 }));
    el.dispatchEvent(ptr("pointermove", x + 20, y, { buttons: 1 }));
    // Out of the surface: right of the text, then over another element below it; released on the document there.
    document.body.dispatchEvent(ptr("pointermove", box.right + 200, y, { buttons: 1 }));
    const pd = o.ctx.selection.pending.current;
    expect([o.ctx.ui.drag.value?.owner, pd.active, pd.endShown]).toEqual(["select", true, l.lines.stop[0] - 1]);
    other.dispatchEvent(ptr("pointermove", x, box.bottom + 150, { buttons: 1 }));
    document.dispatchEvent(ptr("pointerup", x, box.bottom + 150, { button: 0 }));
    const s = o.ctx.selection.current.value.seq;
    const xs: number[] = [];
    if (s) OrderedSet.forEach(s.elements[0].indices, (i: number) => xs.push(i), null);
    // Below the rows: down to the last line in view's position under x, or further as the rows scrolled.
    expect([o.ctx.ui.drag.value, pd.active, xs[0], xs.length > l.lines.stop[0]]).toEqual([null, false, 5, true]);
  });

  it("a real pointer's drag selects its range though an ancestor captures the pointer at its own press", async () => {
    const b = createSeqViewSpec();
    const a = b.data.sequence({ ref: "a", value: "ACDEFGHIKL".repeat(30) });
    b.view({ ref: "v", label: "V", layout: "text", height: "box" }).block({ sequence: a }).track().letters();
    const parent = document.createElement("div");
    document.body.appendChild(parent);
    hosts.push(parent);
    // A host that takes the pointer at a press (a draggable card, a gesture library).
    parent.addEventListener("pointerdown", (e) => parent.setPointerCapture(e.pointerId));
    t = await opened(b.getDocument(), { width: 500, height: 300, parent });
    const o = t;
    const l = o.ctx.layout.get(o.section)?.current;
    const ax = o.ctx.layout.axis(o.section);
    const run = ax?.runs[0];
    if (!l || !ax || !run) throw new Error("no layout");
    const at = (k: number, line: number) => ({
      x: l.dataLeft + runX(run, l.lines.start[line], k, ax.gap) + run.advance / 2,
      y: l.scrolledTop + l.entryTop[2 * line + 1] + l.entryHeight[2 * line + 1] / 2,
    });
    const from = at(5, 0);
    const to = at(l.lines.start[1] + 3, 1);
    await userEvent.dragAndDrop(o.surface.element, o.surface.element, {
      sourcePosition: from,
      targetPosition: to,
    });
    const s = o.ctx.selection.current.value.seq;
    const xs: number[] = [];
    if (s) OrderedSet.forEach(s.elements[0].indices, (i: number) => xs.push(i), null);
    expect([xs[0], xs.at(-1), xs.length]).toEqual([5, l.lines.start[1] + 3, l.lines.start[1] + 3 - 4]);
  });

  it("disposing the surface a press started on cancels the press: the drag ends and ui.drag clears", async () => {
    t = await opened(syntheticFamily(40, 8), { width: 700, height: 300 });
    const o = t;
    o.ctx.viewports.get(o.section)?.setView(14, 0, "api");
    settle(o);
    const l = o.ctx.layout.get(o.section)?.current;
    if (!l) throw new Error("no layout");
    const box = o.surface.element.getBoundingClientRect();
    const x = box.left + l.dataLeft + 120;
    const y = box.top + l.scrolledTop + 40;
    o.surface.element.dispatchEvent(ptr("pointerdown", x, y, { buttons: 1, button: 0 }));
    document.dispatchEvent(ptr("pointermove", x - 40, y, { buttons: 1 }));
    expect(o.ctx.ui.drag.value).not.toBeNull();
    o.surface.dispose();
    expect(o.ctx.ui.drag.value).toBeNull();
    expect(PointerType.cancel).toBe(3);
  });
});
