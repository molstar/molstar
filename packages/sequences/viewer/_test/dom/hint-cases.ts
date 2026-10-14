/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A hint never blocks the pointer (seqspec:interaction#accessibility), the same cases for both UIs (dom/hint.test.ts,
// react/hint.test.tsx):
// it shows after the rest delay, at once for the next trigger in a row and on keyboard focus; it hides when the pointer
// leaves, on a press and on a key; it takes no pointer events, so the control under it gets the hover and the click;
// and it takes the side that covers no control.

import { afterEach, expect, it } from "vitest";
import { userEvent } from "vitest/browser";
import { HINT_DELAY_MS } from "../../src/dom/shared.ts";
import "../../src/skin.css";

/** Mounts hints over `buttons` (one each, labeled "Hint <i>") into `host`; returns the unmount. */
export type MountHints = (host: HTMLElement, buttons: HTMLButtonElement[]) => Promise<() => Promise<void> | void>;

const shown = (b: HTMLElement) => {
  const bubble = b.closest(".seqv-hint")?.querySelector<HTMLElement>(".seqv-hint-bubble");
  return !!bubble && getComputedStyle(bubble).visibility === "visible";
};
const bubbleOf = (b: HTMLElement) => b.closest(".seqv-hint")?.querySelector(".seqv-hint-bubble") as HTMLElement;
const wait = (ms: number) => new Promise((r) => setTimeout(r, ms));
const button = (label: string) => {
  const b = document.createElement("button");
  b.type = "button";
  b.textContent = label;
  b.style.cssText = "font: 12px sans-serif; padding: 2px 6px";
  return b;
};

export function hintCases(mountHints: MountHints) {
  const cleanups: (() => Promise<void> | void)[] = [];
  afterEach(async () => {
    for (const c of cleanups.splice(0)) await c();
  });

  /**
   * A column: a plain control `above`, the triggers in a row under it (`row` px down), and `below` under them; the
   * hints are over `triggers`. `between` px between the rows.
   */
  async function stage(o: { count?: number; between?: number; below?: boolean } = {}) {
    const box = document.createElement("div");
    box.style.cssText = "position:absolute;left:200px;top:200px;display:flex;flex-direction:column;align-items:start";
    box.style.gap = `${o.between ?? 2}px`;
    const above = button("Above, a control under the hint");
    const row = document.createElement("div");
    row.style.cssText = "display:flex;gap:2px";
    const triggers = Array.from({ length: o.count ?? 1 }, (_, i) => button(`Trigger ${i}`));
    const below = button("Below, another control");
    box.append(above, row);
    if (o.below) box.append(below);
    document.body.append(box);
    const unmount = await mountHints(row, triggers);
    cleanups.push(async () => {
      await unmount();
      box.remove();
    });
    await userEvent.unhover(document.body);
    return { above, below, triggers };
  }

  it("shows after the rest delay, not before, and hides at once when the pointer leaves", async () => {
    const { triggers, above } = await stage();
    const [t] = triggers;
    await userEvent.hover(t);
    await wait(HINT_DELAY_MS / 3);
    expect(shown(t)).toBe(false);
    await expect.poll(() => shown(t), { timeout: HINT_DELAY_MS * 3 }).toBe(true);
    expect(bubbleOf(t).textContent).toBe("Hint 0");
    await userEvent.hover(above);
    expect(shown(t)).toBe(false);
  });

  it("takes no pointer events: over a control, the control gets the hover and the click", async () => {
    // A control above and below: no side is free, so the hint goes above, over `above`.
    const { triggers, above } = await stage({ below: true });
    const [t] = triggers;
    await userEvent.hover(t);
    await expect.poll(() => shown(t), { timeout: HINT_DELAY_MS * 3 }).toBe(true);
    const b = bubbleOf(t).getBoundingClientRect();
    const a = above.getBoundingClientRect();
    expect(b.bottom > a.top && b.top < a.bottom).toBe(true);
    // The hint's middle is over `above`: the browser's hit test passes through it.
    const x = Math.min(b.right, a.right) - 4;
    const y = (Math.max(b.top, a.top) + Math.min(b.bottom, a.bottom)) / 2;
    expect(document.elementFromPoint(x, y)).toBe(above);
    let entered = 0;
    let clicked = 0;
    above.addEventListener("pointerenter", () => entered++);
    above.addEventListener("click", () => clicked++);
    // A real pointer at that point (Playwright waits while another element would take the click).
    await userEvent.click(above, { position: { x: x - a.left, y: y - a.top } });
    expect([entered, clicked]).toEqual([1, 1]);
    expect(shown(t)).toBe(false);
  });

  it("covers no control when a side is free: below, when a control is right above", async () => {
    const { triggers, above } = await stage();
    const [t] = triggers;
    await userEvent.hover(t);
    await expect.poll(() => shown(t), { timeout: HINT_DELAY_MS * 3 }).toBe(true);
    const b = bubbleOf(t).getBoundingClientRect();
    expect(b.top).toBeGreaterThanOrEqual(t.getBoundingClientRect().bottom);
    expect(b.top).toBeGreaterThanOrEqual(above.getBoundingClientRect().bottom);
  });

  it("shows the next trigger's at once, moving along a row", async () => {
    const { triggers } = await stage({ count: 3 });
    await userEvent.hover(triggers[0]);
    await expect.poll(() => shown(triggers[0]), { timeout: HINT_DELAY_MS * 3 }).toBe(true);
    await userEvent.hover(triggers[1]);
    expect(shown(triggers[0])).toBe(false);
    await expect.poll(() => shown(triggers[1]), { timeout: HINT_DELAY_MS / 2 }).toBe(true);
  });

  it("hides on a press, and stays hidden while the pointer stays", async () => {
    const { triggers } = await stage();
    const [t] = triggers;
    await userEvent.hover(t);
    await expect.poll(() => shown(t), { timeout: HINT_DELAY_MS * 3 }).toBe(true);
    await userEvent.click(t);
    expect(shown(t)).toBe(false);
    await wait(HINT_DELAY_MS * 1.5);
    expect(shown(t)).toBe(false);
  });

  it("shows at once on keyboard focus, hides on a key, and when focus leaves", async () => {
    const { triggers } = await stage({ count: 2 });
    const before = button("Before");
    triggers[0].closest(".seqv-hint")?.parentElement?.prepend(before);
    before.focus();
    await userEvent.tab();
    expect(document.activeElement).toBe(triggers[0]);
    expect(shown(triggers[0])).toBe(true);
    await userEvent.keyboard("{Shift}");
    expect(shown(triggers[0])).toBe(false);
    await userEvent.tab();
    expect(shown(triggers[0])).toBe(false);
    expect(shown(triggers[1])).toBe(true);
    before.remove();
  });
}
