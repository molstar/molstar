/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { DomHost } from '@molstar/ui/dom/host';

// The headers' DOM gestures on the mounted canvas (_spec/README.md "Behaviors", Headers: Gestures): the secondary click
// (`contextmenu`, prevented when the menu took it: Shift's is the browser's, the behavior says) and the long press of a
// finger (touch's context menu), both in canvas coordinates. Everything it adds is taken away by the function it returns.

/** A long press's slop, CSS px. */
const PressSlop = 8;

/**
 * @internal Listens for the secondary click and the long press on `container`: `context(x, y, e)` answers whether the menu took
 * the click; `pressed(x, y)` runs after `longPressMs` of a finger that did not move.
 */
export function headerGestures(
  container: HTMLElement,
  longPressMs: number,
  context: (x: number, y: number, e: MouseEvent) => boolean,
  pressed: (x: number, y: number) => void,
): () => void {
  const local = (e: MouseEvent) => {
    const r = container.getBoundingClientRect();
    return [e.clientX - r.left, e.clientY - r.top] as const;
  };
  const onContext = (e: MouseEvent) => {
    const [x, y] = local(e);
    if (context(x, y, e)) e.preventDefault();
  };
  let press: { x: number; y: number; timer: unknown } | undefined;
  const endPress = () => {
    if (press) DomHost.clearTimeout(press.timer);
    press = undefined;
  };
  const onDown = (e: PointerEvent) => {
    endPress();
    if (e.pointerType !== 'touch') return;
    const [x, y] = local(e);
    press = {
      x,
      y,
      timer: DomHost.setTimeout(() => {
        press = undefined;
        pressed(x, y);
      }, longPressMs),
    };
  };
  const onMove = (e: PointerEvent) => {
    if (!press) return;
    const [x, y] = local(e);
    if (Math.max(Math.abs(x - press.x), Math.abs(y - press.y)) > PressSlop) endPress();
  };
  const listeners: [string, (e: any) => void][] = [
    ['contextmenu', onContext],
    ['pointerdown', onDown],
    ['pointermove', onMove],
    ['pointerup', endPress],
    ['pointercancel', endPress],
  ];
  for (const [type, fn] of listeners) container.addEventListener(type, fn);
  return () => {
    endPress();
    for (const [type, fn] of listeners) container.removeEventListener(type, fn);
  };
}
