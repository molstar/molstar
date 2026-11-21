/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { modifiersOf } from '@molstar/seqx-viewer/bindings';
import { createDragSession } from '@molstar/seqx-viewer/drag-session';
import type { NavigatorModel } from '@molstar/seqx-viewer/navigator';
import type { SeqViewer } from '@molstar/seqx-viewer/viewer';
import { cancelFrameCallback, requestFrameCallback } from '@molstar/ui/dom/frame';
import { hostPicture } from './compose.js';

// The host-drawn navigator's DOM (the viewer's _spec "Models", NavigatorModel, the contract): the host's element drawn
// from the model's last frame (its domains and windows as elements, placed by `xRange`) and each of its events one of
// the model's methods; a drag ends wherever the button is released (the viewer's drag session, as its own part does).

/** Draws navigator `model` into `nav` (the host's element, `height` px) and maps its events; returns its disposer. */
export function wireHostNavigator(
  nav: HTMLElement,
  viewer: SeqViewer,
  model: NavigatorModel,
  height: number,
  readout?: HTMLElement,
): () => void {
  const doc = nav.ownerDocument;
  nav.tabIndex = 0;
  nav.setAttribute('role', 'slider');
  nav.setAttribute('aria-label', 'Host-drawn navigator');
  Object.assign(nav.style, { position: 'relative', touchAction: 'none', userSelect: 'none' });
  let frame: number | undefined;
  const paint = () => {
    frame = undefined;
    model.resize(nav.clientWidth, height);
    if (!model.draw()) return;
    const pic = hostPicture(model);
    const els: HTMLElement[] = [];
    for (const b of pic.blocks) {
      const el = doc.createElement('div');
      el.textContent = b.label;
      el.dataset.block = b.label;
      el.style.cssText = `position: absolute; left: ${b.x}px; width: ${b.width}px; top: 8px; bottom: 8px; border-radius: 3px; background: color-mix(in srgb, var(--msui-accent) 25%, transparent); color: var(--msui-text); font-size: 10px; line-height: ${height - 16}px; text-align: center; overflow: hidden; pointer-events: none;`;
      els.push(el);
    }
    for (const w of pic.windows) {
      const el = doc.createElement('div');
      el.dataset.window = w.section;
      el.style.cssText = `position: absolute; left: ${w.x}px; width: ${w.width}px; top: 1px; bottom: 1px; box-sizing: border-box; border: 2px solid ${w.idle ? 'var(--msui-line)' : 'var(--msui-accent)'}; border-radius: 4px; background: color-mix(in srgb, var(--msui-accent) 10%, transparent); pointer-events: none;`;
      els.push(el);
    }
    nav.replaceChildren(...els);
    nav.style.cursor = model.cursor;
    if (readout) readout.textContent = model.readout.value ?? 'Hover the navigator';
  };
  const schedule = () => {
    if (frame === undefined) frame = requestFrameCallback(paint);
  };
  const say = () => nav.setAttribute('aria-valuetext', model.values.value[model.focused.value]?.text ?? '');
  const offs = [model.version.listen(schedule), viewer.drawn.listen(schedule), model.values.listen(say)];

  const at = (e: MouseEvent) => {
    const r = nav.getBoundingClientRect();
    return [e.clientX - r.left, e.clientY - r.top] as const;
  };
  const touch = (e: PointerEvent) => ({ touch: e.pointerType === 'touch' });
  const drag = createDragSession(nav, () => model.up());
  const listeners: [string, (e: any) => void, AddEventListenerOptions?][] = [
    ['pointermove', (e: PointerEvent) => model.move(...at(e), touch(e))],
    ['pointerleave', () => model.leave()],
    [
      'pointerdown',
      (e: PointerEvent) => {
        if (e.button !== 0) return;
        drag.end();
        if (model.down(...at(e), touch(e))) {
          // no text selection or native drag starts under the hand
          e.preventDefault();
          nav.setPointerCapture?.(e.pointerId);
          drag.start();
        }
        nav.focus({ preventScroll: true });
      },
    ],
    [
      'pointerup',
      (e: PointerEvent) => {
        drag.end();
        model.up();
        if (nav.hasPointerCapture?.(e.pointerId)) nav.releasePointerCapture(e.pointerId);
      },
    ],
    ['pointercancel', () => (drag.end(), model.cancel())],
    ['click', (e: MouseEvent) => e.button === 0 && model.click(...at(e), e.detail)],
    ['dblclick', (e: MouseEvent) => model.doubleClick(...at(e))],
    ['focus', () => model.setFocus(true)],
    ['blur', () => model.setFocus(false)],
    [
      'wheel',
      (e: WheelEvent) => {
        if (model.wheel(...at(e), e.deltaX, e.deltaY, modifiersOf(e, viewer.input.mac))) e.preventDefault();
      },
      { passive: false },
    ],
    ['keydown', (e: KeyboardEvent) => model.key(e.key, e) && e.preventDefault()],
  ];
  for (const [type, fn, options] of listeners) nav.addEventListener(type, fn, options);
  const resized = typeof ResizeObserver === 'undefined' ? undefined : new ResizeObserver(schedule);
  resized?.observe(nav);
  schedule();
  return () => {
    for (const off of offs) off();
    for (const [type, fn, options] of listeners) nav.removeEventListener(type, fn, options);
    resized?.disconnect();
    drag.dispose();
    if (frame !== undefined) cancelFrameCallback(frame);
  };
}
