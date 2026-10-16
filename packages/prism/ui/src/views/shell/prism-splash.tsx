/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { bind, type MaybeReactive, onActive, onCleanup } from '../../dom/scope.js';
import { classes } from './_generated/prism-splash.module.css.js';
import { PrismField, type PrismFieldOptions } from './prism-field.js';

export type { PrismFieldOptions } from './prism-field.js';
export { PrismFieldDefaults } from './prism-field.js';

function prefersReducedMotion(): boolean {
  return typeof matchMedia === 'function' && matchMedia('(prefers-reduced-motion: reduce)').matches;
}

/**
 * Light through the Mol* star (`prism-field.ts`): a canvas filling what holds it, the light splitting at the star of the
 * wordmark `mark` finds (its `[data-part="star"]`), else at the middle. It answers the pointer over `pointer` (the
 * viewport: the canvas keeps its clicks and drags, the splash only listens), mouse, pen and touch alike. It runs while it
 * is on the screen, the page is visible and its view active; under reduced motion it is one still frame.
 */
export function PrismSplash(props: {
  options?: MaybeReactive<Partial<PrismFieldOptions>>;
  /** The element holding the wordmark the light goes through. */
  mark?: () => Element | null | undefined;
  /** Where the pointer is listened to (passively); none: the splash does not answer it. */
  pointer?: () => EventTarget | null | undefined;
  class?: string;
}): Node {
  const canvas = (
    <canvas class={props.class ? `${classes.splash} ${props.class}` : classes.splash} aria-hidden="true" />
  ) as HTMLCanvasElement;
  const field = new PrismField(canvas, {}, { still: prefersReducedMotion() });
  if (props.options) bind(props.options, (o) => field.setOptions(o));

  // The canvas's place on the page, read when it resizes or the page scrolls (not per pointer event: no layout reads at
  // pointer rate).
  let left = 0;
  let top = 0;
  const locate = () => {
    const box = canvas.getBoundingClientRect();
    left = box.left;
    top = box.top;
  };
  const place = () => {
    field.resize();
    const box = canvas.getBoundingClientRect();
    left = box.left;
    top = box.top;
    const star = props.mark?.()?.querySelector('[data-part="star"]');
    const at = star?.getBoundingClientRect();
    if (at?.width) field.setOrigin(at.left + at.width / 2 - box.left, at.top + at.height / 2 - box.top);
    else field.setOrigin(box.width / 2, box.height / 2);
  };

  let started = false;
  const sizes =
    typeof ResizeObserver === 'function'
      ? new ResizeObserver(() => {
          place();
          if (!started) {
            started = true;
            attachPointer();
            field.start();
          }
        })
      : undefined;
  sizes?.observe(canvas);

  // Shown, the page visible, the view active (a hidden shell part pauses its scope).
  let onScreen = true;
  let active = true;
  const update = () => field.setVisible(onScreen && active && document.visibilityState === 'visible');
  const seen =
    typeof IntersectionObserver === 'function'
      ? new IntersectionObserver(([e]) => {
          onScreen = e?.isIntersecting ?? true;
          update();
        })
      : undefined;
  seen?.observe(canvas);
  document.addEventListener('visibilitychange', update);
  addEventListener('scroll', locate, { capture: true, passive: true });
  onActive((a) => {
    active = a;
    update();
  });

  // The pointer, listened to on the host passively: the canvas below keeps every event.
  let detachPointer = () => {};
  function attachPointer() {
    const host = props.pointer?.();
    if (!host) return;
    const move = (e: Event) => {
      const p = e as PointerEvent;
      field.pointerMove(p.clientX - left, p.clientY - top, p.timeStamp);
    };
    const down = (e: Event) => {
      const p = e as PointerEvent;
      field.pointerDown(p.clientX - left, p.clientY - top, p.timeStamp);
    };
    const leave = (e: Event) => {
      // A touch lifts as it ends; a mouse only when it leaves.
      if (e.type !== 'pointerup' || (e as PointerEvent).pointerType !== 'mouse') field.pointerLeave();
    };
    const opts = { passive: true } as const;
    host.addEventListener('pointermove', move, opts);
    host.addEventListener('pointerdown', down, opts);
    host.addEventListener('pointerleave', leave, opts);
    host.addEventListener('pointerup', leave, opts);
    host.addEventListener('pointercancel', leave, opts);
    detachPointer = () => {
      host.removeEventListener('pointermove', move);
      host.removeEventListener('pointerdown', down);
      host.removeEventListener('pointerleave', leave);
      host.removeEventListener('pointerup', leave);
      host.removeEventListener('pointercancel', leave);
    };
  }

  onCleanup(() => {
    sizes?.disconnect();
    seen?.disconnect();
    document.removeEventListener('visibilitychange', update);
    removeEventListener('scroll', locate, { capture: true });
    detachPointer();
    field.dispose();
  });
  return canvas;
}
