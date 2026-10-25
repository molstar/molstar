/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { bind, type MaybeReactive, onActive, onCleanup, requireScope, runInScope } from '../../dom/scope.js';
import { observeIntersection, observeResize, on } from '../../dom/scoped.js';
import { classes } from './_generated/prism-splash.module.css.js';
import { WordmarkGeometry } from '../../components/icon/wordmark.js';
import { PrismField, type PrismFieldOptions } from './prism-field.js';
import { QuantumField, type QuantumFieldOptions } from './quantum-field.js';

export type { PrismFieldOptions } from './prism-field.js';
export { PrismFieldDefaults } from './prism-field.js';
export type { QuantumFieldOptions } from './quantum-field.js';
export { QuantumFieldDefaults } from './quantum-field.js';

/** What the splash runs: the light (`prism-field.ts`) or the light as a wave (`quantum-field.ts`). */
export type SplashKind = 'light' | 'quantum';

function prefersReducedMotion(): boolean {
  return typeof matchMedia === 'function' && matchMedia('(prefers-reduced-motion: reduce)').matches;
}

/**
 * Light through the Mol* star (`prism-field.ts`; or as a wave, `quantum-field.ts`): a canvas filling what holds it, the light splitting at the star of the
 * wordmark `mark` finds (its `[data-part="star"]`), else at the middle. The canvas draws the wordmark, as glass
 * (`glass.ts`): it assembles from pieces as the splash starts, its pieces lean under the pointer, and a press on it
 * shatters it for a moment; the page's wordmark only places it. It answers the pointer over `pointer` (the
 * viewport: the canvas keeps its clicks and drags, the splash only listens), mouse, pen and touch alike. It runs while it
 * is on the screen, the page is visible and its view active; at rest (no input a while, the glass still, after its
 * intro) the light fades out and it asks for no frame, the word staying drawn, until the pointer or a resize wakes it
 * (`ambient-pace.ts`); under reduced motion it is one still frame.
 */
export function PrismSplash(props: {
  /** The light (default) or the light as a wave; `options` are that field's. */
  field?: SplashKind;
  options?: MaybeReactive<Partial<PrismFieldOptions> | Partial<QuantumFieldOptions>>;
  /** The element holding the wordmark the light goes through. */
  mark?: () => Element | null | undefined;
  /** Where the pointer is listened to (passively); none: the splash does not answer it. */
  pointer?: () => EventTarget | null | undefined;
  class?: string;
}): Node {
  const canvas = (
    <canvas class={props.class ? `${classes.splash} ${props.class}` : classes.splash} aria-hidden="true" />
  ) as HTMLCanvasElement;
  const still = { still: prefersReducedMotion() };
  const field = props.field === 'quantum' ? new QuantumField(canvas, {}, still) : new PrismField(canvas, {}, still);
  if (props.options) bind(props.options, (o) => field.setOptions(o as never));

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
    // The word as glass: where its drawing origin is (its viewBox runs from -1) and px a unit. While the canvas draws it,
    // the page's wordmark only places it (hidden: no hover colors, nothing to hand over to).
    const word = wordOf();
    const w = word?.getBoundingClientRect();
    if (word && w?.width) {
      const scale = w.width / WordmarkGeometry.width;
      if (field.setWord(w.left - box.left + scale, w.top - box.top + scale, scale)) word.style.visibility = 'hidden';
    }
    // Drawn now, in its new place: the resize left the canvas blank until the next frame (a blink).
    field.redraw();
  };
  // The wordmark itself (`role=img`), in or as what `mark` gives.
  const wordOf = (): HTMLElement | undefined => {
    const m = props.mark?.();
    if (!m) return undefined;
    return (m.matches('[role=img]') ? m : m.querySelector('[role=img]')) as HTMLElement | undefined;
  };

  // What it attaches is scoped (dom/scoped.ts): released with it, a listener detached while it is paused.
  const scope = requireScope('PrismSplash');
  let started = false;
  observeResize(canvas, () => {
    place();
    if (!started) {
      started = true;
      // From the observer's callback, outside the scope: attached in it, once.
      if (!scope.disposed) runInScope(scope, attachPointer);
      field.start();
    }
  });

  // Shown, the page visible, the view active (a hidden shell part pauses its scope).
  let onScreen = true;
  let active = true;
  const update = () => field.setVisible(onScreen && active && document.visibilityState === 'visible');
  observeIntersection(canvas, ([e]) => {
    onScreen = e?.isIntersecting ?? true;
    update();
  });
  on(document, 'visibilitychange', update);
  on(window, 'scroll', locate, { capture: true, passive: true });
  onActive((a) => {
    active = a;
    update();
  });

  // The pointer, listened to on the host passively: the canvas below keeps every event.
  function attachPointer() {
    const host = props.pointer?.();
    if (!host) return;
    const opts = { passive: true } as const;
    on<PointerEvent>(
      host,
      'pointermove',
      (p) => field.pointerMove(p.clientX - left, p.clientY - top, p.timeStamp),
      opts,
    );
    on<PointerEvent>(
      host,
      'pointerdown',
      (p) => {
        const x = p.clientX - left;
        const y = p.clientY - top;
        // On the word: it shatters (the light with it); elsewhere, a burst of light.
        if (!field.shatterAt(x, y)) field.pointerDown(x, y, p.timeStamp);
      },
      opts,
    );
    const leave = (e: PointerEvent) => {
      // A touch lifts as it ends; a mouse only when it leaves.
      if (e.type !== 'pointerup' || e.pointerType !== 'mouse') field.pointerLeave();
    };
    on(host, 'pointerleave', leave, opts);
    on(host, 'pointerup', leave, opts);
    on(host, 'pointercancel', leave, opts);
  }

  onCleanup(() => {
    const word = wordOf();
    if (word) word.style.visibility = '';
    field.dispose();
  });
  return canvas;
}
