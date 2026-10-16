/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { bind, type MaybeReactive, onCleanup } from '../../dom/scope.js';
import { classes } from './_generated/prism-splash.module.css.js';
import { PrismField, type PrismFieldOptions } from './prism-field.js';

export type { PrismFieldOptions } from './prism-field.js';
export { PrismFieldDefaults } from './prism-field.js';

function prefersReducedMotion(): boolean {
  return typeof matchMedia === 'function' && matchMedia('(prefers-reduced-motion: reduce)').matches;
}

/**
 * Light through Prism's mark (`prism-field.ts`): a canvas filling what holds it, the light splitting at the prism of the
 * wordmark `mark` finds (its `[data-part="prism"]`), else at the middle. It runs while it is on the screen and the page
 * is visible; under reduced motion it is one still frame.
 */
export function PrismSplash(props: {
  options?: MaybeReactive<Partial<PrismFieldOptions>>;
  /** The element holding the wordmark the light goes through. */
  mark?: () => Element | null | undefined;
  class?: string;
}): Node {
  const canvas = (
    <canvas class={props.class ? `${classes.splash} ${props.class}` : classes.splash} aria-hidden="true" />
  ) as HTMLCanvasElement;
  const field = new PrismField(canvas, {}, { still: prefersReducedMotion() });
  if (props.options) bind(props.options, (o) => field.setOptions(o));

  // The split point: the middle of the prism's triangle, in the canvas's pixels.
  const place = () => {
    field.resize();
    const box = canvas.getBoundingClientRect();
    const prism = props.mark?.()?.querySelector('[data-part="prism"]');
    const at = prism?.getBoundingClientRect();
    if (at && at.width) field.setOrigin(at.left + at.width / 2 - box.left, at.top + at.height / 2 - box.top);
    else field.setOrigin(box.width / 2, box.height / 2);
  };

  let started = false;
  const sizes =
    typeof ResizeObserver === 'function'
      ? new ResizeObserver(() => {
          place();
          if (!started) {
            started = true;
            field.start();
          }
        })
      : undefined;
  sizes?.observe(canvas);
  let onScreen = true;
  const update = () => field.setVisible(onScreen && document.visibilityState === 'visible');
  const seen =
    typeof IntersectionObserver === 'function'
      ? new IntersectionObserver(([e]) => {
          onScreen = e?.isIntersecting ?? true;
          update();
        })
      : undefined;
  seen?.observe(canvas);
  document.addEventListener('visibilitychange', update);
  onCleanup(() => {
    sizes?.disconnect();
    seen?.disconnect();
    document.removeEventListener('visibilitychange', update);
    field.dispose();
  });
  return canvas;
}
