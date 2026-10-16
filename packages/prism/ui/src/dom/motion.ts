/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Motion (Overview › Foundations › Motion): content opens with `data-reveal` (reset.css: it grows to its height and
// fades in, slides in from the side, or only fades in) and closes with `playClose`, timed by the motion tokens; nothing
// moves with reduced motion.

/** How content is revealed: grown to its height, slid in from the side, or only faded (a band drawn on a canvas). */
export type RevealKind = 'grow' | 'side' | 'fade';

/**
 * Plays a reveal back on `el` (collapse and fade, slide out to the side, or fade), timed by the motion tokens;
 * undefined where it should not or cannot run (reduced motion, no Web Animations, not on screen), so the content just
 * goes.
 */
export function playClose(el: Element, kind: RevealKind): Animation | undefined {
  if (typeof (el as HTMLElement).animate !== 'function' || !el.isConnected) return undefined;
  if (globalThis.matchMedia?.('(prefers-reduced-motion: reduce)').matches) return undefined;
  const style = getComputedStyle(el);
  const easing = style.getPropertyValue('--msp-ease-reveal').trim() || 'ease';
  if (kind === 'fade') {
    const duration = toMs(style.getPropertyValue('--msp-motion-fade')) || 360;
    return (el as HTMLElement).animate([{ opacity: 1 }, { opacity: 0 }], { duration, easing, fill: 'forwards' });
  }
  const duration = toMs(style.getPropertyValue('--msp-motion-base')) || 220;
  const frames =
    kind === 'side'
      ? [
          { opacity: 1, translate: '0 0' },
          { opacity: 0, translate: '12px 0' },
        ]
      : [
          { height: `${el.getBoundingClientRect().height}px`, opacity: 1, overflow: 'clip' },
          // The whole box to nothing (padding and margins too, or they are left over and the end jumps).
          { height: '0px', paddingBlock: '0px', marginBlock: '0px', opacity: 0, overflow: 'clip' },
        ];
  return (el as HTMLElement).animate(frames, { duration, easing, fill: 'forwards' });
}

/**
 * Runs `update` (a change of what `scope` shows: one panel replacing another) as a view transition scoped to `scope`
 * where the browser has them (`Element.startViewTransition`): the old content cross-fades into the new. Elsewhere, and
 * under reduced motion, it runs at once. Never the document's view transition: it snapshots the whole page and the 3D
 * canvas would stand still while it runs.
 */
export function swapIn(scope: Element | null | undefined, update: () => void): void {
  const start = (scope as (Element & { startViewTransition?: (cb: () => void) => unknown }) | null | undefined)
    ?.startViewTransition;
  if (
    !scope?.isConnected ||
    typeof start !== 'function' ||
    globalThis.matchMedia?.('(prefers-reduced-motion: reduce)').matches
  ) {
    update();
    return;
  }
  start.call(scope, update);
}

/** A CSS time (`180ms`, `.18s`: the CSS build writes seconds) in milliseconds. */
function toMs(value: string): number {
  const v = value.trim();
  const n = Number.parseFloat(v);
  return v.endsWith('ms') ? n : v.endsWith('s') ? n * 1000 : n;
}
