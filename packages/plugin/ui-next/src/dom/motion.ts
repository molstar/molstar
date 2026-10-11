/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Motion (Overview › Foundations › Motion): content opens with `data-reveal` (reset.css: it grows to its height and
// fades in, or slides in from the side) and closes with `playClose`, timed by the motion tokens; nothing moves with
// reduced motion.

/**
 * Plays a reveal back on `el` (collapse and fade, or slide out to the side), timed by the motion tokens; undefined
 * where it should not or cannot run (reduced motion, no Web Animations, not on screen), so the content just goes.
 */
export function playClose(el: Element, side: boolean): Animation | undefined {
  if (typeof (el as HTMLElement).animate !== 'function' || !el.isConnected) return undefined;
  if (globalThis.matchMedia?.('(prefers-reduced-motion: reduce)').matches) return undefined;
  const style = getComputedStyle(el);
  const duration = toMs(style.getPropertyValue('--msp-motion-base')) || 180;
  const easing = style.getPropertyValue('--msp-ease').trim() || 'ease';
  const frames = side
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

/** A CSS time (`180ms`, `.18s`: the CSS build writes seconds) in milliseconds. */
function toMs(value: string): number {
  const v = value.trim();
  const n = Number.parseFloat(v);
  return v.endsWith('ms') ? n : v.endsWith('s') ? n * 1000 : n;
}
