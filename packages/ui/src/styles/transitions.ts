/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Motion, type RevealKind } from '../dom/index.js';

/**
 * The UI's closings (styles §7.1), installed into the runtime's `Motion` once (with the styles: `adoptStyles`): content
 * that opened with `Show`'s `reveal` plays it back as it goes. A fold shrinks to nothing over `--msui-motion-base`
 * (eased in and out, as it opened); a slid or faded one fades out over `--msui-exit` (fast, eased in). Durations come
 * from the tokens on the element, so reduced motion (every duration 0) or no Web Animations means no animation: the
 * content goes at once. An application that sets `Motion.close` itself keeps its own.
 */
export function installTransitions(): void {
  if (Motion.close) return;
  Motion.close = playClosing;
}

function playClosing(el: Element, kind: RevealKind): Animation | undefined {
  if (typeof (el as HTMLElement).animate !== 'function' || typeof getComputedStyle !== 'function') return undefined;
  const style = getComputedStyle(el);
  const grow = kind === 'grow';
  const duration = milliseconds(style.getPropertyValue(grow ? '--msui-motion-base' : '--msui-motion-fast'));
  if (!duration) return undefined;
  const easing = style.getPropertyValue(grow ? '--msui-ease-in-out' : '--msui-ease-in').trim() || 'ease-in';
  const frames: Keyframe[] = grow
    ? [
        { height: `${(el as HTMLElement).offsetHeight}px`, overflow: 'clip' },
        { height: '0px', paddingBlock: '0px', marginBlock: '0px', opacity: 0, overflow: 'clip' },
      ]
    : kind === 'side'
      ? [{}, { opacity: 0, translate: '12px 0' }]
      : [{}, { opacity: 0 }];
  return (el as HTMLElement).animate(frames, { duration, easing });
}

/** A CSS time in ms (`100ms`, `0.2s`; the duration tokens are registered times: they read back resolved); 0 for none. */
export function milliseconds(time: string): number {
  const t = time.trim();
  if (t.endsWith('ms')) return Number.parseFloat(t) || 0;
  if (t.endsWith('s')) return (Number.parseFloat(t) || 0) * 1000;
  return 0;
}
