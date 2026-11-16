/** @vitest-environment jsdom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { afterEach, describe, expect, it, vi } from 'vitest';
import { Motion } from '../../dom/index.js';
import { installTransitions, milliseconds } from '../transitions.js';

describe('transitions', () => {
  const saved = Motion.close;
  afterEach(() => {
    Motion.close = saved;
    document.body.replaceChildren();
  });

  it('reads CSS times', () => {
    expect(milliseconds('100ms')).toBe(100);
    expect(milliseconds(' 0.2s ')).toBe(200);
    expect(milliseconds('0ms')).toBe(0);
    expect(milliseconds('')).toBe(0);
  });

  it('closes a fold by shrinking it over the base duration, a fade over the exit; none at 0 (reduced motion)', () => {
    Motion.close = undefined;
    installTransitions();
    const el = document.createElement('div');
    el.style.setProperty('--msui-motion-base', '200ms');
    el.style.setProperty('--msui-motion-fast', '100ms');
    document.body.append(el);
    const calls: [Keyframe[], KeyframeAnimationOptions][] = [];
    (el as any).animate = vi.fn((frames: Keyframe[], options: KeyframeAnimationOptions) => {
      calls.push([frames, options]);
      return {} as Animation;
    });
    expect(Motion.close!(el, 'grow')).toBeDefined();
    expect(Motion.close!(el, 'fade')).toBeDefined();
    expect(calls.map(([f, o]) => [f.at(-1)!.height ?? 'same', o.duration])).toEqual([
      ['0px', 200],
      ['same', 100],
    ]);
    el.style.setProperty('--msui-motion-fast', '0ms');
    expect(Motion.close!(el, 'fade')).toBeUndefined();
  });

  it('keeps a closing an application installed', () => {
    const own = () => undefined;
    Motion.close = own;
    installTransitions();
    expect(Motion.close).toBe(own);
  });
});
