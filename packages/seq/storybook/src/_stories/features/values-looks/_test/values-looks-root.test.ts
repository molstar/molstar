// @vitest-environment jsdom
/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { DomHost } from '@molstar/ui/dom';
import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { AllCells } from '../compose.js';
import { Looks } from '../values-looks.stories.js';

// The Values look story's root rendered in a DOM (jsdom) as Storybook renders it: every cell mounts its viewer inside the
// element's render scope (an `onCleanup` outside one throws) and gets a canvas, in both themes.

let mounted: HTMLElement | undefined;

/** jsdom has no 2D context: one that takes every call and measures text by its length (nothing drawn is checked). */
function stubContext(canvas: HTMLCanvasElement): unknown {
  const target: Record<string | symbol, unknown> = {
    canvas,
    measureText: (t: string) => ({ width: 6 * t.length, actualBoundingBoxAscent: 8, actualBoundingBoxDescent: 2 }),
    getImageData: (_x: number, _y: number, w: number, h: number) => ({
      data: new Uint8ClampedArray(4 * w * h),
      width: w,
      height: h,
    }),
    createImageData: (w: number, h: number) => ({ data: new Uint8ClampedArray(4 * w * h), width: w, height: h }),
    getTransform: () => ({ a: 1, b: 0, c: 0, d: 1, e: 0, f: 0 }),
    createLinearGradient: () => ({ addColorStop() {} }),
    createRadialGradient: () => ({ addColorStop() {} }),
    createPattern: () => ({}),
  };
  return new Proxy(target, {
    get: (t, k) => (k in t ? t[k] : () => {}),
    set: (t, k, v) => {
      t[k] = v;
      return true;
    },
  });
}

const raf = DomHost.requestAnimationFrame;
beforeEach(() => {
  vi.spyOn(HTMLCanvasElement.prototype, 'getContext').mockImplementation(function (this: HTMLCanvasElement) {
    return ((this as unknown as { __ctx?: unknown }).__ctx ??= stubContext(this)) as never;
  });
  DomHost.requestAnimationFrame = (fn) => {
    setTimeout(() => fn(performance.now()), 16);
  };
  vi.useFakeTimers({ toFake: ['setTimeout', 'clearTimeout', 'Date', 'performance'] });
});

afterEach(() => {
  DomHost.requestAnimationFrame = raf;
  mounted?.remove();
  mounted = undefined;
  vi.useRealTimers();
});

describe('the values look story in a DOM', () => {
  it('mounts every cell, light and dark, without a throw', async () => {
    const el = (Looks.render as () => HTMLElement)();
    document.body.append(el);
    mounted = el;
    await vi.advanceTimersByTimeAsync(600);
    for (const theme of ['light', 'dark'])
      for (const cell of AllCells) {
        const host = el.querySelector(`[data-msui-theme="${theme}"] [data-cell="${cell.id}"]`);
        expect(host, `${theme} ${cell.id}`).not.toBeNull();
        expect(host!.querySelector('canvas'), `${theme} ${cell.id}`).not.toBeNull();
      }
  });
});
