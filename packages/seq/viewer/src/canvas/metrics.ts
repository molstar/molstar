/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { TextMetricsPort } from '@molstar/seqx-core/context/index';
import { DefaultFonts, type FontStacks } from '@molstar/seqx-core/render';
import { FontCache } from './font.js';

// The text metrics port's browser adapter (_spec/README.md "Canvas2D backend"): widths measured with a canvas context's
// `measureText`, in the same CSS font the backend draws with.

/** What the adapter reads of a 2D context. */
export interface MeasuringContext {
  font: string;
  measureText(text: string): { readonly width: number };
}

/** A metrics port over a context (a canvas the adapter keeps off screen); the font set once per change. */
export function canvasTextMetrics(ctx: MeasuringContext, stacks: FontStacks = DefaultFonts): TextMetricsPort {
  const fonts = new FontCache(stacks);
  let current = '';
  return {
    measure(font, text) {
      const css = fonts.css(font);
      if (css !== current) {
        ctx.font = css;
        current = css;
      }
      return ctx.measureText(text).width;
    },
  };
}

/** The browser's adapter over a new off-screen canvas; `undefined` where there is no DOM (Node: the bundled table). */
export function browserTextMetrics(stacks: FontStacks = DefaultFonts): TextMetricsPort | undefined {
  if (typeof document === 'undefined') return undefined;
  const ctx = document.createElement('canvas').getContext('2d');
  return ctx ? canvasTextMetrics(ctx, stacks) : undefined;
}
