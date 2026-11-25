/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The missing colour (_spec/README.md "Render", The matrix layer, Color): a missing value is told from the ramp's low
// end and from the page by at least 2:1, in both themes (WCAG's relative luminance contrast, from the tokens).

import { describe, expect, it } from 'vitest';
import { DarkTokens, LightTokens } from '../tokens.js';

const lin = (v: number) => (v / 255 <= 0.04045 ? v / 255 / 12.92 : ((v / 255 + 0.055) / 1.055) ** 2.4);
/** A token's (0xRRGGBBAA) relative luminance. */
function luminance(c: number): number {
  return 0.2126 * lin((c >>> 24) & 255) + 0.7152 * lin((c >>> 16) & 255) + 0.0722 * lin((c >>> 8) & 255);
}
function contrast(a: number, b: number): number {
  const [x, y] = [luminance(a), luminance(b)].sort((p, q) => q - p);
  return (x + 0.05) / (y + 0.05);
}

describe('the missing colour', () => {
  for (const tokens of [LightTokens, DarkTokens])
    it(`${tokens.theme}: at least 2:1 against the ramp’s low end and the page`, () => {
      const missing = tokens.colors['seq-missing'];
      expect(contrast(missing, tokens.colors['seq-ramp-0'])).toBeGreaterThanOrEqual(2);
      expect(contrast(missing, tokens.colors['msui-bg'])).toBeGreaterThanOrEqual(2);
    });
});
