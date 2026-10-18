/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// VL-21, LC-06: every renderer folds a color's alpha into its key by one rule (`packAlpha`): an opaque color is its own
// key, and no alpha, 0 included, turns a color into NO_COLOR or a skin token's key (which a near-white color with 255
// in its top byte did).

import { expect, it } from 'vitest';
import { NO_COLOR, packAlpha, SkinToken } from '../src/draw/api.js';

it('a key is the color when opaque, and never NO_COLOR or a skin token', () => {
  const tokens = new Set<number>([NO_COLOR, ...Object.values(SkinToken)]);
  for (const c of [0, 0x123456, 0xfffffe, 0xffffff - 3, 0xffffff])
    for (let a = 0; a <= 255; a++) {
      const k = packAlpha(c, a);
      expect(tokens.has(k), `${c.toString(16)} at ${a}`).toBe(false);
      if (a === 255) expect(k).toBe(c);
    }
});
