/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { Vec4 } from '@molstar/core/math/linear-algebra';
import { flipWindowRect } from '../window-rect.js';

describe('window-rect', () => {
  it('flips to a bottom-left origin', () => {
    const out = Vec4();
    expect(Array.from(flipWindowRect(out, Vec4.create(0.1, 0.2, 0.5, 0.6), 'bounds', 0))).toEqual([0.1, 0.4, 0.5, 0.8]);
    expect(Array.from(flipWindowRect(out, Vec4.create(0.1, 0.2, 0.3, 0.4), 'rect', 0))).toEqual([0.1, 0.4, 0.3, 0.4]);
    expect(Array.from(flipWindowRect(out, Vec4.create(5, 10, 20, 30), 'pixels', 100))).toEqual([5, 60, 20, 30]);
  });
});
