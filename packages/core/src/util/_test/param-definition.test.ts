/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it, vi } from 'vitest';
import { ParamDefinition as PD } from '../param-definition.js';

describe('PD.isParamEqual', () => {
  it('compares array and typed-array values by index, not by keys (perf-audit-2 01/P8)', () => {
    const p = PD.Value<any>(undefined);
    const keys = vi.spyOn(Object, 'keys');
    expect(PD.isParamEqual(p, [1, 2, 3], [1, 2, 3])).toBe(true);
    expect(PD.isParamEqual(p, [1, 2, 3], [1, 2, 4])).toBe(false);
    expect(PD.isParamEqual(p, [1, 2], [1, 2, 3])).toBe(false);
    expect(PD.isParamEqual(p, Float32Array.of(1, 2), Float32Array.of(1, 2))).toBe(true);
    expect(PD.isParamEqual(p, Float32Array.of(1, 2), Float32Array.of(1, 3))).toBe(false);
    expect(keys).not.toHaveBeenCalled();
    keys.mockRestore();
    // objects as before
    expect(PD.isParamEqual(p, { a: 1 }, { a: 1 })).toBe(true);
    expect(PD.isParamEqual(p, { a: 1 }, { a: 2 })).toBe(false);
  });
});
