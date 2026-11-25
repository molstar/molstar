/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Tadej Satler <tadej.satler@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { createPreIntegratedTransferTexture, createTransferFunctionTexture } from '../transfer-function.js';

/** A transfer function of one narrow peak: opacity 0.8 at value 128 only */
function peak() {
  const alpha = new Uint8Array(256);
  alpha[128] = Math.round(0.8 * 255);
  return alpha;
}

describe('pre-integrated transfer function', () => {
  const alpha = peak();
  const { array, width, height } = createPreIntegratedTransferTexture(alpha).ref.value;
  const at = (front: number, back: number) => array[back * width + front];

  it('is a table of 256 by 256, the 1D opacity on its diagonal', () => {
    expect([width, height]).toEqual([256, 256]);
    for (let i = 0; i < 256; ++i) expect(at(i, i)).toBe(alpha[i]);
  });

  it('counts a peak crossed between two samples, which neither sample lands on', () => {
    expect(at(120, 136)).toBeGreaterThan(0);
    expect(at(136, 120)).toBe(at(120, 136));
    // the longer the step in value, the less of it the peak is
    expect(at(100, 156)).toBeLessThan(at(120, 136));
    // a step that does not reach the peak
    expect(at(10, 100)).toBe(0);
  });

  it('is updated in place', () => {
    const cell = createPreIntegratedTransferTexture(alpha);
    const before = cell.ref.value.array;
    createPreIntegratedTransferTexture(new Uint8Array(256), cell);
    expect(cell.ref.value.array).toBe(before);
    expect(cell.ref.value.array.every((v) => v === 0)).toBe(true);
  });
});

describe('transfer function', () => {
  it('stays opaque where the control points are, the spline overshooting 1 between them', () => {
    // a ramp to 1 at 0.5, held at 1 to the end
    const points = [
      { x: 0.2, alpha: 0 },
      { x: 0.5, alpha: 1 },
      { x: 1, alpha: 1 },
    ];
    const { array } = createTransferFunctionTexture(points).ref.value;
    expect(Math.min(...array.slice(128))).toBeGreaterThan(240);
  });
});
