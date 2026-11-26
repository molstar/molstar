/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { calculateBoundingSphere } from '../renderable/util.js';
import { createDrawVersion, createDrawVersionCounter, DrawKey } from '../renderable.js';

describe('renderable', () => {
  it('calculateBoundingSphere', () => {
    const position = new Float32Array([0, 0, 0, 1, 0, 0, -1, 0, 0]);
    const transform = new Float32Array([
      1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,

      1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0,

      1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, -1, 0, 0, 0,
    ]);

    const { boundingSphere, invariantBoundingSphere } = calculateBoundingSphere(
      position,
      position.length / 3,
      transform,
      transform.length / 16,
    );

    expect(invariantBoundingSphere.extrema).toEqual([
      [0, 0, 0],
      [1, 0, 0],
      [-1, 0, 0],
    ]);
    expect(invariantBoundingSphere.radius).toBe(1);
    expect(invariantBoundingSphere.center).toEqual([0, 0, 0]);
    expect(boundingSphere.radius).toBe(2);
    expect(boundingSphere.center).toEqual([0, 0, 0]);
  });
});

describe('createDrawVersion', () => {
  it('keeps a version per key on a shared counter: two keys read in turn advance neither, and never meet', () => {
    const counter = createDrawVersionCounter();
    let a = 1;
    let b = 2;
    const getA = createDrawVersion((out) => DrawKey.add(out, a), counter);
    const getB = createDrawVersion((out) => DrawKey.add(out, b), counter);
    const va = getA();
    const vb = getB();
    expect(va).not.toBe(vb);
    for (let i = 0; i < 3; ++i) {
      expect(getA()).toBe(va);
      expect(getB()).toBe(vb);
    }
    a = 3;
    const va2 = getA();
    expect(va2).not.toBe(va);
    expect(va2).not.toBe(vb);
    b = 4;
    expect(getB()).not.toBe(va2);
  });
});
