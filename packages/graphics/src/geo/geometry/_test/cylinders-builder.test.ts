/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { CylindersBuilder } from '../cylinders/cylinders-builder.js';
import { CylinderStride } from '../cylinders/cylinders.js';

/** cylinder i: start (i, 0, 0), end (i, 1, 0), scale i / 10, caps by i % 4, colorMode i % 3, group 100 + i */
function build(count: number, builder: ReturnType<typeof CylindersBuilder.create>) {
  for (let i = 0; i < count; i++) builder.add(i, 0, 0, i, 1, 0, i / 10, (i & 1) > 0, (i & 2) > 0, i % 3, 100 + i);
  return builder.getCylinders();
}

/** the data of each cylinder: start and group, end and scale, cap and color mode */
function expectCylinders(cylinders: ReturnType<typeof build>, count: number) {
  expect(cylinders.cylinderCount).toBe(count);
  const data = cylinders.cylinderData.ref.value;
  expect(data.count).toBe(count * 3);
  for (let i = 0; i < count; i++) {
    const o = i * CylinderStride;
    expect(Array.from(data.array.subarray(o, o + 4))).toEqual([i, 0, 0, 100 + i]);
    expect(Array.from(data.array.subarray(o + 4, o + 8))).toEqual([i, 1, 0, Math.fround(i / 10)]);
    expect(data.array[o + 8]).toBe((i & 1 ? 1 : 0) + (i & 2 ? 2 : 0));
    expect(data.array[o + 9]).toBe(i % 3);
  }
}

describe('CylindersBuilder.add', () => {
  it('writes the data of a cylinder once, not per vertex', () => {
    expectCylinders(build(100, CylindersBuilder.create(100, 50)), 100);
  });

  it('crosses chunk ends, and refills reused buffers', () => {
    // chunks of 2 cylinders after a first of 1
    const cylinders = build(7, CylindersBuilder.create(1, 2));
    expectCylinders(cylinders, 7);
    expectCylinders(build(5, CylindersBuilder.create(5, 2, cylinders)), 5);
  });
});
