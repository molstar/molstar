/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import type { Location } from '@molstar/model/model/location';
import { LocationIterator } from '../location-iterator.js';

/** A location per (group, instance), `-1` for a voided instance */
const at = (groupIndex: number, instanceIndex: number) =>
  ({ kind: 'test', groupIndex, instanceIndex }) as unknown as Location;

describe('location iterator', () => {
  it("gives an instance's first location by index as moving to it does, also with instances voided", () => {
    for (const voided of [false, true]) {
      const iterator = LocationIterator(3, 4, 1, at, false, (g, i) => (g + i) % 2 === 1);
      if (voided) iterator.voidInstances();
      const moved: [unknown, boolean][] = [];
      while (iterator.hasNext) {
        const { location, isSecondary } = iterator.move();
        moved.push([{ ...location }, isSecondary]);
        iterator.skipInstance();
      }
      const byIndex = Array.from({ length: 4 }, (_, i) => [
        { ...iterator.getInstanceLocation(i) },
        iterator.isInstanceSecondary(i),
      ]);
      expect(byIndex).toEqual(moved);
      expect(moved[1][0]).toEqual({ kind: 'test', groupIndex: 0, instanceIndex: voided ? -1 : 1 });
    }
  });
});
