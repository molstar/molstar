/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { createClipping, createEmptyClipping } from '../clipping-data.js';
import { createEmissive, createEmptyEmissive } from '../emissive-data.js';
import { createEmptyMarkers, createMarkers } from '../marker-data.js';
import { createEmptyOverpaint, createOverpaint } from '../overpaint-data.js';
import { createEmptySubstance, createSubstance } from '../substance-data.js';
import { createEmptyTransparency, createTransparency } from '../transparency-data.js';
import { createEmptyWiggle, createWiggle } from '../wiggle-data.js';

describe('layer data', () => {
  // every object's empty layer the same array: a small layer made from it is its own, not written into it
  const layers = [
    ['overpaint', () => createEmptyOverpaint(), (d: any) => createOverpaint(1, 'instance', d), 'tOverpaint'],
    [
      'transparency',
      () => createEmptyTransparency(),
      (d: any) => createTransparency(1, 'instance', d),
      'tTransparency',
    ],
    ['clipping', () => createEmptyClipping(), (d: any) => createClipping(1, 'instance', d), 'tClipping'],
    ['emissive', () => createEmptyEmissive(), (d: any) => createEmissive(1, 'instance', d), 'tEmissive'],
    ['substance', () => createEmptySubstance(), (d: any) => createSubstance(1, 'instance', d), 'tSubstance'],
    ['wiggle', () => createEmptyWiggle(), (d: any) => createWiggle(1, 'instance', d), 'tWiggle'],
    ['markers', () => createEmptyMarkers(), (d: any) => createMarkers(1, 'instance', d), 'tMarker'],
  ] as const;

  for (const [name, empty, create, key] of layers) {
    it(`makes ${name} of two objects from their empty layers each its own array`, () => {
      const a = create(empty()) as any;
      const b = empty() as any;
      a[key].ref.value.array[0] = 7;
      // the second's empty array untouched, and its own layer cleared, not the first's
      expect(b[key].ref.value.array[0]).toBe(0);
      create(b);
      expect(b[key].ref.value.array).not.toBe(a[key].ref.value.array);
      expect(a[key].ref.value.array[0]).toBe(7);
    });
  }
});
