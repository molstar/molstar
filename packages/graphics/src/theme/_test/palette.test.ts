/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Color } from '@molstar/core/util/color';
import { describe, expect, it } from 'vitest';
import { ColorTheme } from '../color.js';

// A palette from a color list with offsets (the volume-value theme's): interpolated between the stops at their offsets,
// as `ColorScale` does; a set keeps each stop's color up to the middle between stops.

const at = (p: ColorTheme.Palette, t: number) => p.colors[Math.round(t * (p.colors.length - 1))];

describe('ColorTheme.Palette with offsets', () => {
  it('interpolates between the stops, their colors before the first and after the last', () => {
    const p = ColorTheme.Palette(
      [
        [Color(0x000000), 0.2],
        [Color(0xff0000), 0.6],
      ],
      'interpolate',
    );
    expect(p.filter).toBe('linear');
    expect(at(p, 0)).toBe(0x000000);
    expect(at(p, 0.2)).toBe(0x000000);
    const mid = Color.toRgb(at(p, 0.4));
    expect(mid[0]).toBeGreaterThan(120);
    expect(mid[0]).toBeLessThan(136);
    expect(at(p, 0.6)).toBe(0xff0000);
    expect(at(p, 1)).toBe(0xff0000);
  });

  it('keeps a set as steps, and offsets with many decimals do not make a huge palette', () => {
    const p = ColorTheme.Palette(
      [
        [Color(0x0000ff), 0],
        [Color(0x00ff00), 0.123456789],
        [Color(0xff0000), 1],
      ],
      'set',
    );
    expect(p.colors.length).toBeLessThanOrEqual(4096);
    expect(at(p, 0.03)).toBe(0x0000ff);
    expect(at(p, 0.2)).toBe(0x00ff00);
    expect(at(p, 0.9)).toBe(0xff0000);
  });

  it('gives a set with integer offsets, and the part after its last stop, colors (perf-audit-2 09/B3)', () => {
    const ends = ColorTheme.Palette(
      [
        [Color(0x0000ff), 0],
        [Color(0xff0000), 1],
      ],
      'set',
    );
    expect(ends.colors).toEqual([0x0000ff, 0xff0000]);

    const inside = ColorTheme.Palette(
      [
        [Color(0x0000ff), 0.2],
        [Color(0xff0000), 0.6],
      ],
      'set',
    );
    expect(inside.colors.length).toBe(10);
    expect(inside.colors.every((c) => c !== undefined)).toBe(true);
    expect(at(inside, 0)).toBe(0x0000ff);
    expect(at(inside, 1)).toBe(0xff0000);
  });
});
