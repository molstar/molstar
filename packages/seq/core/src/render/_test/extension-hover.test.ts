/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import type { AxisTransform } from '../axis.js';
import { Round } from '../draw-lists.js';
import type { HitMark } from '../layer-type.js';
import { LightTokens } from '../tokens.js';
import { drawMarkOutline } from '../tracks-emphasis.js';
import { DefaultRenderConfig } from '../config.js';

// The outline of an extension's hovered mark (_spec/README.md "Render", Emphasis): outside the mark by the hover's
// offset and concentric, round on the sides the strip does not cut, square where the mark meets its edge.

const axis = { x0: 10, width: 100 } as AxisTransform;

function outline(mark: HitMark) {
  const rects: { args: number[] }[] = [];
  const args = new Float64Array(16);
  const out = { args, rectArgs: () => rects.push({ args: [...args] }) };
  drawMarkOutline(out as never, mark, axis, LightTokens, DefaultRenderConfig);
  return rects[0].args;
}

describe('an extension mark’s hover outline', () => {
  const half = DefaultRenderConfig.outline.hover / 2 + DefaultRenderConfig.emphasis.hoverOffset;

  it('surrounds the mark by the offset, concentric with its radius', () => {
    const g = outline({ x: 40, y: 5, width: 8, height: 8, radius: 4 });
    expect(g.slice(0, 4)).toEqual([40 - half, 5 - half, 8 + 2 * half, 8 + 2 * half]);
    expect(g[7]).toBe(4 + half);
    expect(g[8]).toBe(Round.left | Round.right);
  });

  it('is square on the side the strip cuts', () => {
    expect(outline({ x: 10, y: 0, width: 8, height: 8, radius: 4 })[8]).toBe(Round.right);
    expect(outline({ x: 102, y: 0, width: 8, height: 8, radius: 4 })[8]).toBe(Round.left);
    expect(outline({ x: 10, y: 0, width: 100, height: 8 })[8]).toBe(0);
  });
});
