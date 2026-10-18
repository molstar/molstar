/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import type { Theme } from '@molstar/graphics/theme/theme';
import { GeometryRequest } from '../request.js';

function theme(color: string, size: string) {
  return { color: { granularity: color }, size: { granularity: size } } as unknown as Theme;
}

describe('geometry request', () => {
  it('asks for positions where a color or size theme samples them', () => {
    expect(GeometryRequest.fromTheme(theme('group', 'uniform')).positions).toBe(false);
    expect(GeometryRequest.fromTheme(theme('volume', 'groupInstance')).positions).toBe(false);
    expect(GeometryRequest.fromTheme(theme('vertex', 'uniform')).positions).toBe(true);
    expect(GeometryRequest.fromTheme(theme('vertexInstance', 'group')).positions).toBe(true);
    expect(GeometryRequest.fromTheme(theme('uniform', 'vertex')).positions).toBe(true);
    expect(GeometryRequest.fromTheme(theme('instance', 'vertexInstance')).positions).toBe(true);
  });
});
