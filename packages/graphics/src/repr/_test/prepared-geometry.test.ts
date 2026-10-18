/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import type { Geometry, GeometryUtils } from '@molstar/graphics/geo/geometry/geometry';
import type { GeometryRequest } from '@molstar/graphics/geo/geometry/request';
import type { Theme } from '@molstar/graphics/theme/theme';
import { Visual } from '../visual.js';

const geometry = { kind: 'mesh' } as unknown as Geometry;
const vertexColors = { color: { granularity: 'vertex' }, size: { granularity: 'uniform' } } as unknown as Theme;
const groupColors = { color: { granularity: 'group' }, size: { granularity: 'uniform' } } as unknown as Theme;

describe('Visual.withPreparedGeometry', () => {
  it('updates synchronously without a preparation step', () => {
    const updated: Geometry[] = [];
    const result = Visual.withPreparedGeometry(geometry, {} as GeometryUtils<Geometry>, vertexColors, undefined, (g) =>
      updated.push(g),
    );
    expect(result).toBeUndefined();
    expect(updated).toEqual([geometry]);
  });

  it('updates synchronously when the preparation has nothing to do', () => {
    const requests: GeometryRequest[] = [];
    const utils = { prepare: (_: Geometry, r: GeometryRequest) => void requests.push(r) } as GeometryUtils<Geometry>;
    let updated = false;
    const result = Visual.withPreparedGeometry(geometry, utils, groupColors, undefined, () => (updated = true));
    expect(result).toBeUndefined();
    expect(updated).toBe(true);
    expect(requests).toEqual([{ positions: false }]);
  });

  it('updates once the preparation and an asynchronous creation are done', async () => {
    const order: string[] = [];
    let resolvePrepare!: () => void;
    const utils = {
      prepare: (_: Geometry, r: GeometryRequest) => {
        order.push(`prepare ${r.positions}`);
        return new Promise<void>((resolve) => (resolvePrepare = resolve));
      },
    } as GeometryUtils<Geometry>;
    const result = Visual.withPreparedGeometry(Promise.resolve(geometry), utils, vertexColors, undefined, () =>
      order.push('update'),
    );
    await Promise.resolve();
    await Promise.resolve();
    expect(order).toEqual(['prepare true']);
    resolvePrepare();
    await result;
    expect(order).toEqual(['prepare true', 'update']);
  });
});
