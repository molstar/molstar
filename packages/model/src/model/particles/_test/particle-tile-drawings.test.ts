/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it, vi } from 'vitest';
import { Vec3 } from '@molstar/core/math/linear-algebra/3d/vec3';
import { Clip } from '@molstar/core/util/clip';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { CustomProperties } from '../../custom-property.js';
import type { ParticleList } from '../particle-list.js';
import { ParticleTileDrawings } from '../particle-tile-drawings.js';
import { ParticleTileView } from '../particle-tile-view.js';
import { ParticleTiles } from '../particle-tiles.js';

/** 8 particles on a 2³ lattice, a target of one site */
function createList(): ParticleList {
  const n = 8;
  return {
    label: 'lattice',
    count: n,
    targets: new Int32Array(n).fill(1),
    targetInfo: new Map([[1, { name: 'one' }]]),
    targetMapping: new Map([
      [1, { kind: 'sites', sites: { count: 1, positions: new Float32Array(3), radii: new Float32Array(1).fill(1) } }],
    ]),
    coordinates: Float32Array.from({ length: n * 3 }, (_, k) => ((Math.floor(k / 3) >> (k % 3)) & 1) * 30),
    getParticleLabel: (i) => `particle ${i}`,
    sourceData: { kind: 'test', name: 'test', data: {} },
    customProperties: new CustomProperties(),
    _propertyData: Object.create(null),
  };
}

const sphereAt = (x: number) => ({
  variant: 'pixel' as const,
  objects: [
    { ...PD.getDefaultValues(Clip.Params.objects.element), type: 'sphere' as const, position: Vec3.create(x, 0, 0) },
  ],
});

describe('particle tile drawings', () => {
  it("keeps what draws a set's views, tells of a change (not of the same again), lets go of what is gone", () => {
    const set = ParticleTiles.fromList(createList(), [2, 2, 2]);
    const source = {
      getTile: (id: number, frame: number) => set.provider.getTile(id, frame) as ParticleList,
      getTarget: (id: number) => set.provider.getTarget(id) as never,
    };
    const view = ParticleTileView.create(set, source, [...set.tiles.keys()], 0).list;
    const told = vi.fn();
    const unsubscribe = ParticleTileDrawings.subscribe(set, told);
    const owner = {};
    expect(ParticleTileDrawings.get(set).list).toEqual([]);

    ParticleTileDrawings.set(owner, view, { clip: sphereAt(0) });
    expect(told).toHaveBeenCalledTimes(1);
    const first = ParticleTileDrawings.get(set);
    expect(first.list).toHaveLength(1);
    expect(first.list[0].clip?.count).toBe(1);
    expect(first.list[0].opaque).toBe(true);
    expect(ParticleTileDrawings.isChanging(set)).toBe(true);
    expect(ParticleTileDrawings.isChanging(set, performance.now() + ParticleTileDrawings.SettleMs)).toBe(false);

    // the same again: not told, the version kept
    const version = first.version;
    ParticleTileDrawings.set(owner, view, { clip: sphereAt(0) });
    expect(told).toHaveBeenCalledTimes(1);
    expect(ParticleTileDrawings.get(set).version).toBe(version);

    // moved, or drawn not opaque: told
    ParticleTileDrawings.set(owner, view, { clip: sphereAt(10) });
    ParticleTileDrawings.set(owner, view, { clip: sphereAt(10), alpha: 0.5 });
    expect(told).toHaveBeenCalledTimes(3);
    expect(ParticleTileDrawings.get(set).list[0].opaque).toBe(false);

    // a list not of tiles: drawing none of the set, let go of
    ParticleTileDrawings.set(owner, createList(), { clip: sphereAt(10) });
    expect(ParticleTileDrawings.get(set).list).toEqual([]);
    expect(told).toHaveBeenCalledTimes(4);
    ParticleTileDrawings.set(owner, view, {});
    ParticleTileDrawings.set(owner, undefined, undefined);
    expect(ParticleTileDrawings.get(set).list).toEqual([]);
    unsubscribe();
    ParticleTileDrawings.set(owner, view, {});
    expect(told).toHaveBeenCalledTimes(6);
  });
});
