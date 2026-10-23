/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { Mat4, Tensor, Vec3 } from '@molstar/core/math/linear-algebra';
import { RuntimeContext } from '@molstar/core/task';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { CustomProperties } from '@molstar/model/model/custom-property';
import type { ParticleList } from '@molstar/model/model/particles/particle-list';
import { type Grid, Volume } from '@molstar/model/model/volume';
import { ColorTheme } from '@molstar/graphics/theme/color';
import { BuiltInColorThemes } from '@molstar/graphics/theme/color/catalog';
import { SizeTheme } from '@molstar/graphics/theme/size';
import { BuiltInSizeThemes } from '@molstar/graphics/theme/size/catalog';
import { Theme } from '@molstar/graphics/theme/theme';
import type { GraphicsRenderObject } from '@molstar/graphics/gl/render-object';
import { SpacefillParticlesRepresentationProvider } from '../particles/representation/spacefill.js';
import { IsosurfaceRepresentationProvider } from '../volume/isosurface.js';

const registries = (() => {
  const colorThemeRegistry = ColorTheme.createRegistry();
  for (const p of Object.values(BuiltInColorThemes)) colorThemeRegistry.add(p as ColorTheme.Provider<any, any, any>);
  const sizeThemeRegistry = SizeTheme.createRegistry();
  for (const p of Object.values(BuiltInSizeThemes)) sizeThemeRegistry.add(p as SizeTheme.Provider<any, any>);
  return { colorThemeRegistry, sizeThemeRegistry };
})();

const N = 10;

/** A sphere of radius `r` in an N³ grid: a different `r` is data that is not equivalent (a new render object) */
function createVolume(r: number): Volume {
  const space = Tensor.Space([N, N, N], [2, 1, 0], Float32Array);
  const data = space.create();
  for (let k = 0; k < N; k++) {
    for (let j = 0; j < N; j++) {
      for (let i = 0; i < N; i++) space.set(data, i, j, k, r - Math.hypot(i - 4.5, j - 4.5, k - 4.5));
    }
  }
  return {
    grid: {
      transform: { kind: 'matrix', matrix: Mat4.identity() },
      cells: Tensor.create(space, data),
      stats: { min: r - 8, max: r, mean: 0, sigma: 1 },
    } satisfies Grid,
    instances: [{ transform: Mat4.identity() }],
    sourceData: { kind: 'test', name: 'test', data: {} } as any,
    customProperties: new CustomProperties(),
    _propertyData: Object.create(null),
    _localPropertyData: Object.create(null),
  };
}

function createParticles(): ParticleList {
  return {
    count: 3,
    keys: new Int32Array([0, 1, 2]),
    targets: new Int32Array(3),
    targetInfo: new Map([[0, {}]]),
    coordinates: new Float32Array([0, 0, 0, 4, 0, 0, 8, 0, 0]),
    getParticleLabel: (index) => `${index}`,
    sourceData: { kind: 'test', name: 'test', data: {} },
    customProperties: new CustomProperties(),
    _propertyData: Object.create(null),
  };
}

const transform = Mat4.fromTranslation(Mat4(), Vec3.create(100, 0, 0));

function expectState(renderObjects: ReadonlyArray<GraphicsRenderObject>) {
  expect(renderObjects.length).toBeGreaterThan(0);
  for (const ro of renderObjects) {
    expect(ro.state.visible).toBe(false);
    expect(ro.state.alphaFactor).toBe(0.25);
    expect(Mat4.areEqual(ro.values.matrix.ref.value as Mat4, transform, 1e-6)).toBe(true);
  }
}

describe('new render objects get the representation state (perf-audit-2 07/B1)', () => {
  it('volume: data that is not equivalent keeps visibility, alpha and transform', async () => {
    const provider = IsosurfaceRepresentationProvider;
    const v1 = createVolume(3);
    const repr = provider.factory({ ...registries }, provider.getParams);
    repr.setTheme(
      Theme.create(registries, { volume: v1 }, { colorTheme: { name: 'uniform' }, sizeTheme: { name: 'uniform' } }),
    );
    const params = {
      ...PD.getDefaultValues(provider.getParams(registries, v1)),
      isoValue: Volume.IsoValue.absolute(0),
    };
    await repr.createOrUpdate(params, v1).runInContext(RuntimeContext.Synchronous);
    const before = repr.renderObjects[0];

    repr.setState({ visible: false, alphaFactor: 0.25, transform });
    // a streamed box: new data, a new render object
    await repr.createOrUpdate(params, createVolume(2)).runInContext(RuntimeContext.Synchronous);

    expect(repr.renderObjects[0]).not.toBe(before);
    expectState(repr.renderObjects);
    repr.destroy();
  });

  it('particles: a visual made after the state was set gets it', async () => {
    const provider = SpacefillParticlesRepresentationProvider;
    const particles = createParticles();
    const repr = provider.factory({ ...registries }, provider.getParams);
    repr.setTheme(Theme.create(registries, {}, { colorTheme: { name: 'uniform' }, sizeTheme: { name: 'uniform' } }));
    repr.setState({ visible: false, alphaFactor: 0.25, transform });

    const params = PD.getDefaultValues(provider.getParams(registries, particles));
    await repr.createOrUpdate(params, particles).runInContext(RuntimeContext.Synchronous);

    expectState(repr.renderObjects);
    repr.destroy();
  });
});
