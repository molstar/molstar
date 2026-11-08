/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { Mat4, Vec3 } from '@molstar/core/math/linear-algebra';
import type { Clip } from '@molstar/core/util/clip';
import { CustomProperties } from '@molstar/model/model/custom-property';
import type { ParticleList, ParticleTarget } from '@molstar/model/model/particles/particle-list';
import { ColorTheme } from '@molstar/graphics/theme/color';
import { SizeTheme } from '@molstar/graphics/theme/size';
import { ParticleTargetRepresentationProvider } from '../target/representation.js';

// A clip object moved over particle targets (a clipping gizmo's drag): the clip alone set on what is drawn, merged
// objects included; a change of its defines or of another prop updates the targets in full.

/** a site at the origin, as a target */
const site = (): ParticleTarget => ({
  kind: 'sites',
  sites: { count: 1, positions: new Float32Array([0, 0, 0]), radii: new Float32Array([1]) },
});

/** four particles of two targets (two sphere visuals: merged) */
function particles(): ParticleList {
  return {
    count: 4,
    keys: new Int32Array([0, 1, 2, 3]),
    targets: new Int32Array([0, 0, 1, 1]),
    targetInfo: new Map([
      [0, {}],
      [1, {}],
    ]),
    targetMapping: new Map([
      [0, site()],
      [1, site()],
    ]),
    coordinates: new Float32Array([0, 0, 0, 10, 0, 0, 0, 10, 0, 10, 10, 0]),
    radii: new Float32Array([1, 1, 1, 1]),
    getParticleLabel: (index) => `${index}`,
    sourceData: { kind: 'test', name: 'test', data: {} },
    customProperties: new CustomProperties(),
    _propertyData: Object.create(null),
  };
}

const plane = (x: number) => ({
  type: 'plane' as const,
  invert: false,
  position: Vec3.create(x, 0, 0),
  rotation: { axis: Vec3.create(0, 0, 1), angle: 90 },
  scale: Vec3.create(1, 1, 1),
  transform: Mat4.identity(),
});
const clipOf = (...xs: number[]): Clip.Props => ({ variant: 'pixel', objects: xs.map(plane) });

/** the clip objects' x positions and count of what is drawn */
const drawnClip = (repr: ReturnType<typeof create>) =>
  repr.renderObjects.map((ro) => {
    const v = ro.values as unknown as { clipObjects: { ref: { value: Clip.Objects } } };
    const objects = v.clipObjects.ref.value;
    return Array.from({ length: objects.count }, (_, i) => objects.position[i * 3]);
  });

/** the version of the drawn object's `uLod`: a full update sets it (`BaseGeometry.updateValues`), the clip alone not */
const lodVersion = (repr: ReturnType<typeof create>) =>
  (repr.renderObjects[0].values as unknown as { uLod: { ref: { version: number } } }).uLod.ref.version;

function create() {
  const provider = ParticleTargetRepresentationProvider;
  return provider.factory(
    { colorThemeRegistry: ColorTheme.createRegistry(), sizeThemeRegistry: SizeTheme.createRegistry() },
    provider.getParams,
  );
}

describe('particle targets, a clip object moved', () => {
  it('sets the moved clip on what is drawn, the merged object kept', async () => {
    const list = particles();
    const repr = create();
    await repr.createOrUpdate({ renderMerged: true, gpuCulling: false, clip: clipOf(1) }, list).run();
    expect(repr.renderObjects).toHaveLength(1);
    const merged = repr.renderObjects[0];
    expect(drawnClip(repr)).toEqual([[1]]);
    const lod = lodVersion(repr);

    for (const x of [2, 3, 4]) {
      await repr.createOrUpdate({ clip: clipOf(x) }).run();
      expect(repr.renderObjects).toEqual([merged]);
      expect(drawnClip(repr)).toEqual([[x]]);
      // the targets not updated in full
      expect(lodVersion(repr)).toBe(lod);
    }
    repr.destroy();
  });

  it('updates the targets in full where the clip changes its defines or another prop changes with it', async () => {
    const list = particles();
    const repr = create();
    await repr.createOrUpdate({ renderMerged: true, gpuCulling: false, clip: clipOf(1) }, list).run();

    // a second clip object: another define, a full update
    const lod = lodVersion(repr);
    await repr.createOrUpdate({ clip: clipOf(1, 5) }).run();
    expect(drawnClip(repr)).toEqual([[1, 5]]);
    expect(lodVersion(repr)).not.toBe(lod);

    // another prop with the clip: both applied
    await repr.createOrUpdate({ alpha: 0.5, clip: clipOf(2, 5) }).run();
    expect(drawnClip(repr)).toEqual([[2, 5]]);
    expect((repr.renderObjects[0].values as unknown as { alpha: { ref: { value: number } } }).alpha.ref.value).toBe(
      0.5,
    );
    repr.destroy();
  });
});
