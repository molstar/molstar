/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Ludovic Autin <ludovic.autin@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import fs from 'node:fs';
import path from 'node:path';
import { describe, expect, it } from 'vitest';
import { Mat4, Tensor, Vec3 } from '@molstar/core/math/linear-algebra';
import { RuntimeContext } from '@molstar/core/task';
import { DefaultAssetManager } from '@molstar/core/util/assets';
import type { Clip } from '@molstar/core/util/clip';
import { Color } from '@molstar/core/util/color';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { CIF } from '@molstar/io/reader/cif';
import { trajectoryFromMmCIF } from '@molstar/model/formats/structure/mmcif';
import { CustomProperties } from '@molstar/model/model/custom-property';
import { Shape } from '@molstar/model/model/shape';
import { Structure } from '@molstar/model/model/structure';
import { type Grid, Volume } from '@molstar/model/model/volume';
import type { CustomProperty } from '@molstar/model/props/common/custom-property';
import { Mesh } from '@molstar/graphics/geo/geometry/mesh/mesh';
import { MeshBuilder } from '@molstar/graphics/geo/geometry/mesh/mesh-builder';
import { addSphere } from '@molstar/graphics/geo/geometry/mesh/builder/sphere';
import type { GraphicsRenderObject } from '@molstar/graphics/gl/render-object';
import { ColorTheme } from '@molstar/graphics/theme/color';
import { BuiltInColorThemes } from '@molstar/graphics/theme/color/catalog';
import { SizeTheme } from '@molstar/graphics/theme/size';
import { BuiltInSizeThemes } from '@molstar/graphics/theme/size/catalog';
import { Theme } from '@molstar/graphics/theme/theme';
import type { Representation } from '../representation.js';
import { ShapeRepresentation } from '../shape/representation.js';
import { BallAndStickRepresentationProvider } from '../structure/representation/ball-and-stick.js';
import { CartoonRepresentationProvider } from '../structure/representation/cartoon.js';
import { GaussianSurfaceRepresentationProvider } from '../structure/representation/gaussian-surface.js';
import { MolecularSurfaceRepresentationProvider } from '../structure/representation/molecular-surface.js';
import { IsosurfaceRepresentationProvider } from '../volume/isosurface.js';

// A change of `clip` alone (a clip manager's live drag: Prism's clipping gizmos) updates uniforms
// and the clip defines only: no render object is made anew and no geometry, color or size data changes. Measured at
// 0.01-0.12 ms an update on 1cbs (cartoon, ball-and-stick, molecular and gaussian surface), a volume and a shape.

const registries = (() => {
  const colorThemeRegistry = ColorTheme.createRegistry();
  for (const p of Object.values(BuiltInColorThemes)) colorThemeRegistry.add(p as ColorTheme.Provider<any, any, any>);
  const sizeThemeRegistry = SizeTheme.createRegistry();
  for (const p of Object.values(BuiltInSizeThemes)) sizeThemeRegistry.add(p as SizeTheme.Provider<any, any>);
  return { colorThemeRegistry, sizeThemeRegistry };
})();
const propCtx: CustomProperty.Context = {
  runtime: RuntimeContext.Synchronous,
  assetManager: new DefaultAssetManager(),
};

async function load1cbs() {
  const data = fs.readFileSync(
    path.resolve(import.meta.dirname, '../../../../../data/examples/1cbs_updated.cif'),
    'utf8',
  );
  const parsed = await CIF.parse(data).run();
  if (parsed.isError) throw new Error(parsed.message);
  const trajectory = await trajectoryFromMmCIF(parsed.result.blocks[0], parsed.result).run();
  return Structure.ofModel(trajectory.representative);
}

async function createRepr<D>(provider: any, data: D, ctxData: object) {
  await provider.ensureCustomProperties?.attach(propCtx, data);
  const themeProps = {
    colorTheme: { name: provider.defaultColorTheme.name, params: provider.defaultColorTheme.props ?? {} },
    sizeTheme: { name: provider.defaultSizeTheme.name, params: provider.defaultSizeTheme.props ?? {} },
  };
  await Theme.ensureDependencies(propCtx, registries, ctxData, themeProps);
  const repr = provider.factory({ ...registries }, provider.getParams) as Representation<D>;
  repr.setTheme(Theme.create(registries, ctxData, themeProps));
  const params = PD.getDefaultValues(provider.getParams(registries, data));
  await repr.createOrUpdate(params, data).run();
  return repr;
}

const N = 24;
function createTestVolume(): Volume {
  const space = Tensor.Space([N, N, N], [2, 1, 0], Float32Array);
  const data = space.create();
  for (let k = 0; k < N; k++) {
    for (let j = 0; j < N; j++) {
      for (let i = 0; i < N; i++) space.set(data, i, j, k, 10 - Math.hypot(i - 11.5, j - 11.5, k - 11.5));
    }
  }
  return {
    grid: {
      transform: { kind: 'matrix', matrix: Mat4.identity() },
      cells: Tensor.create(space, data),
      stats: { min: -10, max: 10, mean: 0, sigma: 3 },
    } satisfies Grid,
    instances: [{ transform: Mat4.identity() }],
    sourceData: { kind: 'test', name: 'test', data: {} } as any,
    customProperties: new CustomProperties(),
    _propertyData: Object.create(null),
    _localPropertyData: Object.create(null),
  };
}

function createTestShape() {
  const state = MeshBuilder.createState(1024, 1024);
  for (let g = 0; g < 4; g++) {
    state.currentGroup = g;
    addSphere(state, Vec3.create(g * 3, 0, 0), 1, 2);
  }
  return Shape.create(
    'spheres',
    {},
    MeshBuilder.getMesh(state),
    () => Color(0x3f7fbf),
    () => 1,
    () => '',
    4,
  );
}

function plane(y: number): Clip.Props {
  return {
    variant: 'pixel',
    objects: [
      {
        type: 'plane',
        invert: false,
        position: Vec3.create(0, y, 0),
        rotation: { axis: Vec3.create(1, 0, 0), angle: 30 },
        scale: Vec3.create(1, 1, 1),
        transform: Mat4.identity(),
      },
    ],
  };
}

/** Every value cell of every render object, by its version */
function versions(ros: readonly GraphicsRenderObject[]) {
  return ros.map((ro) => {
    const v = new Map<string, number>();
    for (const [k, cell] of Object.entries(ro.values)) if (cell && 'ref' in cell) v.set(k, (cell as any).ref.version);
    return { ro, v };
  });
}

function changedValues(before: ReturnType<typeof versions>, after: readonly GraphicsRenderObject[]) {
  const changed = new Set<string>();
  after.forEach((ro, i) => {
    for (const [k, cell] of Object.entries(ro.values)) {
      if (cell && 'ref' in cell && before[i].v.get(k) !== (cell as any).ref.version) changed.add(k);
    }
  });
  return [...changed].sort();
}

/** What a clip-only update may touch: uniforms (some are set on every update, `uLod`), the clip defines and objects */
const ClipValues = /^(u\w+|dClipObjects|dClipObjectCount|dClipVariant|clipObjects)$/;

describe('a clip-only update', async () => {
  const structure = await load1cbs();
  const cases: [string, () => Promise<Representation<any>>][] = [
    ['cartoon', () => createRepr(CartoonRepresentationProvider, structure, { structure })],
    ['ball-and-stick', () => createRepr(BallAndStickRepresentationProvider, structure, { structure })],
    ['molecular-surface', () => createRepr(MolecularSurfaceRepresentationProvider, structure, { structure })],
    ['gaussian-surface', () => createRepr(GaussianSurfaceRepresentationProvider, structure, { structure })],
    [
      'volume isosurface',
      async () => {
        const volume = createTestVolume();
        return createRepr(IsosurfaceRepresentationProvider, volume, { volume });
      },
    ],
    [
      'shape mesh',
      async () => {
        const shape = createTestShape();
        const repr = ShapeRepresentation((_: unknown, s: Shape<Mesh>) => s, Mesh.Utils);
        await repr.createOrUpdate({}, shape).run();
        return repr as unknown as Representation<any>;
      },
    ],
  ];

  for (const [name, create] of cases) {
    it(`${name}: touches only clip values`, async () => {
      const repr = await create();
      const ros = [...repr.renderObjects];
      expect(ros.length).toBeGreaterThan(0);

      const before = versions(ros);
      for (let i = 0; i < 3; i++) await repr.createOrUpdate({ ...repr.props, clip: plane(i) } as any).run();

      expect(repr.renderObjects.map((ro) => ro.id)).toEqual(ros.map((ro) => ro.id));
      const changed = changedValues(before, repr.renderObjects);
      expect(changed.filter((k) => !ClipValues.test(k))).toEqual([]);
      expect(changed).toContain('uClipObjectCount');
    });
  }
});
