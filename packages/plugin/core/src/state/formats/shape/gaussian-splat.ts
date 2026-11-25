/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { PluginStateTransform, PluginStateObject as SO, type PluginStateObject } from '@molstar/plugin/state/objects';
import { Task } from '@molstar/core/task';
import { Mat4 } from '@molstar/core/math/linear-algebra';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { shapeFromGaussianSplats } from '@molstar/graphics/formats/shape/gaussian-splat';

/** Gaussian splat files besides PLY (`ShapeFromPly`): antimatter15's `.splat`, Niantic's `.spz`, PlayCanvas's `.sog` */
type GaussianSplatFileKind = 'splat' | 'spz' | 'sog';
import { DataFormatProvider, applyTransformerRaw, rawDataObject } from '@molstar/plugin/state/formats/provider';
import { ShapeFormatCategory } from './category.js';
import type { PluginContext } from '@molstar/plugin/context';
import type { StateObjectRef } from '@molstar/core/state';
import { ShapeRepresentation3D } from '@molstar/plugin/state/transforms/shape/representation';
import type { PluginRegistryEntry } from '@molstar/plugin/spec';

const _identityTransforms = [Mat4.identity()];

export { ShapeFromGaussianSplatFile };
type ShapeFromGaussianSplatFile = typeof ShapeFromGaussianSplatFile;
/** A shape of gaussian splats from the bytes of a `.splat` or `.spz` file (PLY files: `ShapeFromPly`) */
const ShapeFromGaussianSplatFile = PluginStateTransform.BuiltIn({
  name: 'shape-from-gaussian-splat-file',
  display: { name: 'Shape from Gaussian Splats', description: 'Create Shape from a .splat, .spz or .sog file' },
  from: [SO.Data.Binary, SO.Data.Response],
  to: SO.Shape.Provider,
  params(a) {
    return {
      kind: PD.Select<GaussianSplatFileKind>('splat', [
        ['splat', 'SPLAT'],
        ['spz', 'SPZ'],
        ['sog', 'SOG'],
      ]),
      transforms: PD.Optional(PD.Value(_identityTransforms, { isHidden: true })),
      label: PD.Optional(PD.Text('', { isHidden: true })),
    };
  },
})({
  apply({ a, params }) {
    return Task.create('Create shape from gaussian splats', async (ctx) => {
      // the bytes as they are: read, never written (the parsers keep views of them)
      const bytes = SO.Data.readBytes(a);
      const shape = await shapeFromGaussianSplats({
        kind: params.kind,
        bytes,
        transforms: params.transforms,
      }).runInContext(ctx);
      return new SO.Shape.Provider(shape, { label: params.label || 'Gaussian Splats' });
    });
  },
});

function GaussianSplatFileProvider(kind: GaussianSplatFileKind, description: string) {
  return DataFormatProvider({
    name: kind,
    label: kind.toUpperCase(),
    description,
    category: ShapeFormatCategory,
    binaryExtensions: [kind],
    parse: async (plugin, data) => {
      const shape = plugin.state.data.build().to(data).apply(ShapeFromGaussianSplatFile, { kind });
      await shape.commit();
      return { shape: shape.selector };
    },
    parseRaw: async (plugin, ctx, data) => {
      const shape = await applyTransformerRaw(plugin, ctx, ShapeFromGaussianSplatFile, rawDataObject(data), { kind });
      return { shape: shape.data };
    },
    visuals(plugin: PluginContext, data: { shape: StateObjectRef<PluginStateObject.Shape.Provider> }) {
      const repr = plugin.state.data.build().to(data.shape).apply(ShapeRepresentation3D);
      return repr.commit();
    },
  });
}

export const SplatProvider = GaussianSplatFileProvider('splat', 'Gaussian splats (antimatter15 .splat)');
export const SpzProvider = GaussianSplatFileProvider('spz', 'Gaussian splats (Niantic .spz)');
export const SogProvider = GaussianSplatFileProvider('sog', 'Gaussian splats (PlayCanvas .sog)');

/** The Splat data format (antimatter15's gaussian splats; PLY files of splats: `Ply`). */
export const Splat: PluginRegistryEntry = {
  formats: [SplatProvider],
};

/** The Spz data format (Niantic's gaussian splats). */
export const Spz: PluginRegistryEntry = {
  formats: [SpzProvider],
};

/** The Sog data format (PlayCanvas's gaussian splats; its images decoded by the browser). */
export const Sog: PluginRegistryEntry = {
  formats: [SogProvider],
};
