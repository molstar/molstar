/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Neli Fonseca <neli@ebi.ac.uk>
 * @author Ludovic Autin <autin@scripps.edu>
 */

import { PluginStateTransform, PluginStateObject as SO, type PluginStateObject } from '@molstar/plugin/state/objects';
import { Task } from '@molstar/core/task';
import * as PLY from '@molstar/io/reader/ply/parser';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { Mat4 } from '@molstar/core/math/linear-algebra';
import { shapeFromPly } from '@molstar/graphics/formats/shape/ply';
import { DataFormatProvider, applyTransformerRaw, rawDataObject } from '@molstar/plugin/state/formats/provider';
import { ShapeFormatCategory } from './category.js';
import type { PluginContext } from '@molstar/plugin/context';
import type { StateObjectRef } from '@molstar/core/state';
import { ShapeRepresentation3D } from '@molstar/plugin/state/transforms/shape/representation';

export { ParsePly };
type ParsePly = typeof ParsePly;
const ParsePly = PluginStateTransform.BuiltIn({
  name: 'parse-ply',
  display: { name: 'Parse PLY', description: 'Parse PLY from String or Binary data' },
  from: [SO.Data.String, SO.Data.Binary],
  to: SO.Format.Ply,
})({
  apply({ a }) {
    return Task.create('Parse PLY', async (ctx) => {
      const parsed = await PLY.parsePly(a.data).runInContext(ctx);
      if (parsed.isError) throw new Error(parsed.message);
      return new SO.Format.Ply(parsed.result, { label: parsed.result.comments[0] || 'PLY Data' });
    });
  },
});

export { ShapeFromPly };
type ShapeFromPly = typeof ShapeFromPly;
const ShapeFromPly = PluginStateTransform.BuiltIn({
  name: 'shape-from-ply',
  display: { name: 'Shape from PLY', description: 'Create Shape from PLY data' },
  from: SO.Format.Ply,
  to: SO.Shape.Provider,
  params(a) {
    return {
      transforms: PD.Optional(PD.Value([Mat4.identity()], { isHidden: true })),
      label: PD.Optional(PD.Text('', { isHidden: true })),
    };
  },
})({
  apply({ a, params }) {
    return Task.create('Create shape from PLY', async (ctx) => {
      const shape = await shapeFromPly(a.data, params).runInContext(ctx);
      const props = { label: params.label || 'Shape' };
      return new SO.Shape.Provider(shape, props);
    });
  },
});

export const PlyProvider = DataFormatProvider({
  name: 'ply',
  label: 'PLY',
  description: 'PLY',
  category: ShapeFormatCategory,
  stringExtensions: ['ply'],
  binaryExtensions: ['ply'], // binary files have same extension
  parse: async (plugin, data) => {
    const format = plugin.state.data
      .build()
      .to(data)
      .apply(ParsePly, {}, { state: { isGhost: true } });

    const shape = format.apply(ShapeFromPly);

    await format.commit();

    return { format: format.selector, shape: shape.selector };
  },
  parseRaw: async (plugin, ctx, data) => {
    const format = await applyTransformerRaw(plugin, ctx, ParsePly, rawDataObject(data));
    const shape = await applyTransformerRaw(plugin, ctx, ShapeFromPly, format);
    return { shape: shape.data };
  },
  visuals(plugin: PluginContext, data: { shape: StateObjectRef<PluginStateObject.Shape.Provider> }) {
    const repr = plugin.state.data.build().to(data.shape).apply(ShapeRepresentation3D);
    return repr.commit();
  },
});
