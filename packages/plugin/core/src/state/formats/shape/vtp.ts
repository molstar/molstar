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
import * as VTP from '@molstar/io/reader/vtp/parser';
import { Mat4 } from '@molstar/core/math/linear-algebra';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { shapeFromVtp } from '@molstar/graphics/formats/shape/vtp';
import { DataFormatProvider, applyTransformerRaw, rawDataObject } from '@molstar/plugin/state/formats/provider';
import { ShapeFormatCategory } from './category.js';
import type { PluginContext } from '@molstar/plugin/context';
import type { StateObjectRef } from '@molstar/core/state';
import { ShapeRepresentation3D } from '@molstar/plugin/state/transforms/shape/representation';

export { ParseVtp };
type ParseVtp = typeof ParseVtp;
const ParseVtp = PluginStateTransform.BuiltIn({
  name: 'parse-vtp',
  display: { name: 'Parse VTP', description: 'Parse VTP (VTK PolyData) from Binary data' },
  from: [SO.Data.Binary],
  to: SO.Format.Vtp,
})({
  apply({ a }) {
    return Task.create('Parse VTP', async (ctx) => {
      const parsed = await VTP.parseVtp(a.data).runInContext(ctx);
      if (parsed.isError) throw new Error(parsed.message);
      return new SO.Format.Vtp(parsed.result, { label: 'VTP Data' });
    });
  },
});

const _vtpIdentityTransforms = [Mat4.identity()];

export { ShapeFromVtp };
type ShapeFromVtp = typeof ShapeFromVtp;
const ShapeFromVtp = PluginStateTransform.BuiltIn({
  name: 'shape-from-vtp',
  display: { name: 'Shape from VTP', description: 'Create Shape from VTP (VTK PolyData) file' },
  from: SO.Format.Vtp,
  to: SO.Shape.Provider,
  params(a) {
    return {
      transforms: PD.Optional(PD.Value(_vtpIdentityTransforms, { isHidden: true })),
      label: PD.Optional(PD.Text('', { isHidden: true })),
    };
  },
})({
  apply({ a, params }) {
    return Task.create('Create shape from VTP', async (ctx) => {
      const shape = await shapeFromVtp(a.data, params).runInContext(ctx);
      const props = { label: params.label || 'VTP Shape' };
      return new SO.Shape.Provider(shape, props);
    });
  },
});

export const VtpProvider = DataFormatProvider({
  name: 'vtp',
  label: 'VTP',
  description: 'VTK PolyData (VTP)',
  category: ShapeFormatCategory,
  binaryExtensions: ['vtp'],
  parse: async (plugin, data) => {
    const format = plugin.state.data
      .build()
      .to(data)
      .apply(ParseVtp, {}, { state: { isGhost: true } });

    const shape = format.apply(ShapeFromVtp);

    await format.commit();

    return { format: format.selector, shape: shape.selector };
  },
  parseRaw: async (plugin, ctx, data) => {
    const format = await applyTransformerRaw(plugin, ctx, ParseVtp, rawDataObject(data));
    const shape = await applyTransformerRaw(plugin, ctx, ShapeFromVtp, format);
    return { shape: shape.data };
  },
  visuals(plugin: PluginContext, data: { shape: StateObjectRef<PluginStateObject.Shape.Provider> }) {
    const repr = plugin.state.data.build().to(data.shape).apply(ShapeRepresentation3D);
    return repr.commit();
  },
});
