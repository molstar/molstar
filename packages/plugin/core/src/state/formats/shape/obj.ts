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
import * as OBJ from '@molstar/io/reader/obj/parser';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { Mat4 } from '@molstar/core/math/linear-algebra';
import type { PluginContext } from '@molstar/plugin/context';
import { parseMtl } from '@molstar/io/reader/obj/mtl-parser';
import { shapeFromObj } from '@molstar/graphics/formats/shape/obj';
import type { Asset } from '@molstar/core/util/assets';
import { DataFormatProvider, applyTransformerRaw, rawDataObject } from '@molstar/plugin/state/formats/provider';
import { ShapeFormatCategory } from './category.js';
import type { StateObjectRef } from '@molstar/core/state';
import { ShapeRepresentation3D } from '@molstar/plugin/state/transforms/shape/representation';
import type { PluginRegistryEntry } from '@molstar/plugin/spec';

export { ParseObj };
type ParseObj = typeof ParseObj;
const ParseObj = PluginStateTransform.BuiltIn({
  name: 'parse-obj',
  display: { name: 'Parse OBJ', description: 'Parse OBJ from String data' },
  from: [SO.Data.String],
  to: SO.Format.Obj,
})({
  apply({ a }) {
    return Task.create('Parse OBJ', async (ctx) => {
      const parsed = await OBJ.parseObj(a.data).runInContext(ctx);
      if (parsed.isError) throw new Error(parsed.message);
      return new SO.Format.Obj(parsed.result, { label: 'OBJ Data' });
    });
  },
});

export { ShapeFromObj };
type ShapeFromObj = typeof ShapeFromObj;
const ShapeFromObj = PluginStateTransform.BuiltIn({
  name: 'shape-from-obj',
  display: { name: 'Shape from OBJ', description: 'Create Shape from OBJ data' },
  from: SO.Format.Obj,
  to: SO.Shape.Provider,
  params(a) {
    return {
      transforms: PD.Optional(PD.Value([Mat4.identity()], { isHidden: true })),
      label: PD.Optional(PD.Text('', { isHidden: true })),
      mtlFile: PD.Optional(PD.File({ accept: '.mtl', label: 'MTL File' })),
    };
  },
})({
  apply({ a, params, cache }, plugin: PluginContext) {
    return Task.create('Create shape from OBJ', async (ctx) => {
      let mtl;
      if (params.mtlFile) {
        const asset = await plugin.managers.asset.resolve(params.mtlFile, 'string').runInContext(ctx);
        (cache as any).mtlAsset = asset;
        mtl = parseMtl(asset.data as string);
      }
      const shape = await shapeFromObj(a.data, { ...params, mtl }).runInContext(ctx);
      const props = { label: params.label || 'Shape' };
      return new SO.Shape.Provider(shape, props);
    });
  },
  dispose({ cache }) {
    ((cache as any)?.mtlAsset as Asset.Wrapper | undefined)?.dispose();
  },
});

export const ObjProvider = DataFormatProvider({
  name: 'obj',
  label: 'OBJ',
  description: 'OBJ',
  category: ShapeFormatCategory,
  stringExtensions: ['obj'],
  parse: async (plugin, data) => {
    const format = plugin.state.data
      .build()
      .to(data)
      .apply(ParseObj, {}, { state: { isGhost: true } });

    const shape = format.apply(ShapeFromObj);

    await format.commit();

    return { format: format.selector, shape: shape.selector };
  },
  parseRaw: async (plugin, ctx, data) => {
    const format = await applyTransformerRaw(plugin, ctx, ParseObj, rawDataObject(data));
    const shape = await applyTransformerRaw(plugin, ctx, ShapeFromObj, format);
    return { shape: shape.data };
  },
  visuals(plugin: PluginContext, data: { shape: StateObjectRef<PluginStateObject.Shape.Provider> }) {
    const repr = plugin.state.data.build().to(data.shape).apply(ShapeRepresentation3D);
    return repr.commit();
  },
});

/** The Obj data format. */
export const Obj: PluginRegistryEntry = {
  formats: [ObjProvider],
};
