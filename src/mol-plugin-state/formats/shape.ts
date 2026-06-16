/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { StateTransforms } from '../transforms';
import { DataFormatProvider } from './provider';
import { PluginContext } from '../../mol-plugin/context';
import { StateObjectRef } from '../../mol-state';
import { PluginStateObject } from '../objects';

export const ShapeFormatCategory = 'Shape';

export const PlyProvider = DataFormatProvider({
    label: 'PLY',
    description: 'PLY',
    category: ShapeFormatCategory,
    stringExtensions: ['ply'],
    binaryExtensions: ['ply'], // binary files have same extension
    parse: async (plugin, data) => {
        const format = plugin.state.data.build()
            .to(data)
            .apply(StateTransforms.Data.ParsePly, {}, { state: { isGhost: true } });

        const shape = format.apply(StateTransforms.Shape.ShapeFromPly);

        await format.commit();

        return { format: format.selector, shape: shape.selector };
    },
    visuals(plugin: PluginContext, data: { shape: StateObjectRef<PluginStateObject.Shape.Provider> }) {
        const repr = plugin.state.data.build()
            .to(data.shape)
            .apply(StateTransforms.Representation.ShapeRepresentation3D);
        return repr.commit();
    }
});

export const ObjProvider = DataFormatProvider({
    label: 'OBJ',
    description: 'OBJ',
    category: ShapeFormatCategory,
    stringExtensions: ['obj'],
    parse: async (plugin, data) => {
        const format = plugin.state.data.build()
            .to(data)
            .apply(StateTransforms.Data.ParseObj, {}, { state: { isGhost: true } });

        const shape = format.apply(StateTransforms.Shape.ShapeFromObj);

        await format.commit();

        return { format: format.selector, shape: shape.selector };
    },
    visuals(plugin: PluginContext, data: { shape: StateObjectRef<PluginStateObject.Shape.Provider> }) {
        const repr = plugin.state.data.build()
            .to(data.shape)
            .apply(StateTransforms.Representation.ShapeRepresentation3D);
        return repr.commit();
    }
});

export const VtpProvider = DataFormatProvider({
    label: 'VTP',
    description: 'VTK PolyData (VTP)',
    category: ShapeFormatCategory,
    binaryExtensions: ['vtp'],
    parse: async (plugin, data) => {
        const format = plugin.state.data.build()
            .to(data)
            .apply(StateTransforms.Data.ParseVtp, {}, { state: { isGhost: true } });

        const shape = format.apply(StateTransforms.Model.ShapeFromVtp);

        await format.commit();

        return { format: format.selector, shape: shape.selector };
    },
    visuals(plugin: PluginContext, data: { shape: StateObjectRef<PluginStateObject.Shape.Provider> }) {
        const repr = plugin.state.data.build()
            .to(data.shape)
            .apply(StateTransforms.Representation.ShapeRepresentation3D);
        return repr.commit();
    }
});

export const BuiltInShapeFormats = [
    ['ply', PlyProvider] as const,
    ['obj', ObjProvider] as const,
    ['vtp', VtpProvider] as const,
] as const;

export type BuildInShapeFormat = (typeof BuiltInShapeFormats)[number][0]