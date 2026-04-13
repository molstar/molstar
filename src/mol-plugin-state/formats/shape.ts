/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
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
    parse: async (plugin, data) => {
        const format = plugin.state.data.build()
            .to(data)
            .apply(StateTransforms.Data.ParsePly, {}, { state: { isGhost: true } });

        const shape = format.apply(StateTransforms.Model.ShapeFromPly);

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
] as const;

export type BuildInShapeFormat = (typeof BuiltInShapeFormats)[number][0]