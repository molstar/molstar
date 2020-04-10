/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { StateTransforms } from '../transforms';
import { DataFormatProvider } from './provider';

const Category = 'Structure';

export const PsfProvider = DataFormatProvider({
    label: 'PSF',
    description: 'PSF',
    category: Category,
    stringExtensions: ['psf'],
    parse: async (plugin, data) => {
        const format = plugin.state.data.build()
            .to(data)
            .apply(StateTransforms.Data.ParsePsf, {}, { state: { isGhost: true } });
        const topology = format.apply(StateTransforms.Model.TopologyFromPsf);

        await plugin.updateDataState(format);

        return { format: format.selector, topology: topology.selector };
    }
});

export const DcdProvider = DataFormatProvider({
    label: 'DCD',
    description: 'DCD',
    category: Category,
    binaryExtensions: ['dcd'],
    parse: async (plugin, data) => {
        const coordinates = plugin.state.data.build()
            .to(data)
            .apply(StateTransforms.Model.CoordinatesFromDcd);

        await plugin.updateDataState(coordinates);

        return { coordinates: coordinates.selector };
    }
});

export const BuiltInStructureFormats = [
    ['psf', PsfProvider] as const,
    ['dcd', DcdProvider] as const,
] as const

export type BuildInStructureFormat = (typeof BuiltInStructureFormats)[number][0]