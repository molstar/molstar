/**
 * Copyright (c) 2018-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { StateTransforms } from '../transforms';
import { DataFormatProvider } from './provider';

export const TopologyFormatCategory = 'Topology';

export { PsfProvider };
const PsfProvider = DataFormatProvider({
    label: 'PSF',
    description: 'PSF',
    category: TopologyFormatCategory,
    stringExtensions: ['psf'],
    parse: async (plugin, data) => {
        const format = plugin.state.data.build()
            .to(data)
            .apply(StateTransforms.Data.ParsePsf, {}, { state: { isGhost: true } });
        const topology = format.apply(StateTransforms.Model.TopologyFromPsf);

        await format.commit();

        return { format: format.selector, topology: topology.selector };
    }
});
type PsfProvider = typeof PsfProvider;

export type TopologyProvider = PsfProvider;

export const BuiltInTopologyFormats = [
    ['psf', PsfProvider] as const,
] as const;

export type BuiltInTopologyFormat = (typeof BuiltInTopologyFormats)[number][0]
