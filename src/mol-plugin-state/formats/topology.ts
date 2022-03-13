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

export { PrmtopProvider };
const PrmtopProvider = DataFormatProvider({
    label: 'PRMTOP',
    description: 'PRMTOP',
    category: TopologyFormatCategory,
    stringExtensions: ['prmtop', 'parm7'],
    parse: async (plugin, data) => {
        const format = plugin.state.data.build()
            .to(data)
            .apply(StateTransforms.Data.ParsePrmtop, {}, { state: { isGhost: true } });
        const topology = format.apply(StateTransforms.Model.TopologyFromPrmtop);

        await format.commit();

        return { format: format.selector, topology: topology.selector };
    }
});
type PrmtopProvider = typeof PrmtopProvider;

export { TopProvider };
const TopProvider = DataFormatProvider({
    label: 'TOP',
    description: 'TOP',
    category: TopologyFormatCategory,
    stringExtensions: ['top'],
    parse: async (plugin, data) => {
        const format = plugin.state.data.build()
            .to(data)
            .apply(StateTransforms.Data.ParseTop, {}, { state: { isGhost: true } });
        const topology = format.apply(StateTransforms.Model.TopologyFromTop);

        await format.commit();

        return { format: format.selector, topology: topology.selector };
    }
});
type TopProvider = typeof TopProvider;

export type TopologyProvider = PsfProvider;

export const BuiltInTopologyFormats = [
    ['psf', PsfProvider] as const,
    ['prmtop', PrmtopProvider] as const,
    ['top', TopProvider] as const,
] as const;

export type BuiltInTopologyFormat = (typeof BuiltInTopologyFormats)[number][0]
