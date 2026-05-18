/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Ludovic Autin <autin@scripps.edu>
 */

import { StateTransforms } from '../transforms';
import { DataFormatProvider } from './provider';
import { PluginContext } from '../../mol-plugin/context';
import { StateObjectRef } from '../../mol-state';
import { PluginStateObject } from '../objects';

export const ParticlesFormatCategory = 'Particles';

interface ParticleFormatData {
    format: StateObjectRef
    list: StateObjectRef<PluginStateObject.Particle.List>
}

function particleVisuals(plugin: PluginContext, data: ParticleFormatData) {
    const repr = plugin.state.data.build()
        .to(data.list)
        .apply(StateTransforms.Particles.ParticlesRepresentation3D);
    return repr.commit();
}

export const RelionStarParticlesProvider = DataFormatProvider({
    label: 'RELION STAR Particles',
    description: 'RELION STAR Particles',
    category: ParticlesFormatCategory,
    stringExtensions: ['star'],
    parse: async (plugin, data, params?: { label?: string, tomogram?: string }) => {
        const format = plugin.state.data.build()
            .to(data)
            .apply(StateTransforms.Data.ParseCif, void 0, { state: { isGhost: true } });

        const list = format.apply(StateTransforms.Particles.ParticleListFromRelionStar, {
            tomograms: params?.tomogram ? [params.tomogram] : [],
        });

        await format.commit({ revertOnError: true });

        return { format: format.selector, list: list.selector };
    },
    visuals: particleVisuals,
});

export const DynamoTblParticlesProvider = DataFormatProvider({
    label: 'Dynamo TBL Particles',
    description: 'Dynamo TBL Particles',
    category: ParticlesFormatCategory,
    stringExtensions: ['tbl'],
    parse: async (plugin, data, params?: { label?: string, tomo?: number }) => {
        const format = plugin.state.data.build()
            .to(data)
            .apply(StateTransforms.Data.ParseDynamoTbl, void 0, { state: { isGhost: true } });

        const list = format.apply(StateTransforms.Particles.ParticleListFromDynamoTbl, {
            tomos: params?.tomo !== void 0 ? [String(params.tomo)] : [],
        });

        await format.commit({ revertOnError: true });

        return { format: format.selector, list: list.selector };
    },
    visuals: particleVisuals,
});

export const CryoEtDataPortalNdjsonParticlesProvider = DataFormatProvider({
    label: 'CryoET NDJSON Particles',
    description: 'CryoET NDJSON Particles',
    category: ParticlesFormatCategory,
    stringExtensions: ['ndjson', 'jsonl'],
    parse: async (plugin, data, params?: { label?: string, type?: string }) => {
        const format = plugin.state.data.build()
            .to(data)
            .apply(StateTransforms.Data.ParseCryoEtDataPortalNdjson, void 0, { state: { isGhost: true } });

        const list = format.apply(StateTransforms.Particles.ParticleListFromCryoEtDataPortalNdjson, {
            type: params?.type,
        });

        await format.commit({ revertOnError: true });

        return { format: format.selector, list: list.selector };
    },
    visuals: particleVisuals,
});

export const BuiltInParticlesFormats = [
    ['relion_star_particles', RelionStarParticlesProvider] as const,
    ['dynamo_tbl_particles', DynamoTblParticlesProvider] as const,
    ['cryoet_ndjson_particles', CryoEtDataPortalNdjsonParticlesProvider] as const,
] as const;

export type BuiltInParticlesFormat = (typeof BuiltInParticlesFormats)[number][0]
