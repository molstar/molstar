/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Ludovic Autin <autin@scripps.edu>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { StateTransforms } from '../transforms';
import { DataFormatProvider } from './provider';
import { PluginContext } from '../../mol-plugin/context';
import { StateBuilder, StateObjectRef } from '../../mol-state';
import { PluginStateObject } from '../objects';
import { looksLikeMmcifParticles, MmcifVariant } from '../../mol-model-formats/particles/mmcif';
import { getParticleTargetGroups, ParticleList } from '../../mol-model/particles/particle-list';

export const ParticlesFormatCategory = 'Particles';

interface ParticleFormatData {
    format: StateObjectRef
    list: StateObjectRef<PluginStateObject.Particle.List>
}

/** Whether the particle list has target objects and whether some particles are without one. */
function getParticleTargetCoverage(particles: ParticleList | undefined) {
    const mapping = particles?.targetMapping;
    if (!particles || !mapping || mapping.size === 0) {
        return { hasTargets: false, hasUntargeted: !!particles && particles.count > 0 };
    }

    const { targetIds } = getParticleTargetGroups(particles);
    let hasTargets = false;
    let hasUntargeted = false;
    for (let i = 0; i < targetIds.length; ++i) {
        if (mapping.has(targetIds[i])) hasTargets = true;
        else hasUntargeted = true;
    }
    return { hasTargets, hasUntargeted };
}

/**
 * Adds the `target` representation for particles with a reference object and the `spacefill`
 * representation for the remaining ones, which also covers lists without any target.
 */
function addParticleShapeRepresentations(builder: StateBuilder.Root, data: ParticleFormatData, particles: ParticleList | undefined, type: 'spacefill' | 'target', targetProps?: { params?: {}, colorTheme?: { name: string, params: {} }, sizeTheme?: { name: string, params: {} } }) {
    const { hasTargets, hasUntargeted } = getParticleTargetCoverage(particles);
    const useTarget = type === 'target' && hasTargets;

    if (useTarget) {
        builder.to(data.list)
            .apply(StateTransforms.Particles.ParticlesRepresentation3D, {
                type: { name: 'target', params: targetProps?.params ?? {} },
                ...(targetProps?.colorTheme && { colorTheme: targetProps.colorTheme }),
                ...(targetProps?.sizeTheme && { sizeTheme: targetProps.sizeTheme }),
            });
    }
    if (!useTarget || hasUntargeted) {
        builder.to(data.list)
            .apply(StateTransforms.Particles.ParticlesRepresentation3D, { type: { name: 'spacefill', params: { excludeTargets: useTarget } } });
    }
}

function particleVisuals(plugin: PluginContext, data: ParticleFormatData, type: 'spacefill' | 'target' = 'spacefill') {
    const builder = plugin.state.data.build();
    const particleList = StateObjectRef.resolveAndCheck(plugin.state.data, data.list)?.obj?.data;

    addParticleShapeRepresentations(builder, data, particleList, type);

    builder.to(data.list)
        .apply(StateTransforms.Particles.ParticleListUnitcell3D, { attachment: 'center' });

    if (particleList?.fibers && particleList.fibers.count > 0) {
        builder.to(data.list)
            .apply(StateTransforms.Particles.ParticlesRepresentation3D, { type: { name: 'fibers', params: {} } });
    }

    return builder.commit();
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
    stringExtensions: ['ndjson'],
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

export const ArtiatomiEmParticlesProvider = DataFormatProvider({
    label: 'Artiatomi EM Particles',
    description: 'Artiatomi EM Particles',
    category: ParticlesFormatCategory,
    binaryExtensions: ['em'],
    parse: async (plugin, data, params?: { label?: string, tomo?: number, pixelSize?: number }) => {
        const format = plugin.state.data.build()
            .to(data)
            .apply(StateTransforms.Data.ParseArtiatomiEm, void 0, { state: { isGhost: true } });

        const list = format.apply(StateTransforms.Particles.ParticleListFromArtiatomiEm, {
            tomos: params?.tomo !== void 0 ? [String(params.tomo)] : [],
            pixelSize: params?.pixelSize ?? 1,
        });

        await format.commit({ revertOnError: true });

        return { format: format.selector, list: list.selector };
    },
    visuals: particleVisuals,
});

export const MmcifParticlesProvider = DataFormatProvider({
    label: 'mmCIF Particles',
    description: 'mmCIF Particles (CellPack / PetWorld assemblies)',
    category: ParticlesFormatCategory,
    stringExtensions: ['cif'],
    binaryExtensions: ['bcif'],
    /**
     * Higher than the default (0) mmCIF/CifCore trajectory providers so that CellPack/PetWorld
     * assemblies are recognized ahead of them during auto-detection, regardless of registration order.
     */
    priority: 10,
    isApplicable: (info, data) => looksLikeMmcifParticles(info, data),
    parse: async (plugin, data, params?: { label?: string, assemblyId?: string, asymIds?: string[], variant?: MmcifVariant }) => {
        const format = plugin.state.data.build()
            .to(data)
            .apply(StateTransforms.Data.ParseCif, void 0, { state: { isGhost: true } });

        const list = format.apply(StateTransforms.Particles.ParticleListFromMmcifAssembly, {
            assemblyId: params?.assemblyId,
            asymIds: params?.asymIds,
            variant: params?.variant,
            label: params?.label,
        });

        await format.commit({ revertOnError: true });

        return { format: format.selector, list: list.selector };
    },
    visuals: (plugin: PluginContext, data: ParticleFormatData) => {
        const builder = plugin.state.data.build();
        const particleList = StateObjectRef.resolveAndCheck(plugin.state.data, data.list)?.obj?.data;

        addParticleShapeRepresentations(builder, data, particleList, 'target', {
            params: { quality: 'lowest' },
            colorTheme: { name: 'particle-entity', params: {} },
            sizeTheme: { name: 'particle-size', params: { scale: 1 } },
        });

        if (particleList?.fibers && particleList.fibers.count > 0) {
            builder.to(data.list)
                .apply(StateTransforms.Particles.ParticlesRepresentation3D, { type: { name: 'fibers', params: {} } });
        }

        return builder.commit();
    },
});

export const SimulariumParticlesProvider = DataFormatProvider({
    label: 'Simularium Particles',
    description: 'Simularium Particles',
    category: ParticlesFormatCategory,
    binaryExtensions: ['simularium'],
    parse: async (plugin, data, params?: { label?: string, frameIndex?: number, loadGeometries?: boolean }) => {
        const format = plugin.state.data.build()
            .to(data)
            .apply(StateTransforms.Data.ParseSimularium, void 0, { state: { isGhost: false } });

        const trajectory = format.apply(StateTransforms.Particles.ParticleTrajectoryFromSimularium, {
            loadGeometries: params?.loadGeometries ?? true,
        });

        const list = trajectory.apply(StateTransforms.Particles.ParticleListFromTrajectory, {
            frameIndex: params?.frameIndex ?? 0,
        });

        await format.commit({ revertOnError: true });

        return { format: format.selector, trajectory: trajectory.selector, list: list.selector };
    },
    // The `target` representation instances the geometries referenced by the agent types;
    // particles without a geometry are added as spacefill spheres.
    visuals: (plugin: PluginContext, data: ParticleFormatData) => particleVisuals(plugin, data, 'target'),
});

export const BuiltInParticlesFormats = [
    ['relion_star_particles', RelionStarParticlesProvider] as const,
    ['dynamo_tbl_particles', DynamoTblParticlesProvider] as const,
    ['cryoet_ndjson_particles', CryoEtDataPortalNdjsonParticlesProvider] as const,
    ['artiatomi_em_particles', ArtiatomiEmParticlesProvider] as const,
    ['mmcif_particles', MmcifParticlesProvider] as const,
    ['simularium_particles', SimulariumParticlesProvider] as const,
] as const;

export type BuiltInParticlesFormat = (typeof BuiltInParticlesFormats)[number][0]
