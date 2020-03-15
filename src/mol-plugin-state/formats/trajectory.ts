/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { FileInfo } from '../../mol-util/file-info';
import { StateTransforms } from '../transforms';
import { guessCifVariant, DataFormatProvider, DataFormatRegistry } from './registry';
import { StateTransformer, StateObjectRef } from '../../mol-state';
import { PluginStateObject } from '../objects';

export interface TrajectoryFormatProvider<P extends { trajectoryTags?: string | string[] } = { trajectoryTags?: string | string[] }, R extends { trajectory: StateObjectRef<PluginStateObject.Molecule.Trajectory> } = { trajectory: StateObjectRef<PluginStateObject.Molecule.Trajectory> }> 
    extends DataFormatProvider<P, R> {
}

export function TrajectoryFormatRegistry() {
    return new DataFormatRegistry<TrajectoryFormatProvider>(BuildInTrajectoryFormats);
}

export const MmcifProvider: TrajectoryFormatProvider = {
    label: 'mmCIF',
    description: 'mmCIF',
    stringExtensions: ['cif', 'mmcif', 'mcif'],
    binaryExtensions: ['bcif'],
    isApplicable: (info: FileInfo, data: Uint8Array | string) => {
        if (info.ext === 'mmcif' || info.ext === 'mcif') return true
        // assume cif/bcif files that are not DensityServer CIF are mmCIF
        if (info.ext === 'cif' || info.ext === 'bcif') return guessCifVariant(info, data) !== 'dscif'
        return false
    },
    parse: async (plugin, data, params) => {
        const state = plugin.state.dataState;
        const trajectory = state.build().to(data)
            .apply(StateTransforms.Data.ParseCif, void 0, { state: { isGhost: true } })
            .apply(StateTransforms.Model.TrajectoryFromMmCif, void 0, { tags: params?.trajectoryTags })
        await plugin.runTask(state.updateTree(trajectory, { revertOnError: true }));
        return { trajectory: trajectory.selector };
    }
}

function directTrajectory(transformer: StateTransformer<PluginStateObject.Data.String | PluginStateObject.Data.Binary, PluginStateObject.Molecule.Trajectory>): TrajectoryFormatProvider['parse'] {
    return async (plugin, data, params) => {
        const state = plugin.state.dataState;
        const trajectory = state.build().to(data)
            .apply(transformer, void 0, { tags: params?.trajectoryTags })
        await plugin.runTask(state.updateTree(trajectory, { revertOnError: true }));
        return { trajectory: trajectory.selector };
    }
}

export const PdbProvider: TrajectoryFormatProvider = {
    label: 'PDB',
    description: 'PDB',
    stringExtensions: ['pdb', 'ent'],
    binaryExtensions: [],
    isApplicable: (info: FileInfo, data: string) => {
        return info.ext === 'pdb' || info.ext === 'ent'
    },
    parse: directTrajectory(StateTransforms.Model.TrajectoryFromPDB)
}

export const GroProvider: TrajectoryFormatProvider = {
    label: 'GRO',
    description: 'GRO',
    stringExtensions: ['gro'],
    binaryExtensions: [],
    isApplicable: (info: FileInfo, data: string) => {
        return info.ext === 'gro'
    },
    parse: directTrajectory(StateTransforms.Model.TrajectoryFromGRO)
}

export const Provider3dg: TrajectoryFormatProvider = {
    label: '3DG',
    description: '3DG',
    stringExtensions: ['3dg'],
    binaryExtensions: [],
    isApplicable: (info: FileInfo, data: string) => {
        return info.ext === '3dg'
    },
    parse: directTrajectory(StateTransforms.Model.TrajectoryFrom3DG)
}

export const BuildInTrajectoryFormats = [
    ['mmcif', MmcifProvider] as const,
    ['pdb', PdbProvider] as const,
    ['gro', GroProvider] as const,
    ['3dg', Provider3dg] as const,
] as const

export type BuiltInTrajectoryFormat = (typeof BuildInTrajectoryFormats)[number][0]