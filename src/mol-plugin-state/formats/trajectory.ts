/**
 * Copyright (c) 2018-2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Ludovic Autin <ludovic.autin@gmail.com>
 */

import { StateTransforms } from '../transforms';
import { guessCifVariant, DataFormatProvider } from './provider';
import { StateTransformer, StateObjectRef } from '../../mol-state';
import { PluginStateObject } from '../objects';
import { PluginContext } from '../../mol-plugin/context';

export interface TrajectoryFormatProvider<P extends { trajectoryTags?: string | string[] } = { trajectoryTags?: string | string[] }, R extends { trajectory: StateObjectRef<PluginStateObject.Molecule.Trajectory> } = { trajectory: StateObjectRef<PluginStateObject.Molecule.Trajectory> }>
    extends DataFormatProvider<P, R> {
}

export const TrajectoryFormatCategory = 'Trajectory';

function defaultVisuals(plugin: PluginContext, data: { trajectory: StateObjectRef<PluginStateObject.Molecule.Trajectory> }) {
    return plugin.builders.structure.hierarchy.applyPreset(data.trajectory, 'default');
}

export const MmcifProvider: TrajectoryFormatProvider = {
    label: 'mmCIF',
    description: 'mmCIF',
    category: TrajectoryFormatCategory,
    stringExtensions: ['cif', 'mmcif', 'mcif'],
    binaryExtensions: ['bcif'],
    isApplicable: (info, data) => {
        if (info.ext === 'mmcif' || info.ext === 'mcif') return true;
        // assume undetermined cif/bcif files are mmCIF
        if (info.ext === 'cif' || info.ext === 'bcif') return guessCifVariant(info, data) === -1;
        return false;
    },
    parse: async (plugin, data, params) => {
        const state = plugin.state.data;
        const cif = state.build().to(data)
            .apply(StateTransforms.Data.ParseCif, void 0, { state: { isGhost: true } });
        const trajectory = await cif
            .apply(StateTransforms.Model.TrajectoryFromMmCif, void 0, { tags: params?.trajectoryTags })
            .commit({ revertOnError: true });

        if ((cif.selector.cell?.obj?.data.blocks.length || 0) > 1) {
            plugin.state.data.updateCellState(cif.ref, { isGhost: false });
        }

        return { trajectory };
    },
    visuals: defaultVisuals
};

export const CifCoreProvider: TrajectoryFormatProvider = {
    label: 'cifCore',
    description: 'CIF Core',
    category: TrajectoryFormatCategory,
    stringExtensions: ['cif'],
    isApplicable: (info, data) => {
        if (info.ext === 'cif') return guessCifVariant(info, data) === 'coreCif';
        return false;
    },
    parse: async (plugin, data, params) => {
        const state = plugin.state.data;
        const cif = state.build().to(data)
            .apply(StateTransforms.Data.ParseCif, void 0, { state: { isGhost: true } });
        const trajectory = await cif
            .apply(StateTransforms.Model.TrajectoryFromCifCore, void 0, { tags: params?.trajectoryTags })
            .commit({ revertOnError: true });

        if ((cif.selector.cell?.obj?.data.blocks.length || 0) > 1) {
            plugin.state.data.updateCellState(cif.ref, { isGhost: false });
        }
        return { trajectory };
    },
    visuals: defaultVisuals
};

function directTrajectory<P extends {}>(transformer: StateTransformer<PluginStateObject.Data.String | PluginStateObject.Data.Binary, PluginStateObject.Molecule.Trajectory, P>, transformerParams?: P): TrajectoryFormatProvider['parse'] {
    return async (plugin, data, params) => {
        const state = plugin.state.data;
        const trajectory = await state.build().to(data)
            .apply(transformer, transformerParams, { tags: params?.trajectoryTags })
            .commit({ revertOnError: true });
        return { trajectory };
    };
}

export const PdbProvider: TrajectoryFormatProvider = {
    label: 'PDB',
    description: 'PDB',
    category: TrajectoryFormatCategory,
    stringExtensions: ['pdb', 'ent'],
    parse: directTrajectory(StateTransforms.Model.TrajectoryFromPDB),
    visuals: defaultVisuals
};

export const PdbqtProvider: TrajectoryFormatProvider = {
    label: 'PDBQT',
    description: 'PDBQT',
    category: TrajectoryFormatCategory,
    stringExtensions: ['pdbqt'],
    parse: directTrajectory(StateTransforms.Model.TrajectoryFromPDB, { isPdbqt: true }),
    visuals: defaultVisuals
};

export const XyzProvider: TrajectoryFormatProvider = {
    label: 'XYZ',
    description: 'XYZ',
    category: TrajectoryFormatCategory,
    stringExtensions: ['xyz'],
    parse: directTrajectory(StateTransforms.Model.TrajectoryFromXYZ),
    visuals: defaultVisuals
};

export const LammpsDataProvider: TrajectoryFormatProvider = {
    label: 'Lammps Data',
    description: 'Lammps Data',
    category: TrajectoryFormatCategory,
    stringExtensions: ['data'],
    parse: directTrajectory(StateTransforms.Model.TrajectoryFromLammpsData),
    visuals: defaultVisuals
};

export const LammpsTrajectoryDataProvider: TrajectoryFormatProvider = {
    label: 'Lammps Trajectory Data',
    description: 'Lammps Trajectory Data',
    category: TrajectoryFormatCategory,
    stringExtensions: ['lammpstrj'],
    parse: directTrajectory(StateTransforms.Model.TrajectoryFromLammpsTrajData),
    visuals: defaultVisuals
};

export const GroProvider: TrajectoryFormatProvider = {
    label: 'GRO',
    description: 'GRO',
    category: TrajectoryFormatCategory,
    stringExtensions: ['gro'],
    binaryExtensions: [],
    parse: directTrajectory(StateTransforms.Model.TrajectoryFromGRO),
    visuals: defaultVisuals
};

export const MolProvider: TrajectoryFormatProvider = {
    label: 'MOL',
    description: 'MOL',
    category: TrajectoryFormatCategory,
    stringExtensions: ['mol'],
    parse: directTrajectory(StateTransforms.Model.TrajectoryFromMOL),
    visuals: defaultVisuals
};

export const SdfProvider: TrajectoryFormatProvider = {
    label: 'SDF',
    description: 'SDF',
    category: TrajectoryFormatCategory,
    stringExtensions: ['sdf', 'sd'],
    parse: directTrajectory(StateTransforms.Model.TrajectoryFromSDF),
    visuals: defaultVisuals
};

export const Mol2Provider: TrajectoryFormatProvider = {
    label: 'MOL2',
    description: 'MOL2',
    category: TrajectoryFormatCategory,
    stringExtensions: ['mol2'],
    parse: directTrajectory(StateTransforms.Model.TrajectoryFromMOL2),
    visuals: defaultVisuals
};

export const BuiltInTrajectoryFormats = [
    ['mmcif', MmcifProvider] as const,
    ['cifCore', CifCoreProvider] as const,
    ['pdb', PdbProvider] as const,
    ['pdbqt', PdbqtProvider] as const,
    ['gro', GroProvider] as const,
    ['xyz', XyzProvider] as const,
    ['lammps_data', LammpsDataProvider] as const,
    ['lammps_traj_data', LammpsTrajectoryDataProvider] as const,
    ['mol', MolProvider] as const,
    ['sdf', SdfProvider] as const,
    ['mol2', Mol2Provider] as const,
] as const;

export type BuiltInTrajectoryFormat = (typeof BuiltInTrajectoryFormats)[number][0]