/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { DataFormatProvider } from './provider';
import { PluginContext } from '../../../mol-plugin/context';
import { FileInfo } from '../../../mol-util/file-info';
import { StateTransforms } from '../../transforms';

export const MmcifFormatProvider = DataFormatProvider({
    id: 'mmcif',
    display: { name: 'mmCIF', group: 'Molecule' },
    extensions: { text: ['cif', 'mmcif', 'mcif'], binary: ['bcif'] }
})({
    isApplicable(plugin: PluginContext, data: string | Uint8Array, info?: FileInfo): boolean {
        // TODO: check CIF variants
        return true;
    },
    async apply({ plugin, state }, data) {
        const dictionary = state.build().to(data).apply(StateTransforms.Data.ParseCif, void 0, { state: { isGhost: true } });
        const trajectory = dictionary.apply(StateTransforms.Model.TrajectoryFromMmCif);
        await plugin.runTask(state.updateTree(trajectory));
        return { dictionary: dictionary.selector, trajectory: trajectory.selector };
    }
});

export const PdbFormatProvider = DataFormatProvider({
    id: 'pdb',
    display: { name: 'PDB', group: 'Molecule' },
    extensions: { text: ['pdb', 'ent'] }
})({
    async apply({ plugin, state }, data) {
        const trajectory = state.build().to(data).apply(StateTransforms.Model.TrajectoryFromPDB);
        await plugin.runTask(state.updateTree(trajectory));
        return { trajectory: trajectory.selector };
    }
});

export const BuiltInDataFormats = {
    'mmcif': MmcifFormatProvider,
    'pdb': PdbFormatProvider
}
export type BuiltInDataFormats = typeof BuiltInDataFormats

// export const TrajectoryFormatProviders = {
//     'mmcif': MmcifFormatProvider,
//     'pdb': PdbFormatProvider
// }