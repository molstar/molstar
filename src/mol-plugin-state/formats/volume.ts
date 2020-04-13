/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { StateTransforms } from '../transforms';
import { DataFormatProvider, guessCifVariant } from './provider';
import { PluginContext } from '../../mol-plugin/context';
import { StateObjectSelector } from '../../mol-state';
import { PluginStateObject } from '../objects';
import { VolumeRepresentation3DHelpers } from '../transforms/representation';
import { ColorNames } from '../../mol-util/color/names';
import { VolumeIsoValue } from '../../mol-model/volume';

const Category = 'Volume';

async function defaultVisuals(plugin: PluginContext, data: { volume: StateObjectSelector<PluginStateObject.Volume.Data> }) {
    const visual = plugin.build().to(data.volume).apply(StateTransforms.Representation.VolumeRepresentation3D);
    return [await visual.commit()];
}

export const Ccp4Provider = DataFormatProvider({
    label: 'CCP4/MRC/BRIX',
    description: 'CCP4/MRC/BRIX',
    category: Category,
    binaryExtensions: ['ccp4', 'mrc', 'map'],
    parse: async (plugin, data) => {
        const format = plugin.build()
            .to(data)
            .apply(StateTransforms.Data.ParseCcp4, {}, { state: { isGhost: true } });

        const volume = format.apply(StateTransforms.Volume.VolumeFromCcp4);

        await format.commit({ revertOnError: true });

        return { format: format.selector, volume: volume.selector };
    },
    visuals: defaultVisuals
});

export const Dsn6Provider = DataFormatProvider({
    label: 'DSN6/BRIX',
    description: 'DSN6/BRIX',
    category: Category,
    binaryExtensions: ['dsn6', 'brix'],
    parse: async (plugin, data) => {
        const format = plugin.build()
            .to(data)
            .apply(StateTransforms.Data.ParseDsn6, {}, { state: { isGhost: true } });

        const volume = format.apply(StateTransforms.Volume.VolumeFromDsn6);

        await format.commit({ revertOnError: true });

        return { format: format.selector, volume: volume.selector };
    },
    visuals: defaultVisuals
});

export const DscifProvider = DataFormatProvider({
    label: 'DensityServer CIF',
    description: 'DensityServer CIF',
    category: Category,
    stringExtensions: ['cif'],
    binaryExtensions: ['bcif'],
    isApplicable: (info, data) => {
        return guessCifVariant(info, data) === 'dscif' ? true : false;
    },
    parse: async (plugin, data) => {
        const cifCell = await plugin.build().to(data).apply(StateTransforms.Data.ParseCif).commit();
        const b = plugin.build().to(cifCell);
        const blocks = cifCell.obj!.data.blocks.slice(1); // zero block contains query meta-data

        if (blocks.length !== 1 && blocks.length !== 2) throw new Error('unknown number of blocks');

        const volumes: StateObjectSelector<PluginStateObject.Volume.Data>[] = [];
        for (const block of blocks) {
            volumes.push(b.apply(StateTransforms.Volume.VolumeFromDensityServerCif, { blockHeader: block.header }).selector);
        }

        await b.commit();

        return { volumes };
    },
    visuals: async (plugin, data: { volumes: StateObjectSelector<PluginStateObject.Volume.Data>[] }) => {
        const { volumes } = data;
        const tree = plugin.build();
        const visuals: StateObjectSelector<PluginStateObject.Volume.Representation3D>[] = [];

        if (volumes.length > 0) {
            visuals[0] = tree
                .to(volumes[0])
                .apply(StateTransforms.Representation.VolumeRepresentation3D, VolumeRepresentation3DHelpers.getDefaultParamsStatic(plugin, 'isosurface', { isoValue: VolumeIsoValue.relative(1.5), alpha: 0.3 }, 'uniform', { value: ColorNames.teal }))
                .selector;
        }

        if (volumes.length > 1) {
            const posParams = VolumeRepresentation3DHelpers.getDefaultParamsStatic(plugin, 'isosurface', { isoValue: VolumeIsoValue.relative(3), alpha: 0.3 }, 'uniform', { value: ColorNames.green });
            const negParams = VolumeRepresentation3DHelpers.getDefaultParamsStatic(plugin, 'isosurface', { isoValue: VolumeIsoValue.relative(-3), alpha: 0.3 }, 'uniform', { value: ColorNames.red });
            visuals[visuals.length] = tree.to(volumes[1]).apply(StateTransforms.Representation.VolumeRepresentation3D, posParams).selector;
            visuals[visuals.length] = tree.to(volumes[1]).apply(StateTransforms.Representation.VolumeRepresentation3D, negParams).selector;
        }

        await tree.commit();

        return visuals;
    }
});

export const BuiltInVolumeFormats = [
    ['ccp4', Ccp4Provider] as const,
    ['dns6', Dsn6Provider] as const,
    ['dscif', DscifProvider] as const,
] as const;

export type BuildInVolumeFormat = (typeof BuiltInVolumeFormats)[number][0]