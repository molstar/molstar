/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { CIF } from '../../mol-io/reader/cif';
import { Vec3 } from '../../mol-math/linear-algebra';
import { volumeFromCcp4 } from '../../mol-model-formats/volume/ccp4';
import { volumeFromDensityServerData } from '../../mol-model-formats/volume/density-server';
import { volumeFromDsn6 } from '../../mol-model-formats/volume/dsn6';
import { Task } from '../../mol-task';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { PluginStateObject as SO, PluginStateTransform } from '../objects';

export { VolumeFromCcp4 };
export { VolumeFromDsn6 };
export { VolumeFromDensityServerCif };
type VolumeFromCcp4 = typeof VolumeFromCcp4
const VolumeFromCcp4 = PluginStateTransform.BuiltIn({
    name: 'volume-from-ccp4',
    display: { name: 'Volume from CCP4/MRC/MAP', description: 'Create Volume from CCP4/MRC/MAP data' },
    from: SO.Format.Ccp4,
    to: SO.Volume.Data,
    params(a) {
        return {
            voxelSize: PD.Vec3(Vec3.create(1, 1, 1)),
            offset: PD.Vec3(Vec3.create(0, 0, 0))
        };
    }
})({
    apply({ a, params }) {
        return Task.create('Create volume from CCP4/MRC/MAP', async ctx => {
            const volume = await volumeFromCcp4(a.data, params).runInContext(ctx);
            const props = { label: 'Volume' };
            return new SO.Volume.Data(volume, props);
        });
    }
});

type VolumeFromDsn6 = typeof VolumeFromDsn6
const VolumeFromDsn6 = PluginStateTransform.BuiltIn({
    name: 'volume-from-dsn6',
    display: { name: 'Volume from DSN6/BRIX', description: 'Create Volume from DSN6/BRIX data' },
    from: SO.Format.Dsn6,
    to: SO.Volume.Data,
    params(a) {
        return {
            voxelSize: PD.Vec3(Vec3.create(1, 1, 1))
        };
    }
})({
    apply({ a, params }) {
        return Task.create('Create volume from DSN6/BRIX', async ctx => {
            const volume = await volumeFromDsn6(a.data, params).runInContext(ctx);
            const props = { label: 'Volume' };
            return new SO.Volume.Data(volume, props);
        });
    }
});

type VolumeFromDensityServerCif = typeof VolumeFromDensityServerCif
const VolumeFromDensityServerCif = PluginStateTransform.BuiltIn({
    name: 'volume-from-density-server-cif',
    display: { name: 'Volume from density-server CIF', description: 'Identify and create all separate models in the specified CIF data block' },
    from: SO.Format.Cif,
    to: SO.Volume.Data,
    params(a) {
        if (!a) {
            return {
                blockHeader: PD.Optional(PD.Text(void 0, { description: 'Header of the block to parse. If none is specifed, the 1st data block in the file is used.' }))
            };
        }
        const blocks = a.data.blocks.slice(1); // zero block contains query meta-data
        return {
            blockHeader: PD.Optional(PD.Select(blocks[0] && blocks[0].header, blocks.map(b => [b.header, b.header] as [string, string]), { description: 'Header of the block to parse' }))
        };
    }
})({
    isApplicable: a => a.data.blocks.length > 0,
    apply({ a, params }) {
        return Task.create('Parse density-server CIF', async ctx => {
            const header = params.blockHeader || a.data.blocks[1].header; // zero block contains query meta-data
            const block = a.data.blocks.find(b => b.header === header);
            if (!block) throw new Error(`Data block '${[header]}' not found.`);
            const densityServerCif = CIF.schema.densityServer(block);
            const volume = await volumeFromDensityServerData(densityServerCif).runInContext(ctx);
            const props = { label: densityServerCif.volume_data_3d_info.name.value(0), description: `${densityServerCif.volume_data_3d_info.name.value(0)}` };
            return new SO.Volume.Data(volume, props);
        });
    }
});