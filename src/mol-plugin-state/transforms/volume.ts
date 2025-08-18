/**
 * Copyright (c) 2018-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Yakov Pechersky <ffxen158@gmail.com>
 */

import { CIF } from '../../mol-io/reader/cif';
import { Mat4, Vec3 } from '../../mol-math/linear-algebra';
import { volumeFromCcp4 } from '../../mol-model-formats/volume/ccp4';
import { volumeFromDensityServerData } from '../../mol-model-formats/volume/density-server';
import { volumeFromDsn6 } from '../../mol-model-formats/volume/dsn6';
import { Task } from '../../mol-task';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { PluginStateObject as SO, PluginStateTransform } from '../objects';
import { volumeFromCube } from '../../mol-model-formats/volume/cube';
import { volumeFromDx } from '../../mol-model-formats/volume/dx';
import { Grid, Volume } from '../../mol-model/volume';
import { PluginContext } from '../../mol-plugin/context';
import { StateSelection } from '../../mol-state';
import { volumeFromSegmentationData } from '../../mol-model-formats/volume/segmentation';
import { getTransformFromParams, TransformParam, transformParamsNeedCentroid } from './helpers';

export { VolumeFromCcp4 };
export { VolumeFromDsn6 };
export { VolumeFromCube };
export { VolumeFromDx };
export { AssignColorVolume };
export { VolumeFromDensityServerCif };
export { VolumeFromSegmentationCif };

type VolumeFromCcp4 = typeof VolumeFromCcp4
const VolumeFromCcp4 = PluginStateTransform.BuiltIn({
    name: 'volume-from-ccp4',
    display: { name: 'Volume from CCP4/MRC/MAP', description: 'Create Volume from CCP4/MRC/MAP data' },
    from: SO.Format.Ccp4,
    to: SO.Volume.Data,
    params(a) {
        return {
            voxelSize: PD.Vec3(Vec3.create(1, 1, 1)),
            offset: PD.Vec3(Vec3.create(0, 0, 0)),
            entryId: PD.Text(''),
        };
    }
})({
    apply({ a, params }) {
        return Task.create('Create volume from CCP4/MRC/MAP', async ctx => {
            const volume = await volumeFromCcp4(a.data, { ...params, label: a.data.name || a.label }).runInContext(ctx);
            const props = { label: volume.label || 'Volume', description: `Volume ${a.data.header.NX}\u00D7${a.data.header.NX}\u00D7${a.data.header.NX}` };
            return new SO.Volume.Data(volume, props);
        });
    },
    dispose({ b }) {
        b?.data.customProperties.dispose();
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
            voxelSize: PD.Vec3(Vec3.create(1, 1, 1)),
            entryId: PD.Text(''),
        };
    }
})({
    apply({ a, params }) {
        return Task.create('Create volume from DSN6/BRIX', async ctx => {
            const volume = await volumeFromDsn6(a.data, { ...params, label: a.data.name || a.label }).runInContext(ctx);
            const props = { label: volume.label || 'Volume', description: `Volume ${a.data.header.xExtent}\u00D7${a.data.header.yExtent}\u00D7${a.data.header.zExtent}` };
            return new SO.Volume.Data(volume, props);
        });
    },
    dispose({ b }) {
        b?.data.customProperties.dispose();
    }
});

type VolumeFromCube = typeof VolumeFromCube
const VolumeFromCube = PluginStateTransform.BuiltIn({
    name: 'volume-from-cube',
    display: { name: 'Volume from Cube', description: 'Create Volume from Cube data' },
    from: SO.Format.Cube,
    to: SO.Volume.Data,
    params(a) {
        const dataIndex = a ? PD.Select(0, a.data.header.dataSetIds.map((id, i) => [i, `${id}`] as const)) : PD.Numeric(0);
        return {
            dataIndex,
            entryId: PD.Text(''),
        };
    }
})({
    apply({ a, params }) {
        return Task.create('Create volume from Cube', async ctx => {
            const volume = await volumeFromCube(a.data, { ...params, label: a.data.name || a.label }).runInContext(ctx);
            const props = { label: volume.label || 'Volume', description: `Volume ${a.data.header.dim[0]}\u00D7${a.data.header.dim[1]}\u00D7${a.data.header.dim[2]}` };
            return new SO.Volume.Data(volume, props);
        });
    },
    dispose({ b }) {
        b?.data.customProperties.dispose();
    }
});

type VolumeFromDx = typeof VolumeFromDx
const VolumeFromDx = PluginStateTransform.BuiltIn({
    name: 'volume-from-dx',
    display: { name: 'Parse DX', description: 'Create volume from DX data.' },
    from: SO.Format.Dx,
    to: SO.Volume.Data
})({
    apply({ a }) {
        return Task.create('Parse DX', async ctx => {
            const volume = await volumeFromDx(a.data, { label: a.data.name || a.label }).runInContext(ctx);
            const props = { label: volume.label || 'Volume', description: `Volume ${a.data.header.dim[0]}\u00D7${a.data.header.dim[1]}\u00D7${a.data.header.dim[2]}` };
            return new SO.Volume.Data(volume, props);
        });
    },
    dispose({ b }) {
        b?.data.customProperties.dispose();
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
                blockHeader: PD.Optional(PD.Text(void 0, { description: 'Header of the block to parse. If none is specifed, the 1st data block in the file is used.' })),
                entryId: PD.Text(''),
            };
        }
        const blocks = a.data.blocks.slice(1); // zero block contains query meta-data
        return {
            blockHeader: PD.Optional(PD.Select(blocks[0] && blocks[0].header, blocks.map(b => [b.header, b.header] as [string, string]), { description: 'Header of the block to parse' })),
            entryId: PD.Text(''),
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
            const volume = await volumeFromDensityServerData(densityServerCif, { entryId: params.entryId }).runInContext(ctx);
            const [x, y, z] = volume.grid.cells.space.dimensions;
            const props = { label: params.entryId ?? densityServerCif.volume_data_3d_info.name.value(0), description: `Volume ${x}\u00D7${y}\u00D7${z}` };
            return new SO.Volume.Data(volume, props);
        });
    },
    dispose({ b }) {
        b?.data.customProperties.dispose();
    }
});

type VolumeFromSegmentationCif = typeof VolumeFromSegmentationCif
const VolumeFromSegmentationCif = PluginStateTransform.BuiltIn({
    name: 'volume-from-segmentation-cif',
    display: { name: 'Volume from Segmentation CIF' },
    from: SO.Format.Cif,
    to: SO.Volume.Data,
    params(a) {
        const blocks = a?.data.blocks.slice(1);
        const blockHeaderParam = blocks ?
            PD.Optional(PD.Select(blocks[0] && blocks[0].header, blocks.map(b => [b.header, b.header] as [string, string]), { description: 'Header of the block to parse' }))
            : PD.Optional(PD.Text(void 0, { description: 'Header of the block to parse. If none is specifed, the 1st data block in the file is used.' }));
        return {
            blockHeader: blockHeaderParam,
            segmentLabels: PD.ObjectList({ id: PD.Numeric(-1), label: PD.Text('') }, s => `${s.id} = ${s.label}`, { description: 'Mapping of segment IDs to segment labels' }),
            ownerId: PD.Text('', { isHidden: true, description: 'Reference to the object which manages this volume' }),
        };
    }
})({
    isApplicable: a => a.data.blocks.length > 0,
    apply({ a, params }) {
        return Task.create('Parse segmentation CIF', async ctx => {
            const header = params.blockHeader || a.data.blocks[1].header; // zero block contains query meta-data
            const block = a.data.blocks.find(b => b.header === header);
            if (!block) throw new Error(`Data block '${[header]}' not found.`);
            const segmentationCif = CIF.schema.segmentation(block);
            const segmentLabels: { [id: number]: string } = {};
            for (const segment of params.segmentLabels) segmentLabels[segment.id] = segment.label;
            const volume = await volumeFromSegmentationData(segmentationCif, { segmentLabels, ownerId: params.ownerId }).runInContext(ctx);
            const [x, y, z] = volume.grid.cells.space.dimensions;
            const props = { label: segmentationCif.volume_data_3d_info.name.value(0), description: `Segmentation ${x}\u00D7${y}\u00D7${z}` };
            return new SO.Volume.Data(volume, props);
        });
    },
    dispose({ b }) {
        b?.data.customProperties.dispose();
    }
});

type AssignColorVolume = typeof AssignColorVolume
const AssignColorVolume = PluginStateTransform.BuiltIn({
    name: 'assign-color-volume',
    display: { name: 'Assign Color Volume', description: 'Assigns another volume to be available for coloring.' },
    from: SO.Volume.Data,
    to: SO.Volume.Data,
    isDecorator: true,
    params(a, plugin: PluginContext) {
        if (!a) return { ref: PD.Text() };
        const cells = plugin.state.data.select(StateSelection.Generators.root.subtree().ofType(SO.Volume.Data).filter(cell => !!cell.obj && !cell.obj?.data.colorVolume && cell.obj !== a));
        if (cells.length === 0) return { ref: PD.Text('', { isHidden: true }) };
        return { ref: PD.Select(cells[0].transform.ref, cells.map(c => [c.transform.ref, c.obj!.label])) };
    }
})({
    apply({ a, params, dependencies }) {
        return Task.create('Assign Color Volume', async ctx => {
            if (!dependencies || !dependencies[params.ref]) {
                throw new Error('Dependency not available.');
            }
            const colorVolume = dependencies[params.ref].data as Volume;
            const volume: Volume = {
                ...a.data,
                colorVolume
            };
            const props = { label: a.label, description: 'Volume + Colors' };
            return new SO.Volume.Data(volume, props);
        });
    }
});

export type VolumeTransform = typeof VolumeTransform;
export const VolumeTransform = PluginStateTransform.BuiltIn({
    name: 'volume-transform',
    display: { name: 'Transform Volume' },
    isDecorator: true,
    from: SO.Volume.Data,
    to: SO.Volume.Data,
    params: {
        transform: TransformParam,
    },
})({
    canAutoUpdate() {
        return false;
    },
    apply({ a, params }) {
        // similar to StateTransforms.Model.TransformStructureConformation;
        const center = transformParamsNeedCentroid(params.transform) ? Grid.getBoundingSphere(a.data.grid).center : Vec3.unit;
        const transform = getTransformFromParams(params.transform, center);
        const gridTransform = {
            kind: 'matrix' as const,
            matrix: Mat4.mul(Mat4(), transform, Grid.getGridToCartesianTransform(a.data.grid)),
        };
        return new SO.Volume.Data({
            ...a.data,
            grid: {
                ...a.data.grid,
                transform: gridTransform,
            },
        }, {
            label: a.label,
            description: `${a.description} [Transformed]`,
        });
    },
});

export type VolumeInstances = typeof VolumeInstances;
export const VolumeInstances = PluginStateTransform.BuiltIn({
    name: 'volume-instances',
    display: { name: 'Volume Instances' },
    isDecorator: true,
    from: SO.Volume.Data,
    to: SO.Volume.Data,
    params: {
        transforms: PD.ObjectList({ transform: TransformParam }, () => 'Transform')
    },
})({
    canAutoUpdate() {
        return true;
    },
    apply({ a, params }) {
        const center = params.transforms.some(t => transformParamsNeedCentroid(t.transform)) ? Grid.getBoundingSphere(a.data.grid).center : Vec3.unit;
        const instances = params.transforms.map(t => ({ transform: getTransformFromParams(t.transform, center) }));
        if (!instances.length) {
            return a;
        }
        return new SO.Volume.Data({
            ...a.data,
            instances,
        }, {
            label: a.label,
            description: `${a.description} [Instanced]`,
        });
    },
});
