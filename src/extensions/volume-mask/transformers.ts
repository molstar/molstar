/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 * @author Tadej Satler <tadej.satler@gmail.com>
 */

import { Volume } from '../../mol-model/volume';
import { PluginStateObject as SO, PluginStateTransform } from '../../mol-plugin-state/objects';
import { Task } from '../../mol-task';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { StateTransformer } from '../../mol-state';
import { computeVolumeMask, computeStructureMask, computeSoftEdge, buildMaskVolume, dilate3D } from './internal/mask-compute';
import type { ViewMask, MaskSource } from './types';

export { MaskVolumeFromSource };

type MaskVolumeFromSource = typeof MaskVolumeFromSource;
const MaskVolumeFromSource = PluginStateTransform.BuiltIn({
    name: 'volume-mask-from-source',
    display: { name: 'Volume Mask', description: 'Binary voxel mask from 2D polygon projections, structure proximity, or both' },
    from: SO.Volume.Data,
    to: SO.Volume.Data,
    params: {
        viewMasks: PD.Value<ViewMask[]>([], { isHidden: true }),
        threshold: PD.Value<Volume.IsoValue>(Volume.IsoValue.relative(1), { isHidden: true }),
        dilation: PD.Numeric(6, { min: 0, max: 30, step: 1 }, { description: 'Dilation radius in voxels' }),
        maskSource: PD.Value<MaskSource>('volume', { isHidden: true }),
        atomPositions: PD.Value<Float32Array | undefined>(undefined, { isHidden: true }),
        proteinRadius: PD.Numeric(1, { min: 0.5, max: 20, step: 0.5 }, { description: 'Atom proximity radius in Å' }),
        softEdge: PD.Numeric(3, { min: 0, max: 20, step: 1 }, { description: 'Soft edge width in voxels' }),
        invertOutput: PD.Boolean(false, { isHidden: true }),
    },
})({
    apply({ a, params }) {
        return Task.create('Create Volume Mask', async ctx => {
            let maskData: Uint8Array;

            if (params.maskSource === 'volume') {
                maskData = await computeVolumeMask(a.data, { ...params, invertPolygons: params.invertOutput }, ctx);
            } else {
                // 'structure': proximity mask, optionally AND'd with polygon spatial filter
                await ctx.update({ message: 'Computing structure mask…' });
                maskData = computeStructureMask(a.data, params.atomPositions!, params.proteinRadius);

                if (params.viewMasks.length > 0) {
                    await ctx.update({ message: 'Applying polygon filter…' });
                    // When inverted: select voxels inside polygons but NOT near structure
                    const polyMask = await computeVolumeMask(a.data, { ...params, skipThreshold: true, dilation: 0 }, ctx);
                    if (params.invertOutput) {
                        for (let i = 0; i < maskData.length; i++) maskData[i] = (polyMask[i] && !maskData[i]) ? 1 : 0;
                    } else {
                        for (let i = 0; i < maskData.length; i++) maskData[i] = maskData[i] & polyMask[i];
                    }
                } else if (params.invertOutput) {
                    // No polygons: invert flips the structure mask
                    for (let i = 0; i < maskData.length; i++) maskData[i] = maskData[i] ? 0 : 1;
                }

                if (params.dilation > 0) {
                    await ctx.update({ message: 'Dilating mask…' });
                    const { cells: { space } } = a.data.grid;
                    const [nx, ny, nz] = space.dimensions as [number, number, number];
                    maskData = dilate3D(maskData, nx, ny, nz, params.dilation, space);
                }
            }

            const { cells: { space: sp } } = a.data.grid;
            const [nx, ny, nz] = sp.dimensions as [number, number, number];
            let count = 0;
            for (let i = 0; i < maskData.length; i++) if (maskData[i]) count++;

            let finalData: Uint8Array | Float32Array = maskData;
            if (params.softEdge > 0) {
                await ctx.update({ message: 'Applying soft edge…' });
                finalData = computeSoftEdge(maskData, nx, ny, nz, params.softEdge, sp);
            }

            const maskVolume = buildMaskVolume(a.data, finalData);
            return new SO.Volume.Data(maskVolume, {
                label: 'Mask',
                description: `${count.toLocaleString()} voxels (${nx}\u00D7${ny}\u00D7${nz})`,
            });
        });
    },
    update() {
        return StateTransformer.UpdateResult.Recreate;
    },
    dispose({ b }) {
        b?.data.customProperties.dispose();
    },
});
