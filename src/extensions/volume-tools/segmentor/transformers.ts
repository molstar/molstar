/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Tadej Satler <tadej.satler@gmail.com>
 */

import { Volume } from '../../../mol-model/volume';
import { PluginStateObject as SO, PluginStateTransform } from '../../../mol-plugin-state/objects';
import { StateTransformer } from '../../../mol-state';
import { Task } from '../../../mol-task';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { buildCroppedMaskVolume, computeBodyMask } from './internal/mask-compute';
import { BodyLabels } from './labels';
import { MaxBodyId } from './types';

export { BodyMaskFromLabels };

export const BodyMaskFromLabelsTag = 'volume-segmentor-mask';

type BodyMaskFromLabels = typeof BodyMaskFromLabels
const BodyMaskFromLabels = PluginStateTransform.BuiltIn({
    name: 'volume-body-mask-from-labels',
    display: { name: 'Body Mask', description: 'Soft-edged mask of one body, derived from the voxel labels attached to the source volume' },
    from: SO.Volume.Data,
    to: SO.Volume.Data,
    params: {
        bodyId: PD.Numeric(1, { min: 1, max: MaxBodyId, step: 1 }, { isHidden: true }),
        extend: PD.Numeric(6, { min: 0, max: 30, step: 1 }, { description: 'Extend (dilate) the body by this many voxels' }),
        softEdge: PD.Numeric(2, { min: 0, max: 30, step: 1 }, { description: 'Width of the cosine soft edge in voxels' }),
        pruneBelowThreshold: PD.Boolean(true, { isHidden: true }),
        threshold: PD.Value<Volume.IsoValue>(Volume.IsoValue.relative(1), { isHidden: true }),
        /** Mirrors `LabelStore.version`; bumping it recomputes the mask after labels change. */
        version: PD.Numeric(0, {}, { isHidden: true }),
    },
})({
    apply({ a, params }) {
        return Task.create('Body Mask', async ctx => {
            const store = BodyLabels.get(a.data);
            if (!store) throw new Error('No body labels are attached to the source volume.');

            const body = BodyLabels.getBody(store, params.bodyId);
            const label = body?.name ?? `Body ${params.bodyId}`;
            const thresholdAbs = Volume.IsoValue.toAbsolute(params.threshold, a.data.grid.stats).absoluteValue;

            const result = await computeBodyMask(a.data, store.labels, params.bodyId, params, thresholdAbs, ctx);
            if (!result) throw new Error(`${label} has no voxels.`);

            const volume = buildCroppedMaskVolume(a.data, result, label);
            return new SO.Volume.Data(volume, {
                label,
                description: `${result.voxelCount.toLocaleString()} voxels, extend ${params.extend}, soft edge ${params.softEdge}`,
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
