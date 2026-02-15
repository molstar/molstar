/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { CustomProperty } from '../common/custom-property';
import { CustomPropertyDescriptor } from '../../mol-model/custom-property';
import { CustomVolumeProperty } from '../common/custom-volume-property';
import { Volume } from '../../mol-model/volume/volume';
import { calculateBasicStreamlines, BasicStreamlineCalculationParams } from './streamlines/basic';
import { Streamlines } from './streamlines/shared';

export const StreamlinesParams = {
    type: PD.MappedStatic('basic', {
        'basic': PD.Group(BasicStreamlineCalculationParams, { isFlat: true }),
    })
};

export type StreamlinesParams = typeof StreamlinesParams
export type StreamlinesProps = PD.Values<StreamlinesParams>
export type StreamlinesValue = Streamlines

export const StreamlinesProvider: CustomVolumeProperty.Provider<StreamlinesParams, StreamlinesValue> = CustomVolumeProperty.createProvider({
    label: 'Streamlines',
    descriptor: CustomPropertyDescriptor({
        name: 'molstar_streamlines',
        // TODO `cifExport` and `symbol`
    }),
    defaultParams: StreamlinesParams,
    getParams: (data: Volume) => StreamlinesParams,
    isApplicable: (data: Volume) => !Volume.Segmentation.get(data),
    obtain: async (ctx: CustomProperty.Context, data: Volume, props: Partial<StreamlinesProps>) => {
        const p = { ...PD.getDefaultValues(StreamlinesParams), ...props };
        switch (p.type.name) {
            case 'basic': return { value: await calculateBasicStreamlines(ctx, data, p.type.params) };
        }
    }
});
