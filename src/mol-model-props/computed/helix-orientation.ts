/**
 * Copyright (c) 2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { Model } from '../../mol-model/structure';
import { CustomPropertyDescriptor } from '../../mol-model/custom-property';
import { CustomModelProperty } from '../common/custom-model-property';
import { calcHelixOrientation, HelixOrientation } from './helix-orientation/helix-orientation';

export const HelixOrientationParams = { };
export type HelixOrientationParams = typeof HelixOrientationParams
export type HelixOrientationProps = PD.Values<HelixOrientationParams>

export type HelixOrientationValue = HelixOrientation;

export const HelixOrientationProvider: CustomModelProperty.Provider<HelixOrientationParams, HelixOrientationValue> = CustomModelProperty.createProvider({
    label: 'Helix Orientation',
    descriptor: CustomPropertyDescriptor({
        name: 'molstar_helix_orientation'
    }),
    type: 'dynamic',
    defaultParams: {},
    getParams: () => ({}),
    isApplicable: (data: Model) => true,
    obtain: async (ctx, data) => {
        return { value: calcHelixOrientation(data) };
    }
});