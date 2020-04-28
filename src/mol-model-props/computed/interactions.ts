/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Structure } from '../../mol-model/structure';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { computeInteractions, Interactions, InteractionsParams as _InteractionsParams } from './interactions/interactions';
import { CustomStructureProperty } from '../common/custom-structure-property';
import { CustomProperty } from '../common/custom-property';
import { CustomPropertyDescriptor } from '../../mol-model/custom-property';

export const InteractionsParams = {
    ..._InteractionsParams
};
export type InteractionsParams = typeof InteractionsParams
export type InteractionsProps = PD.Values<InteractionsParams>

export type InteractionsValue = Interactions

export const InteractionsProvider: CustomStructureProperty.Provider<InteractionsParams, InteractionsValue> = CustomStructureProperty.createProvider({
    label: 'Interactions',
    descriptor: CustomPropertyDescriptor({
        name: 'molstar_computed_interactions',
        // TODO `cifExport` and `symbol`
    }),
    type: 'local',
    defaultParams: InteractionsParams,
    getParams: (data: Structure) => InteractionsParams,
    isApplicable: (data: Structure) => true,
    obtain: async (ctx: CustomProperty.Context, data: Structure, props: Partial<InteractionsProps>) => {
        const p = { ...PD.getDefaultValues(InteractionsParams), ...props };
        return { value: await computeInteractions(ctx, data, p) };
    }
});