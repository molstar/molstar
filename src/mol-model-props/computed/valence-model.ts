/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Structure } from '../../mol-model/structure';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { calcValenceModel, ValenceModel, ValenceModelParams as _ValenceModelParams } from './chemistry/valence-model';
import { CustomStructureProperty } from '../common/custom-structure-property';
import { CustomProperty } from '../common/custom-property';
import { CustomPropertyDescriptor } from '../../mol-model/custom-property';

export const ValenceModelParams = {
    ..._ValenceModelParams
};
export type ValenceModelParams = typeof ValenceModelParams
export type ValenceModelProps = PD.Values<ValenceModelParams>

export type ValenceModelValue = Map<number, ValenceModel>

export const ValenceModelProvider: CustomStructureProperty.Provider<ValenceModelParams, ValenceModelValue> = CustomStructureProperty.createProvider({
    label: 'Valence Model',
    descriptor: CustomPropertyDescriptor({
        name: 'molstar_computed_valence_model',
        // TODO `cifExport` and `symbol`
    }),
    type: 'local',
    defaultParams: ValenceModelParams,
    getParams: (data: Structure) => ValenceModelParams,
    isApplicable: (data: Structure) => true,
    obtain: async (ctx: CustomProperty.Context, data: Structure, props: Partial<ValenceModelProps>) => {
        const p = { ...PD.getDefaultValues(ValenceModelParams), ...props };
        return { value: await calcValenceModel(ctx.runtime, data, p) };
    }
});