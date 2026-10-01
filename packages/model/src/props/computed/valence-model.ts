/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { Structure } from '@molstar/model/model/structure';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { calcValenceModel, type ValenceModel, ValenceModelParams as _ValenceModelParams } from './chemistry/valence-model.js';
import { CustomStructureProperty } from '../common/custom-structure-property.js';
import type { CustomProperty } from '../common/custom-property.js';
import { CustomPropertyDescriptor } from '@molstar/model/model/custom-property';

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