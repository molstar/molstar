/**
 * Copyright (c) 2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 */

import { CustomModelProperty } from '../../../mol-model-props/common/custom-model-property';
import { CustomProperty } from '../../../mol-model-props/common/custom-property';
import { CustomPropertyDescriptor } from '../../../mol-model/custom-property';
import { Model, Structure } from '../../../mol-model/structure';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';


/** Parameter definition for custom model property "Is MVS" */
export type IsMVSModelParams = typeof IsMVSModelParams
export const IsMVSModelParams = {
    isMvs: PD.Boolean(false, { description: 'Flag this model as managed by MolViewSpec and enable MolViewSpec features' }),
};

/** Parameter values for custom model property "Is MVS" */
export type IsMVSModelProps = PD.Values<IsMVSModelParams>

/** Provider for custom model property "Is MVS" */
export const IsMVSModelProvider: CustomModelProperty.Provider<IsMVSModelParams, {}> = CustomModelProperty.createProvider({
    label: 'MVS',
    descriptor: CustomPropertyDescriptor({
        name: 'mvs-is-mvs-model',
    }),
    type: 'static',
    defaultParams: IsMVSModelParams,
    getParams: (data: Model) => IsMVSModelParams,
    isApplicable: (data: Model) => true,
    obtain: async (ctx: CustomProperty.Context, data: Model, props: Partial<IsMVSModelProps>) => ({ value: {} }),
});

/** Decide if the model is flagged as managed by MolViewSpec */
export function isMVSModel(model: Model): boolean {
    return !!IsMVSModelProvider.props(model)?.isMvs;
}
/** Decide if the structure is flagged as managed by MolViewSpec */
export function isMVSStructure(structure: Structure): boolean {
    return structure.models.some(isMVSModel);
}
