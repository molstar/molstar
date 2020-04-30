/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { CustomStructureProperty } from '../../mol-model-props/common/custom-structure-property';
import { Structure } from '../../mol-model/structure';
import { CustomProperty } from '../../mol-model-props/common/custom-property';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { Color } from '../../mol-util/color';
import { CustomPropertyDescriptor } from '../../mol-model/custom-property';

export type CellPackInfoValue = {
    packingsCount: number
    packingIndex: number
    colors?: Color[]
}

const CellPackInfoParams = {
    info: PD.Value<CellPackInfoValue>({ packingsCount: 1, packingIndex: 0, colors: undefined }, { isHidden: true })
};
type CellPackInfoParams = PD.Values<typeof CellPackInfoParams>

export const CellPackInfoProvider: CustomStructureProperty.Provider<typeof CellPackInfoParams, CellPackInfoValue> = CustomStructureProperty.createProvider({
    label: 'CellPack Info',
    descriptor: CustomPropertyDescriptor({ name: 'cellpack-info' }),
    type: 'root',
    defaultParams: CellPackInfoParams,
    getParams: (data: Structure) => CellPackInfoParams,
    isApplicable: (data: Structure) => true,
    obtain: async (ctx: CustomProperty.Context, data: Structure, props: CellPackInfoParams) => {
        return {
            value: { ...CellPackInfoParams.info.defaultValue, ...props.info }
        };
    }
});