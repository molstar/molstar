/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { CarbohydrateSymbolVisual, CarbohydrateSymbolParams } from '../visual/carbohydrate-symbol-mesh';
import { CarbohydrateLinkVisual, CarbohydrateLinkParams } from '../visual/carbohydrate-link-cylinder';
import { SizeThemeName, SizeThemeOptions } from 'mol-theme/size';
import { ParamDefinition as PD } from 'mol-util/param-definition';
import { ComplexRepresentation } from '../complex-representation';
import { StructureRepresentation } from '../index';
import { Representation } from 'mol-repr';

export const CarbohydrateParams = {
    ...CarbohydrateSymbolParams,
    ...CarbohydrateLinkParams,
    sizeTheme: PD.Select<SizeThemeName>('Size Theme', '', 'uniform', SizeThemeOptions),
    sizeValue: PD.Numeric('Size Value', '', 1, 0, 0.1, 20),
    sizeFactor: PD.Numeric('Size Factor', '', 1, 0, 10, 0.1),
}
export const DefaultCarbohydrateProps = PD.getDefaultValues(CarbohydrateParams)
export type CarbohydrateProps = typeof DefaultCarbohydrateProps

export type CarbohydrateRepresentation = StructureRepresentation<CarbohydrateProps>

export function CarbohydrateRepresentation(): CarbohydrateRepresentation {
    return Representation.createMulti('Carbohydrate', CarbohydrateParams, DefaultCarbohydrateProps, [
        ComplexRepresentation('Carbohydrate symbol mesh', CarbohydrateSymbolVisual),
        ComplexRepresentation('Carbohydrate link cylinder', CarbohydrateLinkVisual)
    ] as StructureRepresentation<CarbohydrateProps>[])
}