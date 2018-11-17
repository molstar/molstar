/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { CarbohydrateSymbolVisual, CarbohydrateSymbolParams } from '../visual/carbohydrate-symbol-mesh';
import { CarbohydrateLinkVisual, CarbohydrateLinkParams } from '../visual/carbohydrate-link-cylinder';
import { CarbohydrateTerminalLinkParams, CarbohydrateTerminalLinkVisual } from '../visual/carbohydrate-terminal-link-cylinder';
import { ParamDefinition as PD } from 'mol-util/param-definition';
import { ComplexRepresentation } from '../complex-representation';
import { StructureRepresentation, StructureRepresentationProvider } from '../representation';
import { Representation, RepresentationParamsGetter } from 'mol-repr/representation';
import { ThemeRegistryContext } from 'mol-theme/theme';
import { Structure } from 'mol-model/structure';
import { BuiltInColorThemeOptions, getBuiltInColorThemeParams } from 'mol-theme/color';

const CarbohydrateVisuals = {
    'carbohydrate-symbol': (getParams: RepresentationParamsGetter<Structure, CarbohydrateSymbolParams>) => ComplexRepresentation('Carbohydrate symbol mesh', getParams, CarbohydrateSymbolVisual),
    'carbohydrate-link': (getParams: RepresentationParamsGetter<Structure, CarbohydrateLinkParams>) => ComplexRepresentation('Carbohydrate link cylinder', getParams, CarbohydrateLinkVisual),
    'carbohydrate-terminal-link': (getParams: RepresentationParamsGetter<Structure, CarbohydrateTerminalLinkParams>) => ComplexRepresentation('Carbohydrate terminal link cylinder', getParams, CarbohydrateTerminalLinkVisual),
}
type CarbohydrateVisualName = keyof typeof CarbohydrateVisuals
const CarbohydrateVisualOptions = Object.keys(CarbohydrateVisuals).map(name => [name, name] as [CarbohydrateVisualName, string])

export const CarbohydrateParams = {
    ...CarbohydrateSymbolParams,
    ...CarbohydrateLinkParams,
    ...CarbohydrateTerminalLinkParams,
    colorTheme: PD.Mapped('carbohydrate-symbol', BuiltInColorThemeOptions, getBuiltInColorThemeParams),
    visuals: PD.MultiSelect<CarbohydrateVisualName>(['carbohydrate-symbol', 'carbohydrate-link', 'carbohydrate-terminal-link'], CarbohydrateVisualOptions),
}
PD.getDefaultValues(CarbohydrateParams).colorTheme.name
export type CarbohydrateParams = typeof CarbohydrateParams
export function getCarbohydrateParams(ctx: ThemeRegistryContext, structure: Structure) {
    return PD.clone(CarbohydrateParams)
}

export type CarbohydrateRepresentation = StructureRepresentation<CarbohydrateParams>
export function CarbohydrateRepresentation(getParams: RepresentationParamsGetter<Structure, CarbohydrateParams>): CarbohydrateRepresentation {
    return Representation.createMulti('Carbohydrate', getParams, CarbohydrateVisuals as unknown as Representation.Def<Structure, CarbohydrateParams>)
}

export const CarbohydrateRepresentationProvider: StructureRepresentationProvider<CarbohydrateParams> = {
    label: 'Carbohydrate',
    description: 'Displays carbohydrate symbols (3D SNFG).',
    factory: CarbohydrateRepresentation,
    getParams: getCarbohydrateParams,
    defaultValues: PD.getDefaultValues(CarbohydrateParams)
}