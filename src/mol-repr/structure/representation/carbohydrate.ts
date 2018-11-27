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
import { Representation, RepresentationParamsGetter, RepresentationContext } from 'mol-repr/representation';
import { ThemeRegistryContext } from 'mol-theme/theme';
import { Structure } from 'mol-model/structure';

const CarbohydrateVisuals = {
    'carbohydrate-symbol': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, CarbohydrateSymbolParams>) => ComplexRepresentation('Carbohydrate symbol mesh', ctx, getParams, CarbohydrateSymbolVisual),
    'carbohydrate-link': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, CarbohydrateLinkParams>) => ComplexRepresentation('Carbohydrate link cylinder', ctx, getParams, CarbohydrateLinkVisual),
    'carbohydrate-terminal-link': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, CarbohydrateTerminalLinkParams>) => ComplexRepresentation('Carbohydrate terminal link cylinder', ctx, getParams, CarbohydrateTerminalLinkVisual),
}
type CarbohydrateVisualName = keyof typeof CarbohydrateVisuals
const CarbohydrateVisualOptions = Object.keys(CarbohydrateVisuals).map(name => [name, name] as [CarbohydrateVisualName, string])

export const CarbohydrateParams = {
    ...CarbohydrateSymbolParams,
    ...CarbohydrateLinkParams,
    ...CarbohydrateTerminalLinkParams,
    visuals: PD.MultiSelect<CarbohydrateVisualName>(['carbohydrate-symbol', 'carbohydrate-link', 'carbohydrate-terminal-link'], CarbohydrateVisualOptions),
}
export type CarbohydrateParams = typeof CarbohydrateParams
export function getCarbohydrateParams(ctx: ThemeRegistryContext, structure: Structure) {
    return PD.clone(CarbohydrateParams)
}

export type CarbohydrateRepresentation = StructureRepresentation<CarbohydrateParams>
export function CarbohydrateRepresentation(ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, CarbohydrateParams>): CarbohydrateRepresentation {
    return Representation.createMulti('Carbohydrate', ctx, getParams, CarbohydrateVisuals as unknown as Representation.Def<Structure, CarbohydrateParams>)
}

export const CarbohydrateRepresentationProvider: StructureRepresentationProvider<CarbohydrateParams> = {
    label: 'Carbohydrate',
    description: 'Displays carbohydrate symbols (3D SNFG).',
    factory: CarbohydrateRepresentation,
    getParams: getCarbohydrateParams,
    defaultValues: PD.getDefaultValues(CarbohydrateParams)
}