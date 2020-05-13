/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Structure, Model } from '../../../mol-model/structure';
import { Representation, RepresentationContext, RepresentationParamsGetter } from '../../../mol-repr/representation';
import { ThemeRegistryContext } from '../../../mol-theme/theme';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { ComplexRepresentation } from '../complex-representation';
import { StructureRepresentation, StructureRepresentationProvider, StructureRepresentationStateBuilder } from '../representation';
import { CarbohydrateLinkParams, CarbohydrateLinkVisual } from '../visual/carbohydrate-link-cylinder';
import { CarbohydrateSymbolParams, CarbohydrateSymbolVisual } from '../visual/carbohydrate-symbol-mesh';
import { CarbohydrateTerminalLinkParams, CarbohydrateTerminalLinkVisual } from '../visual/carbohydrate-terminal-link-cylinder';

const CarbohydrateVisuals = {
    'carbohydrate-symbol': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, CarbohydrateSymbolParams>) => ComplexRepresentation('Carbohydrate symbol mesh', ctx, getParams, CarbohydrateSymbolVisual),
    'carbohydrate-link': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, CarbohydrateLinkParams>) => ComplexRepresentation('Carbohydrate link cylinder', ctx, getParams, CarbohydrateLinkVisual),
    'carbohydrate-terminal-link': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, CarbohydrateTerminalLinkParams>) => ComplexRepresentation('Carbohydrate terminal link cylinder', ctx, getParams, CarbohydrateTerminalLinkVisual),
};

export const CarbohydrateParams = {
    ...CarbohydrateSymbolParams,
    ...CarbohydrateLinkParams,
    ...CarbohydrateTerminalLinkParams,
    visuals: PD.MultiSelect(['carbohydrate-symbol', 'carbohydrate-link', 'carbohydrate-terminal-link'], PD.objectToOptions(CarbohydrateVisuals)),
};
export type CarbohydrateParams = typeof CarbohydrateParams
export function getCarbohydrateParams(ctx: ThemeRegistryContext, structure: Structure) {
    return PD.clone(CarbohydrateParams);
}

export type CarbohydrateRepresentation = StructureRepresentation<CarbohydrateParams>
export function CarbohydrateRepresentation(ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, CarbohydrateParams>): CarbohydrateRepresentation {
    return Representation.createMulti('Carbohydrate', ctx, getParams, StructureRepresentationStateBuilder, CarbohydrateVisuals as unknown as Representation.Def<Structure, CarbohydrateParams>);
}

export const CarbohydrateRepresentationProvider = StructureRepresentationProvider({
    name: 'carbohydrate',
    label: 'Carbohydrate',
    description: 'Displays carbohydrate symbols (3D SNFG).',
    factory: CarbohydrateRepresentation,
    getParams: getCarbohydrateParams,
    defaultValues: PD.getDefaultValues(CarbohydrateParams),
    defaultColorTheme: { name: 'carbohydrate-symbol' },
    defaultSizeTheme: { name: 'uniform' },
    isApplicable: (structure: Structure) => {
        return structure.models.some(m => Model.hasCarbohydrate(m));
    }
});