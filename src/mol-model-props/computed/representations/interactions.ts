/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { Representation, RepresentationParamsGetter, RepresentationContext } from '../../../mol-repr/representation';
import { ThemeRegistryContext } from '../../../mol-theme/theme';
import { Structure } from '../../../mol-model/structure';
import { UnitsRepresentation, StructureRepresentation, StructureRepresentationStateBuilder, StructureRepresentationProvider, ComplexRepresentation } from '../../../mol-repr/structure/representation';
import { InteractionsIntraUnitParams, InteractionsIntraUnitVisual } from './interactions-intra-unit-cylinder';
import { InteractionsProvider } from '../interactions';
import { InteractionsInterUnitParams, InteractionsInterUnitVisual } from './interactions-inter-unit-cylinder';
import { CustomProperty } from '../../common/custom-property';
import { getUnitKindsParam } from '../../../mol-repr/structure/params';

const InteractionsVisuals = {
    'intra-unit': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, InteractionsIntraUnitParams>) => UnitsRepresentation('Intra-unit interactions cylinder', ctx, getParams, InteractionsIntraUnitVisual),
    'inter-unit': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, InteractionsInterUnitParams>) => ComplexRepresentation('Inter-unit interactions cylinder', ctx, getParams, InteractionsInterUnitVisual),
};

export const InteractionsParams = {
    ...InteractionsIntraUnitParams,
    ...InteractionsInterUnitParams,
    unitKinds: getUnitKindsParam(['atomic']),
    sizeFactor: PD.Numeric(0.2, { min: 0.01, max: 1, step: 0.01 }),
    visuals: PD.MultiSelect(['intra-unit', 'inter-unit'], PD.objectToOptions(InteractionsVisuals)),
};
export type InteractionsParams = typeof InteractionsParams
export function getInteractionParams(ctx: ThemeRegistryContext, structure: Structure) {
    return PD.clone(InteractionsParams);
}

export type InteractionRepresentation = StructureRepresentation<InteractionsParams>
export function InteractionRepresentation(ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, InteractionsParams>): InteractionRepresentation {
    return Representation.createMulti('Interactions', ctx, getParams, StructureRepresentationStateBuilder, InteractionsVisuals as unknown as Representation.Def<Structure, InteractionsParams>);
}

export const InteractionsRepresentationProvider = StructureRepresentationProvider({
    name: 'interactions',
    label: 'Non-covalent Interactions',
    description: 'Displays non-covalent interactions as dashed cylinders.',
    factory: InteractionRepresentation,
    getParams: getInteractionParams,
    defaultValues: PD.getDefaultValues(InteractionsParams),
    defaultColorTheme: { name: 'interaction-type' },
    defaultSizeTheme: { name: 'uniform' },
    isApplicable: (structure: Structure) => structure.elementCount > 0,
    ensureCustomProperties: {
        attach: (ctx: CustomProperty.Context, structure: Structure) => InteractionsProvider.attach(ctx, structure, void 0, true),
        detach: (data) => InteractionsProvider.ref(data, false)
    }
});