/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { Representation, RepresentationParamsGetter, RepresentationContext } from '../../representation';
import { ThemeRegistryContext } from '../../../mol-theme/theme';
import { Structure } from '../../../mol-model/structure';
import { UnitsRepresentation, StructureRepresentation, StructureRepresentationStateBuilder, StructureRepresentationProvider } from '../representation';
import { InteractionsIntraUnitParams, InteractionsIntraUnitVisual } from '../visual/interactions-intra-unit-cylinder';
import { UnitKindOptions, UnitKind } from '../visual/util/common';
import { InteractionsProvider } from '../../../mol-model-props/computed/interactions';

const InteractionsVisuals = {
    'intra-unit': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, InteractionsIntraUnitParams>) => UnitsRepresentation('Intra-unit interactions cylinder', ctx, getParams, InteractionsIntraUnitVisual),
}

export const InteractionsParams = {
    ...InteractionsIntraUnitParams,
    unitKinds: PD.MultiSelect<UnitKind>(['atomic'], UnitKindOptions),
    sizeFactor: PD.Numeric(0.3, { min: 0.01, max: 10, step: 0.01 }),
    sizeAspectRatio: PD.Numeric(2/3, { min: 0.01, max: 3, step: 0.01 }),
    visuals: PD.MultiSelect(['intra-unit'], PD.objectToOptions(InteractionsVisuals)),
}
export type InteractionsParams = typeof InteractionsParams
export function getInteractionParams(ctx: ThemeRegistryContext, structure: Structure) {
    return PD.clone(InteractionsParams)
}

export type InteractionRepresentation = StructureRepresentation<InteractionsParams>
export function InteractionRepresentation(ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, InteractionsParams>): InteractionRepresentation {
    return Representation.createMulti('Interactions', ctx, getParams, StructureRepresentationStateBuilder, InteractionsVisuals as unknown as Representation.Def<Structure, InteractionsParams>)
}

export const InteractionsRepresentationProvider: StructureRepresentationProvider<InteractionsParams> = {
    label: 'Non-covalent Interactions',
    description: 'Displays non-covalent interactions as dashed cylinders.',
    factory: InteractionRepresentation,
    getParams: getInteractionParams,
    defaultValues: PD.getDefaultValues(InteractionsParams),
    defaultColorTheme: 'interaction-type',
    defaultSizeTheme: 'uniform',
    isApplicable: (structure: Structure) => structure.elementCount > 0,
    ensureDependencies: (structure: Structure) => {
        return InteractionsProvider.attach(structure.root)
    }
}