/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { CrossLinkRestraintVisual, CrossLinkRestraintParams } from '../visual/cross-link-restraint-cylinder';
import { ParamDefinition as PD } from 'mol-util/param-definition';
import { ComplexRepresentation } from '../complex-representation';
import { StructureRepresentation, StructureRepresentationProvider, StructureRepresentationStateBuilder } from '../representation';
import { Representation, RepresentationContext, RepresentationParamsGetter } from 'mol-repr/representation';
import { ThemeRegistryContext } from 'mol-theme/theme';
import { Structure } from 'mol-model/structure';
import { UnitKind, UnitKindOptions } from '../visual/util/common';

const DistanceRestraintVisuals = {
    'distance-restraint': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, CrossLinkRestraintParams>) => ComplexRepresentation('Cross-link restraint', ctx, getParams, CrossLinkRestraintVisual),
}

export const DistanceRestraintParams = {
    ...CrossLinkRestraintParams,
    unitKinds: PD.MultiSelect<UnitKind>(['atomic', 'spheres'], UnitKindOptions),
}
export type DistanceRestraintParams = typeof DistanceRestraintParams
export function getDistanceRestraintParams(ctx: ThemeRegistryContext, structure: Structure) {
    return PD.clone(DistanceRestraintParams)
}

export type DistanceRestraintRepresentation = StructureRepresentation<DistanceRestraintParams>
export function DistanceRestraintRepresentation(ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, DistanceRestraintParams>): DistanceRestraintRepresentation {
    return Representation.createMulti('DistanceRestraint', ctx, getParams, StructureRepresentationStateBuilder, DistanceRestraintVisuals as unknown as Representation.Def<Structure, DistanceRestraintParams>)
}

export const DistanceRestraintRepresentationProvider: StructureRepresentationProvider<DistanceRestraintParams> = {
    label: 'Distance Restraint',
    description: 'Displays cross-link distance restraints.',
    factory: DistanceRestraintRepresentation,
    getParams: getDistanceRestraintParams,
    defaultValues: PD.getDefaultValues(DistanceRestraintParams),
    defaultColorTheme: 'cross-link',
    defaultSizeTheme: 'uniform'
}