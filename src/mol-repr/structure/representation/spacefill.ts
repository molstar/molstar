/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { getElementSphereVisual, ElementSphereParams } from '../visual/element-sphere';
import { UnitsRepresentation } from '../units-representation';
import { ParamDefinition as PD } from 'mol-util/param-definition';
import { StructureRepresentation, StructureRepresentationProvider } from '../representation';
import { RepresentationParamsGetter, RepresentationContext, Representation } from 'mol-repr/representation';
import { ThemeRegistryContext } from 'mol-theme/theme';
import { Structure } from 'mol-model/structure';
import { UnitKind, UnitKindOptions } from '../visual/util/common';

const SpacefillVisuals = {
    'element-sphere': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, ElementSphereParams>) => UnitsRepresentation('Sphere mesh', ctx, getParams, getElementSphereVisual(ctx.webgl)),
}

export const SpacefillParams = {
    ...ElementSphereParams,
    unitKinds: PD.MultiSelect<UnitKind>(['atomic', 'spheres'], UnitKindOptions),
}
export type SpacefillParams = typeof SpacefillParams
export function getSpacefillParams(ctx: ThemeRegistryContext, structure: Structure) {
    return PD.clone(SpacefillParams)
}

export type SpacefillRepresentation = StructureRepresentation<SpacefillParams>
export function SpacefillRepresentation(ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, SpacefillParams>): SpacefillRepresentation {
    return Representation.createMulti('Spacefill', ctx, getParams, SpacefillVisuals as unknown as Representation.Def<Structure, SpacefillParams>)
}

export const SpacefillRepresentationProvider: StructureRepresentationProvider<SpacefillParams> = {
    label: 'Spacefill',
    description: 'Displays atomic/coarse elements as spheres.',
    factory: SpacefillRepresentation,
    getParams: getSpacefillParams,
    defaultValues: PD.getDefaultValues(SpacefillParams),
    defaultColorTheme: 'element-symbol',
    defaultSizeTheme: 'physical'
}