/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { getElementSphereVisual, ElementSphereParams } from '../visual/element-sphere';
import { IntraUnitLinkVisual, IntraUnitLinkParams } from '../visual/intra-unit-link-cylinder';
import { InterUnitLinkVisual, InterUnitLinkParams } from '../visual/inter-unit-link-cylinder';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { UnitsRepresentation } from '../units-representation';
import { ComplexRepresentation } from '../complex-representation';
import { StructureRepresentation, StructureRepresentationProvider, StructureRepresentationStateBuilder } from '../representation';
import { Representation, RepresentationParamsGetter, RepresentationContext } from '../../../mol-repr/representation';
import { ThemeRegistryContext } from '../../../mol-theme/theme';
import { Structure } from '../../../mol-model/structure';
import { UnitKind, UnitKindOptions } from '../visual/util/common';

const BallAndStickVisuals = {
    'element-sphere': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, ElementSphereParams>) => UnitsRepresentation('Element sphere mesh', ctx, getParams, getElementSphereVisual(ctx.webgl)),
    'intra-link': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, IntraUnitLinkParams>) => UnitsRepresentation('Intra-unit link cylinder', ctx, getParams, IntraUnitLinkVisual),
    'inter-link': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, InterUnitLinkParams>) => ComplexRepresentation('Inter-unit link cylinder', ctx, getParams, InterUnitLinkVisual),
}
type BallAndStickVisualName = keyof typeof BallAndStickVisuals
const BallAndStickVisualOptions = Object.keys(BallAndStickVisuals).map(name => [name, name] as [BallAndStickVisualName, string])

export const BallAndStickParams = {
    ...ElementSphereParams,
    ...IntraUnitLinkParams,
    ...InterUnitLinkParams,
    unitKinds: PD.MultiSelect<UnitKind>(['atomic'], UnitKindOptions),
    sizeFactor: PD.Numeric(0.15, { min: 0.01, max: 10, step: 0.01 }),
    sizeAspectRatio: PD.Numeric(2/3, { min: 0.01, max: 3, step: 0.01 }),
    visuals: PD.MultiSelect<BallAndStickVisualName>(['element-sphere', 'intra-link', 'inter-link'], BallAndStickVisualOptions),
}
export type BallAndStickParams = typeof BallAndStickParams
export function getBallAndStickParams(ctx: ThemeRegistryContext, structure: Structure) {
    return PD.clone(BallAndStickParams)
}

export type BallAndStickRepresentation = StructureRepresentation<BallAndStickParams>
export function BallAndStickRepresentation(ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, BallAndStickParams>): BallAndStickRepresentation {
    return Representation.createMulti('Ball & Stick', ctx, getParams, StructureRepresentationStateBuilder, BallAndStickVisuals as unknown as Representation.Def<Structure, BallAndStickParams>)
}

export const BallAndStickRepresentationProvider: StructureRepresentationProvider<BallAndStickParams> = {
    label: 'Ball & Stick',
    description: 'Displays atoms as spheres and bonds as cylinders.',
    factory: BallAndStickRepresentation,
    getParams: getBallAndStickParams,
    defaultValues: PD.getDefaultValues(BallAndStickParams),
    defaultColorTheme: 'element-symbol',
    defaultSizeTheme: 'physical',
    isApplicable: (structure: Structure) => structure.elementCount > 0
}