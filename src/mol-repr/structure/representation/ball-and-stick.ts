/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { getElementSphereVisual, ElementSphereParams } from '../visual/element-sphere';
import { IntraUnitBondVisual, IntraUnitBondParams } from '../visual/bond-intra-unit-cylinder';
import { InterUnitBondVisual, InterUnitBondParams } from '../visual/bond-inter-unit-cylinder';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { UnitsRepresentation } from '../units-representation';
import { ComplexRepresentation } from '../complex-representation';
import { StructureRepresentation, StructureRepresentationProvider, StructureRepresentationStateBuilder } from '../representation';
import { Representation, RepresentationParamsGetter, RepresentationContext } from '../../../mol-repr/representation';
import { ThemeRegistryContext } from '../../../mol-theme/theme';
import { Structure } from '../../../mol-model/structure';
import { getUnitKindsParam } from '../params';

const BallAndStickVisuals = {
    'element-sphere': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, ElementSphereParams>) => UnitsRepresentation('Element sphere mesh', ctx, getParams, getElementSphereVisual(ctx.webgl)),
    'intra-bond': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, IntraUnitBondParams>) => UnitsRepresentation('Intra-unit bond cylinder', ctx, getParams, IntraUnitBondVisual),
    'inter-bond': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, InterUnitBondParams>) => ComplexRepresentation('Inter-unit bond cylinder', ctx, getParams, InterUnitBondVisual),
};

export const BallAndStickParams = {
    ...ElementSphereParams,
    ...IntraUnitBondParams,
    ...InterUnitBondParams,
    unitKinds: getUnitKindsParam(['atomic']),
    sizeFactor: PD.Numeric(0.15, { min: 0.01, max: 10, step: 0.01 }),
    sizeAspectRatio: PD.Numeric(2 / 3, { min: 0.01, max: 3, step: 0.01 }),
    visuals: PD.MultiSelect(['element-sphere', 'intra-bond', 'inter-bond'], PD.objectToOptions(BallAndStickVisuals))
};
export type BallAndStickParams = typeof BallAndStickParams
export function getBallAndStickParams(ctx: ThemeRegistryContext, structure: Structure) {
    return PD.clone(BallAndStickParams);
}

export type BallAndStickRepresentation = StructureRepresentation<BallAndStickParams>
export function BallAndStickRepresentation(ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, BallAndStickParams>): BallAndStickRepresentation {
    return Representation.createMulti('Ball & Stick', ctx, getParams, StructureRepresentationStateBuilder, BallAndStickVisuals as unknown as Representation.Def<Structure, BallAndStickParams>);
}

export const BallAndStickRepresentationProvider = StructureRepresentationProvider({
    name: 'ball-and-stick',
    label: 'Ball & Stick',
    description: 'Displays atoms as spheres and bonds as cylinders.',
    factory: BallAndStickRepresentation,
    getParams: getBallAndStickParams,
    defaultValues: PD.getDefaultValues(BallAndStickParams),
    defaultColorTheme: { name: 'element-symbol' },
    defaultSizeTheme: { name: 'physical' },
    isApplicable: (structure: Structure) => structure.elementCount > 0
});