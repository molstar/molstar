/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { IntraUnitBondCylinderVisual, IntraUnitBondCylinderParams, type StructureIntraUnitBondCylinderParams, StructureIntraUnitBondCylinderVisual } from '../visual/bond-intra-unit-cylinder.js';
import { InterUnitBondCylinderParams, InterUnitBondCylinderVisual } from '../visual/bond-inter-unit-cylinder.js';
import { ElementSphereVisual, ElementSphereParams, StructureElementSphereVisual } from '../visual/element-sphere.js';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { UnitsRepresentation } from '../units-representation.js';
import { ComplexRepresentation } from '../complex-representation.js';
import { type StructureRepresentation, StructureRepresentationProvider, StructureRepresentationStateBuilder } from '../representation.js';
import { Representation, type RepresentationParamsGetter, type RepresentationContext } from '../../representation.js';
import type { ThemeRegistryContext } from '@molstar/graphics/theme/theme';
import { Structure } from '@molstar/model/model/structure';
import { getUnitKindsParam } from '../params.js';
import { BaseGeometry } from '@molstar/graphics/geo/geometry/base';

const BallAndStickVisuals = {
    'element-sphere': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, ElementSphereParams>) => UnitsRepresentation('Element sphere', ctx, getParams, ElementSphereVisual),
    'intra-bond': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, IntraUnitBondCylinderParams>) => UnitsRepresentation('Intra-unit bond cylinder', ctx, getParams, IntraUnitBondCylinderVisual),
    'inter-bond': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, InterUnitBondCylinderParams>) => ComplexRepresentation('Inter-unit bond cylinder', ctx, getParams, InterUnitBondCylinderVisual),
    'structure-element-sphere': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, ElementSphereParams>) => ComplexRepresentation('Structure element sphere', ctx, getParams, StructureElementSphereVisual),
    'structure-intra-bond': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, StructureIntraUnitBondCylinderParams>) => ComplexRepresentation('Structure intra-unit bond cylinder', ctx, getParams, StructureIntraUnitBondCylinderVisual),
};

export const BallAndStickParams = {
    ...ElementSphereParams,
    traceOnly: PD.Boolean(false, { isHidden: true }), // not useful here
    ...IntraUnitBondCylinderParams,
    ...InterUnitBondCylinderParams,
    includeParent: PD.Boolean(false),
    unitKinds: getUnitKindsParam(['atomic']),
    sizeFactor: PD.Numeric(0.15, { min: 0.01, max: 10, step: 0.01 }),
    sizeAspectRatio: PD.Numeric(2 / 3, { min: 0.01, max: 3, step: 0.01 }),
    visuals: PD.MultiSelect(['element-sphere', 'intra-bond', 'inter-bond'], PD.objectToOptions(BallAndStickVisuals)),
    bumpFrequency: PD.Numeric(5, { min: 0, max: 10, step: 0.1 }, BaseGeometry.ShadingCategory),
    density: PD.Numeric(0.1, { min: 0, max: 1, step: 0.01 }, BaseGeometry.ShadingCategory),
};
export type BallAndStickParams = typeof BallAndStickParams
export function getBallAndStickParams(ctx: ThemeRegistryContext, structure: Structure) {
    let params = BallAndStickParams;
    const size = Structure.getSize(structure);
    if (size >= Structure.Size.Huge) {
        params = PD.clone(params);
        params.visuals.defaultValue = ['element-sphere', 'intra-bond'];
    }
    return params;
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
    isApplicable: (structure: Structure) => structure.elementCount > 0,
    getData: (structure: Structure, props: PD.Values<BallAndStickParams>) => {
        return props.includeParent ? structure.asParent() : structure;
    },
    mustRecreate: (oldProps: PD.Values<BallAndStickParams>, newProps: PD.Values<BallAndStickParams>) => {
        return oldProps.includeParent !== newProps.includeParent;
    }
});