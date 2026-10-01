/**
 * Copyright (c) 2019-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { type RepresentationParamsGetter, type RepresentationContext, Representation } from '../../representation.js';
import type { ThemeRegistryContext } from '@molstar/graphics/theme/theme';
import { Structure } from '@molstar/model/model/structure';
import { UnitsRepresentation, type StructureRepresentation, StructureRepresentationStateBuilder, StructureRepresentationProvider, ComplexRepresentation } from '../representation.js';
import { EllipsoidMeshParams, EllipsoidMeshVisual, type StructureEllipsoidMeshParams, StructureEllipsoidMeshVisual } from '../visual/ellipsoid-mesh.js';
import { AtomSiteAnisotrop } from '@molstar/model/formats/structure/property/anisotropic';
import { IntraUnitBondCylinderParams, IntraUnitBondCylinderVisual, type StructureIntraUnitBondCylinderParams, StructureIntraUnitBondCylinderVisual } from '../visual/bond-intra-unit-cylinder.js';
import { InterUnitBondCylinderVisual, InterUnitBondCylinderParams } from '../visual/bond-inter-unit-cylinder.js';
import { getUnitKindsParam } from '../params.js';
import { BaseGeometry } from '@molstar/graphics/geo/geometry/base';

const EllipsoidVisuals = {
    'ellipsoid-mesh': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, EllipsoidMeshParams>) => UnitsRepresentation('Ellipsoid Mesh', ctx, getParams, EllipsoidMeshVisual),
    'intra-bond': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, IntraUnitBondCylinderParams>) => UnitsRepresentation('Intra-unit bond cylinder', ctx, getParams, IntraUnitBondCylinderVisual),
    'inter-bond': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, InterUnitBondCylinderParams>) => ComplexRepresentation('Inter-unit bond cylinder', ctx, getParams, InterUnitBondCylinderVisual),
    'structure-ellipsoid-mesh': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, StructureEllipsoidMeshParams>) => ComplexRepresentation('Structure Ellipsoid Mesh', ctx, getParams, StructureEllipsoidMeshVisual),
    'structure-intra-bond': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, StructureIntraUnitBondCylinderParams>) => ComplexRepresentation('Structure intra-unit bond cylinder', ctx, getParams, StructureIntraUnitBondCylinderVisual),
};

export const EllipsoidParams = {
    ...EllipsoidMeshParams,
    ...IntraUnitBondCylinderParams,
    ...InterUnitBondCylinderParams,
    includeParent: PD.Boolean(false),
    adjustCylinderLength: PD.Boolean(false, { isHidden: true }), // not useful here
    unitKinds: getUnitKindsParam(['atomic']),
    sizeFactor: PD.Numeric(1, { min: 0.01, max: 10, step: 0.01 }),
    sizeAspectRatio: PD.Numeric(0.1, { min: 0.01, max: 3, step: 0.01 }),
    linkCap: PD.Boolean(true),
    visuals: PD.MultiSelect(['ellipsoid-mesh', 'intra-bond', 'inter-bond'], PD.objectToOptions(EllipsoidVisuals)),
    bumpFrequency: PD.Numeric(0, { min: 0, max: 10, step: 0.1 }, BaseGeometry.ShadingCategory),
};
export type EllipsoidParams = typeof EllipsoidParams
export function getEllipsoidParams(ctx: ThemeRegistryContext, structure: Structure) {
    let params = EllipsoidParams;
    const size = Structure.getSize(structure);
    if (size >= Structure.Size.Huge) {
        params = PD.clone(params);
        params.visuals.defaultValue = ['ellipsoid-mesh', 'intra-bond'];
    }
    return params;
}

export type EllipsoidRepresentation = StructureRepresentation<EllipsoidParams>
export function EllipsoidRepresentation(ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, EllipsoidParams>): EllipsoidRepresentation {
    return Representation.createMulti('Ellipsoid', ctx, getParams, StructureRepresentationStateBuilder, EllipsoidVisuals as unknown as Representation.Def<Structure, EllipsoidParams>);
}

export const EllipsoidRepresentationProvider = StructureRepresentationProvider({
    name: 'ellipsoid',
    label: 'Ellipsoid',
    description: 'Displays anisotropic displacement ellipsoids of atomic elements plus bonds as cylinders.',
    factory: EllipsoidRepresentation,
    getParams: getEllipsoidParams,
    defaultValues: PD.getDefaultValues(EllipsoidParams),
    defaultColorTheme: { name: 'element-symbol' },
    defaultSizeTheme: { name: 'uniform' },
    isApplicable: (structure: Structure) => structure.elementCount > 0 && structure.models.some(m => AtomSiteAnisotrop.Provider.isApplicable(m)),
    getData: (structure: Structure, props: PD.Values<EllipsoidParams>) => {
        return props.includeParent ? structure.asParent() : structure;
    },
    mustRecreate: (oldProps: PD.Values<EllipsoidParams>, newProps: PD.Values<EllipsoidParams>) => {
        return oldProps.includeParent !== newProps.includeParent;
    }
});