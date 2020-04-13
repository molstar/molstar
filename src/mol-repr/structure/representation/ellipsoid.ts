/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { RepresentationParamsGetter, RepresentationContext, Representation } from '../../../mol-repr/representation';
import { ThemeRegistryContext } from '../../../mol-theme/theme';
import { Structure } from '../../../mol-model/structure';
import { UnitsRepresentation, StructureRepresentation, StructureRepresentationStateBuilder, StructureRepresentationProvider, ComplexRepresentation } from '../../../mol-repr/structure/representation';
import { EllipsoidMeshParams, EllipsoidMeshVisual } from '../visual/ellipsoid-mesh';
import { AtomSiteAnisotrop } from '../../../mol-model-formats/structure/property/anisotropic';
import { IntraUnitBondParams, IntraUnitBondVisual } from '../visual/bond-intra-unit-cylinder';
import { InterUnitBondParams, InterUnitBondVisual } from '../visual/bond-inter-unit-cylinder';
import { getUnitKindsParam } from '../params';

const EllipsoidVisuals = {
    'ellipsoid-mesh': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, EllipsoidMeshParams>) => UnitsRepresentation('Ellipsoid Mesh', ctx, getParams, EllipsoidMeshVisual),
    'intra-bond': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, IntraUnitBondParams>) => UnitsRepresentation('Intra-unit bond cylinder', ctx, getParams, IntraUnitBondVisual),
    'inter-bond': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, InterUnitBondParams>) => ComplexRepresentation('Inter-unit bond cylinder', ctx, getParams, InterUnitBondVisual),
};

export const EllipsoidParams = {
    ...EllipsoidMeshParams,
    ...IntraUnitBondParams,
    ...InterUnitBondParams,
    unitKinds: getUnitKindsParam(['atomic']),
    sizeFactor: PD.Numeric(1, { min: 0.01, max: 10, step: 0.01 }),
    sizeAspectRatio: PD.Numeric(0.1, { min: 0.01, max: 3, step: 0.01 }),
    linkCap: PD.Boolean(true),
    visuals: PD.MultiSelect(['ellipsoid-mesh', 'intra-bond', 'inter-bond'], PD.objectToOptions(EllipsoidVisuals)),
};
export type EllipsoidParams = typeof EllipsoidParams
export function getEllipsoidParams(ctx: ThemeRegistryContext, structure: Structure) {
    return PD.clone(EllipsoidParams);
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
    isApplicable: (structure: Structure) => structure.elementCount > 0 && structure.models.some(m => AtomSiteAnisotrop.Provider.isApplicable(m))
});