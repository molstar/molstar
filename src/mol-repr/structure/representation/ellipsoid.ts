/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { RepresentationParamsGetter, RepresentationContext, Representation } from '../../../mol-repr/representation';
import { ThemeRegistryContext } from '../../../mol-theme/theme';
import { Structure } from '../../../mol-model/structure';
import { UnitsRepresentation, StructureRepresentation, StructureRepresentationStateBuilder, StructureRepresentationProvider } from '../../../mol-repr/structure/representation';
import { EllipsoidMeshParams, EllipsoidMeshVisual } from '../visual/ellipsoid-mesh';
import { UnitKind, UnitKindOptions } from '../../../mol-repr/structure/visual/util/common';
import { AtomSiteAnisotrop } from '../../../mol-model-formats/structure/mmcif/anisotropic';

const EllipsoidVisuals = {
    'ellipsoid-mesh': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, EllipsoidParams>) => UnitsRepresentation('Ellipsoid Mesh', ctx, getParams, EllipsoidMeshVisual),
}

export const EllipsoidParams = {
    ...EllipsoidMeshParams,
    unitKinds: PD.MultiSelect<UnitKind>(['atomic'], UnitKindOptions),
}
export type EllipsoidParams = typeof EllipsoidMeshParams
export function getEllipsoidParams(ctx: ThemeRegistryContext, structure: Structure) {
    return PD.clone(EllipsoidMeshParams)
}

export type EllipsoidRepresentation = StructureRepresentation<EllipsoidParams>
export function EllipsoidRepresentation(ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, EllipsoidParams>): EllipsoidRepresentation {
    return Representation.createMulti('Ellipsoid', ctx, getParams, StructureRepresentationStateBuilder, EllipsoidVisuals as unknown as Representation.Def<Structure, EllipsoidParams>)
}

export const EllipsoidRepresentationProvider: StructureRepresentationProvider<EllipsoidParams> = {
    label: 'Ellipsoid',
    description: 'Displays anisotropic displacement ellipsoids of atomic elements.',
    factory: EllipsoidRepresentation,
    getParams: getEllipsoidParams,
    defaultValues: PD.getDefaultValues(EllipsoidMeshParams),
    defaultColorTheme: 'element-symbol',
    defaultSizeTheme: 'uniform',
    isApplicable: (structure: Structure) => structure.elementCount > 0 && structure.models.some(m => m.customProperties.has(AtomSiteAnisotrop.Descriptor))
}