/**
 * Copyright (c) 2019-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { UnitsRepresentation } from '../units-representation.js';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { type StructureRepresentation, StructureRepresentationProvider, StructureRepresentationStateBuilder } from '../representation.js';
import { Representation, type RepresentationParamsGetter, type RepresentationContext } from '../../representation.js';
import type { ThemeRegistryContext } from '@molstar/graphics/theme/theme';
import type { Structure } from '@molstar/model/model/structure';
import { OrientationEllipsoidMeshParams, OrientationEllipsoidMeshVisual } from '../visual/orientation-ellipsoid-mesh.js';
import { BaseGeometry } from '@molstar/graphics/geo/geometry/base';

const OrientationVisuals = {
    'orientation-ellipsoid-mesh': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, OrientationEllipsoidMeshParams>) => UnitsRepresentation('Orientation ellipsoid mesh', ctx, getParams, OrientationEllipsoidMeshVisual),
};

export const OrientationParams = {
    ...OrientationEllipsoidMeshParams,
    visuals: PD.MultiSelect(['orientation-ellipsoid-mesh'], PD.objectToOptions(OrientationVisuals)),
    bumpFrequency: PD.Numeric(1, { min: 0, max: 10, step: 0.1 }, BaseGeometry.ShadingCategory),
};
export type OrientationParams = typeof OrientationParams
export function getOrientationParams(ctx: ThemeRegistryContext, structure: Structure) {
    return OrientationParams;
}

export type OrientationRepresentation = StructureRepresentation<OrientationParams>
export function OrientationRepresentation(ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, OrientationParams>): OrientationRepresentation {
    return Representation.createMulti('Orientation', ctx, getParams, StructureRepresentationStateBuilder, OrientationVisuals as unknown as Representation.Def<Structure, OrientationParams>);
}

export const OrientationRepresentationProvider = StructureRepresentationProvider({
    name: 'orientation',
    label: 'Orientation',
    description: 'Displays orientation ellipsoids for polymer chains.',
    factory: OrientationRepresentation,
    getParams: getOrientationParams,
    defaultValues: PD.getDefaultValues(OrientationParams),
    defaultColorTheme: { name: 'chain-id' },
    defaultSizeTheme: { name: 'uniform' },
    isApplicable: (structure: Structure) => structure.elementCount > 0
});