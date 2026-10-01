/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Gianluca Tomasello <giagitom@gmail.com>
 */

import { GaussianSurfaceMeshParams, type StructureGaussianSurfaceMeshParams, StructureGaussianSurfaceVisual, GaussianSurfaceVisual } from '../visual/gaussian-surface-mesh.js';
import { UnitsRepresentation } from '../units-representation.js';
import { GaussianWireframeVisual, GaussianWireframeParams, StructureGaussianWireframeVisual } from '../visual/gaussian-surface-wireframe.js';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { type StructureRepresentation, StructureRepresentationProvider, StructureRepresentationStateBuilder, ComplexRepresentation } from '../representation.js';
import { Representation, type RepresentationParamsGetter, type RepresentationContext } from '../../representation.js';
import type { ThemeRegistryContext } from '@molstar/graphics/theme/theme';
import type { Structure } from '@molstar/model/model/structure';
import { BaseGeometry } from '@molstar/graphics/geo/geometry/base';

const GaussianSurfaceVisuals = {
    'gaussian-surface-mesh': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, GaussianSurfaceMeshParams>) => UnitsRepresentation('Gaussian surface mesh', ctx, getParams, GaussianSurfaceVisual),
    'structure-gaussian-surface-mesh': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, StructureGaussianSurfaceMeshParams>) => ComplexRepresentation('Structure-Gaussian surface mesh', ctx, getParams, StructureGaussianSurfaceVisual),
    'gaussian-surface-wireframe': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, GaussianWireframeParams>) => UnitsRepresentation('Gaussian surface wireframe', ctx, getParams, GaussianWireframeVisual),
    'structure-gaussian-surface-wireframe': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, GaussianWireframeParams>) => ComplexRepresentation('Structure-Gaussian surface wireframe', ctx, getParams, StructureGaussianWireframeVisual),
};

export const GaussianSurfaceParams = {
    ...GaussianSurfaceMeshParams,
    ...GaussianWireframeParams,
    visuals: PD.MultiSelect(['gaussian-surface-mesh'], PD.objectToOptions(GaussianSurfaceVisuals)),
    bumpFrequency: PD.Numeric(1, { min: 0, max: 10, step: 0.1 }, BaseGeometry.ShadingCategory),
    density: PD.Numeric(0.5, { min: 0, max: 1, step: 0.01 }, BaseGeometry.ShadingCategory),
};
export type GaussianSurfaceParams = typeof GaussianSurfaceParams
export function getGaussianSurfaceParams(ctx: ThemeRegistryContext, structure: Structure) {
    return GaussianSurfaceParams;
}

export type GaussianSurfaceRepresentation = StructureRepresentation<GaussianSurfaceParams>
export function GaussianSurfaceRepresentation(ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, GaussianSurfaceParams>): GaussianSurfaceRepresentation {
    return Representation.createMulti('Gaussian Surface', ctx, getParams, StructureRepresentationStateBuilder, GaussianSurfaceVisuals as unknown as Representation.Def<Structure, GaussianSurfaceParams>);
}

export const GaussianSurfaceRepresentationProvider = StructureRepresentationProvider({
    name: 'gaussian-surface',
    label: 'Gaussian Surface',
    description: 'Displays a gaussian molecular surface.',
    factory: GaussianSurfaceRepresentation,
    getParams: getGaussianSurfaceParams,
    defaultValues: PD.getDefaultValues(GaussianSurfaceParams),
    defaultColorTheme: { name: 'chain-id' },
    defaultSizeTheme: { name: 'physical' },
    isApplicable: (structure: Structure) => structure.elementCount > 0
});