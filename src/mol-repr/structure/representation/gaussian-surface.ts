/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { GaussianSurfaceMeshVisual, GaussianSurfaceTextureMeshVisual, GaussianSurfaceMeshParams } from '../visual/gaussian-surface-mesh';
import { UnitsRepresentation } from '../units-representation';
import { GaussianWireframeVisual, GaussianWireframeParams } from '../visual/gaussian-surface-wireframe';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { StructureRepresentation, StructureRepresentationProvider, StructureRepresentationStateBuilder } from '../representation';
import { Representation, RepresentationParamsGetter, RepresentationContext } from '../../../mol-repr/representation';
import { ThemeRegistryContext } from '../../../mol-theme/theme';
import { Structure } from '../../../mol-model/structure';

const GaussianSurfaceVisuals = {
    'gaussian-surface-mesh': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, GaussianSurfaceMeshParams>) => UnitsRepresentation('Gaussian surface', ctx, getParams, GaussianSurfaceMeshVisual),
    'gaussian-surface-texture-mesh': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, GaussianSurfaceMeshParams>) => UnitsRepresentation('Gaussian surface', ctx, getParams, GaussianSurfaceTextureMeshVisual),
    'gaussian-wireframe': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, GaussianWireframeParams>) => UnitsRepresentation('Gaussian wireframe', ctx, getParams, GaussianWireframeVisual),
}
type GaussianSurfaceVisualName = keyof typeof GaussianSurfaceVisuals
const GaussianSurfaceVisualOptions = Object.keys(GaussianSurfaceVisuals).map(name => [name, name] as [GaussianSurfaceVisualName, string])

export const GaussianSurfaceParams = {
    ...GaussianSurfaceMeshParams,
    ...GaussianWireframeParams,
    visuals: PD.MultiSelect<GaussianSurfaceVisualName>(['gaussian-surface-mesh'], GaussianSurfaceVisualOptions),
}
export type GaussianSurfaceParams = typeof GaussianSurfaceParams
export function getGaussianSurfaceParams(ctx: ThemeRegistryContext, structure: Structure) {
    return PD.clone(GaussianSurfaceParams)
}

export type GaussianSurfaceRepresentation = StructureRepresentation<GaussianSurfaceParams>
export function GaussianSurfaceRepresentation(ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, GaussianSurfaceParams>): GaussianSurfaceRepresentation {
    return Representation.createMulti('Gaussian Surface', ctx, getParams, StructureRepresentationStateBuilder, GaussianSurfaceVisuals as unknown as Representation.Def<Structure, GaussianSurfaceParams>)
}

export const GaussianSurfaceRepresentationProvider: StructureRepresentationProvider<GaussianSurfaceParams> = {
    label: 'Gaussian Surface',
    description: 'Displays a gaussian molecular surface.',
    factory: GaussianSurfaceRepresentation,
    getParams: getGaussianSurfaceParams,
    defaultValues: PD.getDefaultValues(GaussianSurfaceParams),
    defaultColorTheme: 'polymer-id',
    defaultSizeTheme: 'uniform',
    isApplicable: (structure: Structure) => structure.elementCount > 0
}