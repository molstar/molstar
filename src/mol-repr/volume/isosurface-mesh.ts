/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { VolumeData } from 'mol-model/volume'
import { VolumeVisual, VolumeRepresentation, VolumeRepresentationProvider } from './representation';
import { EmptyLoci } from 'mol-model/loci';
import { ParamDefinition as PD } from 'mol-util/param-definition';
import { Mesh } from 'mol-geo/geometry/mesh/mesh';
import { computeMarchingCubesMesh, computeMarchingCubesLines } from 'mol-geo/util/marching-cubes/algorithm';
import { LocationIterator } from 'mol-geo/util/location-iterator';
import { VisualUpdateState } from 'mol-repr/util';
import { RepresentationContext, RepresentationParamsGetter, Representation } from 'mol-repr/representation';
import { Theme, ThemeRegistryContext } from 'mol-theme/theme';
import { VisualContext } from 'mol-repr/visual';
import { NullLocation } from 'mol-model/location';
import { Lines } from 'mol-geo/geometry/lines/lines';

interface VolumeIsosurfaceProps {
    isoValue: number
}

//

export async function createVolumeIsosurfaceMesh(ctx: VisualContext, volume: VolumeData, theme: Theme, props: VolumeIsosurfaceProps, mesh?: Mesh) {
    ctx.runtime.update({ message: 'Marching cubes...' });

    const surface = await computeMarchingCubesMesh({
        isoLevel: props.isoValue,
        scalarField: volume.data
    }, mesh).runAsChild(ctx.runtime);

    const transform = VolumeData.getGridToCartesianTransform(volume);
    ctx.runtime.update({ message: 'Transforming mesh...' });
    Mesh.transformImmediate(surface, transform);
    Mesh.computeNormalsImmediate(surface)

    return surface;
}

export const IsosurfaceMeshParams = {
    ...Mesh.Params,
    isoValue: PD.Numeric(0.22, { min: -1, max: 1, step: 0.01 }),
}
export type IsosurfaceMeshParams = typeof IsosurfaceMeshParams

export function IsosurfaceMeshVisual(): VolumeVisual<IsosurfaceMeshParams> {
    return VolumeVisual<Mesh, IsosurfaceMeshParams>({
        defaultProps: PD.getDefaultValues(IsosurfaceMeshParams),
        createGeometry: createVolumeIsosurfaceMesh,
        createLocationIterator: (volume: VolumeData) => LocationIterator(1, 1, () => NullLocation),
        getLoci: () => EmptyLoci,
        mark: () => false,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<IsosurfaceMeshParams>, currentProps: PD.Values<IsosurfaceMeshParams>) => {
            if (newProps.isoValue !== currentProps.isoValue) state.createGeometry = true
        },
        geometryUtils: Mesh.Utils
    })
}

//

export async function createVolumeIsosurfaceWireframe(ctx: VisualContext, volume: VolumeData, theme: Theme, props: VolumeIsosurfaceProps, lines?: Lines) {
    ctx.runtime.update({ message: 'Marching cubes...' });

    const params = {
        isoLevel: props.isoValue,
        scalarField: volume.data
    }
    const wireframe = await computeMarchingCubesLines(params, lines).runAsChild(ctx.runtime)

    const transform = VolumeData.getGridToCartesianTransform(volume);
    Lines.transformImmediate(wireframe, transform)

    return wireframe;
}

export const IsosurfaceWireframeParams = {
    ...Lines.Params,
    isoValue: PD.Numeric(0.22, { min: -1, max: 1, step: 0.01 }),
}
export type IsosurfaceWireframeParams = typeof IsosurfaceWireframeParams

export function IsosurfaceWireframeVisual(): VolumeVisual<IsosurfaceWireframeParams> {
    return VolumeVisual<Lines, IsosurfaceWireframeParams>({
        defaultProps: PD.getDefaultValues(IsosurfaceWireframeParams),
        createGeometry: createVolumeIsosurfaceWireframe,
        createLocationIterator: (volume: VolumeData) => LocationIterator(1, 1, () => NullLocation),
        getLoci: () => EmptyLoci,
        mark: () => false,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<IsosurfaceWireframeParams>, currentProps: PD.Values<IsosurfaceWireframeParams>) => {
            if (newProps.isoValue !== currentProps.isoValue) state.createGeometry = true
        },
        geometryUtils: Lines.Utils
    })
}

//

const IsosurfaceVisuals = {
    'solid': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<VolumeData, IsosurfaceMeshParams>) => VolumeRepresentation('Isosurface mesh', ctx, getParams, IsosurfaceMeshVisual),
    'wireframe': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<VolumeData, IsosurfaceWireframeParams>) => VolumeRepresentation('Isosurface wireframe', ctx, getParams, IsosurfaceWireframeVisual),
}
type IsosurfaceVisualName = keyof typeof IsosurfaceVisuals
const IsosurfaceVisualOptions = Object.keys(IsosurfaceVisuals).map(name => [name, name] as [IsosurfaceVisualName, string])

export const IsosurfaceParams = {
    ...IsosurfaceMeshParams,
    ...IsosurfaceWireframeParams,
    visuals: PD.MultiSelect<IsosurfaceVisualName>(['solid'], IsosurfaceVisualOptions),
}
export type IsosurfaceParams = typeof IsosurfaceParams
export function getIsosurfaceParams(ctx: ThemeRegistryContext, volume: VolumeData) {
    return PD.clone(IsosurfaceParams)
}

export type IsosurfaceRepresentation = VolumeRepresentation<IsosurfaceParams>
export function IsosurfaceRepresentation(ctx: RepresentationContext, getParams: RepresentationParamsGetter<VolumeData, IsosurfaceParams>): IsosurfaceRepresentation {
    return Representation.createMulti('Isosurface', ctx, getParams, Representation.StateBuilder, IsosurfaceVisuals as unknown as Representation.Def<VolumeData, IsosurfaceParams>)
}

export const IsosurfaceRepresentationProvider: VolumeRepresentationProvider<IsosurfaceParams> = {
    label: 'Isosurface',
    description: 'Displays an isosurface of volumetric data.',
    factory: IsosurfaceRepresentation,
    getParams: getIsosurfaceParams,
    defaultValues: PD.getDefaultValues(IsosurfaceParams),
    defaultColorTheme: 'uniform',
    defaultSizeTheme: 'uniform'
}