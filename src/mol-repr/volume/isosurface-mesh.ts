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
import { computeMarchingCubesMesh } from 'mol-geo/util/marching-cubes/algorithm';
import { LocationIterator } from 'mol-geo/util/location-iterator';
import { VisualUpdateState } from 'mol-repr/util';
import { RepresentationContext, RepresentationParamsGetter } from 'mol-repr/representation';
import { Theme, ThemeRegistryContext } from 'mol-theme/theme';
import { VisualContext } from 'mol-repr/visual';
import { NullLocation } from 'mol-model/location';

interface VolumeIsosurfaceProps {
    isoValue: number
}

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
export function getIsosurfaceParams(ctx: ThemeRegistryContext, volume: VolumeData) {
    return PD.clone(IsosurfaceMeshParams)
}

export function IsosurfaceVisual(): VolumeVisual<IsosurfaceMeshParams> {
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

export function IsosurfaceRepresentation(ctx: RepresentationContext, getParams: RepresentationParamsGetter<VolumeData, IsosurfaceMeshParams>): VolumeRepresentation<IsosurfaceMeshParams> {
    return VolumeRepresentation('Isosurface ', ctx, getParams, IsosurfaceVisual)
}

export const IsosurfaceRepresentationProvider: VolumeRepresentationProvider<IsosurfaceMeshParams> = {
    label: 'Isosurface',
    description: 'Displays an isosurface of volumetric data.',
    factory: IsosurfaceRepresentation,
    getParams: getIsosurfaceParams,
    defaultValues: PD.getDefaultValues(IsosurfaceMeshParams),
    defaultColorTheme: 'uniform',
    defaultSizeTheme: 'uniform'
}