/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { VolumeData } from 'mol-model/volume'
import { VolumeVisual, VolumeRepresentation, VolumeRepresentationProvider } from './representation';
import { createRenderObject } from 'mol-gl/render-object';
import { EmptyLoci } from 'mol-model/loci';
import { ParamDefinition as PD } from 'mol-util/param-definition';
import { Mesh } from 'mol-geo/geometry/mesh/mesh';
import { computeMarchingCubesMesh } from 'mol-geo/util/marching-cubes/algorithm';
import { LocationIterator } from 'mol-geo/util/location-iterator';
import { createIdentityTransform } from 'mol-geo/geometry/transform-data';
import { VisualUpdateState } from 'mol-repr/util';
import { RepresentationContext, RepresentationParamsGetter } from 'mol-repr/representation';
import { Theme, ThemeRegistryContext } from 'mol-theme/theme';
import { VisualContext } from 'mol-repr/visual';

interface VolumeIsosurfaceProps {
    isoValue: number
}

export async function createVolumeIsosurface(ctx: VisualContext, volume: VolumeData, props: VolumeIsosurfaceProps, mesh?: Mesh) {
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

export const IsosurfaceParams = {
    ...Mesh.Params,
    isoValue: PD.Numeric(0.22, { min: -1, max: 1, step: 0.01 }),
}
export type IsosurfaceParams = typeof IsosurfaceParams
export function getIsosurfaceParams(ctx: ThemeRegistryContext, volume: VolumeData) {
    return PD.clone(IsosurfaceParams)
}

export function IsosurfaceVisual(): VolumeVisual<IsosurfaceParams> {
    return VolumeVisual<IsosurfaceParams>({
        defaultProps: PD.getDefaultValues(IsosurfaceParams),
        createGeometry: createVolumeIsosurface,
        getLoci: () => EmptyLoci,
        mark: () => false,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<IsosurfaceParams>, currentProps: PD.Values<IsosurfaceParams>) => {
            if (newProps.isoValue !== currentProps.isoValue) state.createGeometry = true
        },
        createRenderObject: (geometry: Mesh, locationIt: LocationIterator, theme: Theme, props: PD.Values<IsosurfaceParams>) => {
            const transform = createIdentityTransform()
            const values = Mesh.Utils.createValues(geometry, transform, locationIt, theme, props)
            const state = Mesh.Utils.createRenderableState(props)
            return createRenderObject('mesh', values, state)
        },
        updateValues: Mesh.Utils.updateValues,
        updateBoundingSphere: Mesh.Utils.updateBoundingSphere,
        updateRenderableState: Mesh.Utils.updateRenderableState
    })
}

export function IsosurfaceRepresentation(ctx: RepresentationContext, getParams: RepresentationParamsGetter<VolumeData, IsosurfaceParams>): VolumeRepresentation<IsosurfaceParams> {
    return VolumeRepresentation('Isosurface', ctx, getParams, IsosurfaceVisual)
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