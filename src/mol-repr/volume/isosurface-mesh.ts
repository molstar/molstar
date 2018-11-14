/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { VolumeData } from 'mol-model/volume'
import { VolumeVisual, VolumeRepresentation } from './representation';
import { createMeshRenderObject } from 'mol-gl/render-object';
import { EmptyLoci } from 'mol-model/loci';
import { ParamDefinition as PD } from 'mol-util/param-definition';
import { Mesh } from 'mol-geo/geometry/mesh/mesh';
import { computeMarchingCubesMesh } from 'mol-geo/util/marching-cubes/algorithm';
import { LocationIterator } from 'mol-geo/util/location-iterator';
import { createIdentityTransform } from 'mol-geo/geometry/transform-data';
import { createRenderableState } from 'mol-geo/geometry/geometry';
import { VisualUpdateState } from 'mol-repr/util';
import { VisualContext } from 'mol-repr/representation';
import { Theme, ThemeRegistryContext } from 'mol-theme/theme';

interface VolumeIsosurfaceProps {
    isoValueAbsolute: number
}

export async function createVolumeIsosurface(ctx: VisualContext, volume: VolumeData, props: VolumeIsosurfaceProps, mesh?: Mesh) {
    ctx.runtime.update({ message: 'Marching cubes...' });

    const surface = await computeMarchingCubesMesh({
        isoLevel: props.isoValueAbsolute,
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
    isoValueAbsolute: PD.Range('Iso Value Absolute', '', 0.22, -1, 1, 0.01),
    isoValueRelative: PD.Range('Iso Value Relative', '', 2, -10, 10, 0.1),
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
        setUpdateState: (state: VisualUpdateState, newProps: PD.DefaultValues<IsosurfaceParams>, currentProps: PD.DefaultValues<IsosurfaceParams>) => {
            if (newProps.isoValueAbsolute !== currentProps.isoValueAbsolute) state.createGeometry = true
        },
        createRenderObject: async (ctx: VisualContext, geometry: Mesh, locationIt: LocationIterator, theme: Theme, props: PD.DefaultValues<IsosurfaceParams>) => {
            const transform = createIdentityTransform()
            const values = await Mesh.createValues(ctx.runtime, geometry, transform, locationIt, theme, props)
            const state = createRenderableState(props)
            return createMeshRenderObject(values, state)
        },
        updateValues: Mesh.updateValues
    })
}

export function IsosurfaceRepresentation(): VolumeRepresentation<IsosurfaceParams> {
    return VolumeRepresentation('Isosurface', getIsosurfaceParams, IsosurfaceVisual)
}