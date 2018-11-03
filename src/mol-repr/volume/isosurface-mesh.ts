/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { VolumeData } from 'mol-model/volume'
import { RuntimeContext } from 'mol-task'
import { VolumeVisual, VolumeRepresentation } from './index';
import { createMeshRenderObject } from 'mol-gl/render-object';
import { Loci, EmptyLoci } from 'mol-model/loci';
import { paramDefaultValues, RangeParam } from 'mol-util/parameter';
import { Mesh } from 'mol-geo/geometry/mesh/mesh';
import { computeMarchingCubesMesh } from 'mol-geo/util/marching-cubes/algorithm';
import { LocationIterator } from 'mol-geo/util/location-iterator';
import { createIdentityTransform } from 'mol-geo/geometry/transform-data';
import { createRenderableState } from 'mol-geo/geometry/geometry';
import { PickingId } from 'mol-geo/geometry/picking';
import { MarkerAction } from 'mol-geo/geometry/marker-data';
import { VisualUpdateState } from 'mol-repr/util';

interface VolumeIsosurfaceProps {
    isoValueAbsolute: number
}

export async function createVolumeIsosurface(ctx: RuntimeContext, volume: VolumeData, props: VolumeIsosurfaceProps, mesh?: Mesh) {
    ctx.update({ message: 'Marching cubes...' });

    const surface = await computeMarchingCubesMesh({
        isoLevel: props.isoValueAbsolute,
        scalarField: volume.data
    }, mesh).runAsChild(ctx);

    const transform = VolumeData.getGridToCartesianTransform(volume);
    ctx.update({ message: 'Transforming mesh...' });
    Mesh.transformImmediate(surface, transform);
    Mesh.computeNormalsImmediate(surface)

    return surface;
}

export const IsosurfaceParams = {
    ...Mesh.Params,
    isoValueAbsolute: RangeParam('Iso Value Absolute', '', 0.22, -1, 1, 0.01),
    isoValueRelative: RangeParam('Iso Value Relative', '', 2, -10, 10, 0.1),
}
export const DefaultIsosurfaceProps = paramDefaultValues(IsosurfaceParams)
export type IsosurfaceProps = typeof DefaultIsosurfaceProps

export function IsosurfaceVisual(): VolumeVisual<IsosurfaceProps> {
    return VolumeVisual<IsosurfaceProps>({
        defaultProps: DefaultIsosurfaceProps,
        createGeometry: createVolumeIsosurface,
        getLoci: () => EmptyLoci,
        mark: () => false,
        setUpdateState: (state: VisualUpdateState, newProps: IsosurfaceProps, currentProps: IsosurfaceProps) => {
            if (newProps.isoValueAbsolute !== currentProps.isoValueAbsolute) state.createGeometry = true
        },
        createRenderObject: async (ctx: RuntimeContext, geometry: Mesh, locationIt: LocationIterator, props: IsosurfaceProps) => {
            const transform = createIdentityTransform()
            const values = await Mesh.createValues(ctx, geometry, transform, locationIt, props)
            const state = createRenderableState(props)
            return createMeshRenderObject(values, state)
        },
        updateValues: Mesh.updateValues
    })
}

export function IsosurfaceRepresentation(): VolumeRepresentation<IsosurfaceProps> {
    let currentProps: IsosurfaceProps
    const volumeRepr = VolumeRepresentation(IsosurfaceVisual)
    return {
        label: 'Isosurface',
        params: IsosurfaceParams,
        get renderObjects() {
            return [ ...volumeRepr.renderObjects ]
        },
        get props() {
            return { ...volumeRepr.props }
        },
        createOrUpdate: (props: Partial<IsosurfaceProps> = {}, volume?: VolumeData) => {
            currentProps = Object.assign({}, DefaultIsosurfaceProps, currentProps, props)
            return volumeRepr.createOrUpdate(currentProps, volume)
        },
        getLoci: (pickingId: PickingId) => {
            return volumeRepr.getLoci(pickingId)
        },
        mark: (loci: Loci, action: MarkerAction) => {
            return volumeRepr.mark(loci, action)
        },
        destroy() {
            volumeRepr.destroy()
        }
    }
}