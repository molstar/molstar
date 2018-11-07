/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { VolumeData } from 'mol-model/volume'
import { VolumeVisual, VolumeRepresentation } from './index';
import { createMeshRenderObject } from 'mol-gl/render-object';
import { Loci, EmptyLoci } from 'mol-model/loci';
import { ParamDefinition as PD } from 'mol-util/param-definition';
import { Mesh } from 'mol-geo/geometry/mesh/mesh';
import { computeMarchingCubesMesh } from 'mol-geo/util/marching-cubes/algorithm';
import { LocationIterator } from 'mol-geo/util/location-iterator';
import { createIdentityTransform } from 'mol-geo/geometry/transform-data';
import { createRenderableState, Theme } from 'mol-geo/geometry/geometry';
import { PickingId } from 'mol-geo/geometry/picking';
import { MarkerAction } from 'mol-geo/geometry/marker-data';
import { VisualUpdateState } from 'mol-repr/util';
import { RepresentationContext, VisualContext } from 'mol-repr';

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
export const DefaultIsosurfaceProps = PD.getDefaultValues(IsosurfaceParams)
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
        createRenderObject: async (ctx: VisualContext, geometry: Mesh, locationIt: LocationIterator, theme: Theme, props: IsosurfaceProps) => {
            const transform = createIdentityTransform()
            const values = await Mesh.createValues(ctx.runtime, geometry, transform, locationIt, theme, props)
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
        createOrUpdate: (ctx: RepresentationContext, props: Partial<IsosurfaceProps> = {}, volume?: VolumeData) => {
            currentProps = Object.assign({}, DefaultIsosurfaceProps, currentProps, props)
            return volumeRepr.createOrUpdate(ctx, currentProps, volume)
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