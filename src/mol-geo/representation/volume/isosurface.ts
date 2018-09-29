/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { VolumeData, VolumeIsoValue } from 'mol-model/volume'
import { Task, RuntimeContext } from 'mol-task'
import { computeMarchingCubesMesh } from '../../util/marching-cubes/algorithm';
import { Mesh } from '../../geometry/mesh/mesh';
import { VolumeVisual } from '.';
import { createMeshRenderObject, MeshRenderObject } from 'mol-gl/render-object';
import { PickingId } from '../../geometry/picking';
import { MarkerAction } from '../../geometry/marker-data';
import { Loci, EmptyLoci } from 'mol-model/loci';
import { LocationIterator } from '../../util/location-iterator';
import { NullLocation } from 'mol-model/location';
import { createIdentityTransform } from '../../geometry/transform-data';
import { createRenderableState } from '../../geometry/geometry';
import { paramDefaultValues, NumberParam } from 'mol-view/parameter';

export function computeVolumeSurface(volume: VolumeData, isoValue: VolumeIsoValue) {
    return Task.create<Mesh>('Volume Surface', async ctx => {
        ctx.update({ message: 'Marching cubes...' });

        const mesh = await computeMarchingCubesMesh({
            isoLevel: VolumeIsoValue.toAbsolute(isoValue).absoluteValue,
            scalarField: volume.data
        }).runAsChild(ctx);

        const transform = VolumeData.getGridToCartesianTransform(volume);
        ctx.update({ message: 'Transforming mesh...' });
        Mesh.transformImmediate(mesh, transform);

        return mesh;
    });
}

export const IsosurfaceParams = {
    ...Mesh.Params,
    isoValue: NumberParam('Iso Value', '', 2, -5, 5, 0.01),
}
export const DefaultIsosurfaceProps = paramDefaultValues(IsosurfaceParams)
export type IsosurfaceProps = typeof DefaultIsosurfaceProps

export default function IsosurfaceVisual(): VolumeVisual<IsosurfaceProps> {
    let renderObject: MeshRenderObject
    let currentProps = DefaultIsosurfaceProps

    return {
        get renderObject () { return renderObject },
        async createOrUpdate(ctx: RuntimeContext, props: Partial<IsosurfaceProps> = {}, volume?: VolumeData) {
            currentProps = { ...DefaultIsosurfaceProps, ...props }

            console.log('MOINMOIN')
            if (!volume) return

            const mesh = await computeVolumeSurface(volume,  VolumeIsoValue.relative(volume.dataStats, currentProps.isoValue)).runAsChild(ctx)
            if (!props.flatShaded) {
                Mesh.computeNormalsImmediate(mesh)
            }

            const locationIt = LocationIterator(1, 1, () => NullLocation)
            const transform = createIdentityTransform()

            const values = await Mesh.createValues(ctx, mesh, transform, locationIt, currentProps)
            const state = createRenderableState(currentProps)

            renderObject = createMeshRenderObject(values, state)
            console.log('renderObject', renderObject)
        },
        getLoci(pickingId: PickingId) {
            // TODO
            return EmptyLoci
        },
        mark(loci: Loci, action: MarkerAction) {
            // TODO
            return false
        },
        destroy() {
            // TODO
        }
    }
}
