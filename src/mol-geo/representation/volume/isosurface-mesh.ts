/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { VolumeData, VolumeIsoValue } from 'mol-model/volume'
import { RuntimeContext } from 'mol-task'
import { computeMarchingCubesMesh } from '../../util/marching-cubes/algorithm';
import { Mesh } from '../../geometry/mesh/mesh';
import { VolumeVisual, VolumeRepresentation } from '.';
import { createMeshRenderObject, MeshRenderObject } from 'mol-gl/render-object';
import { PickingId } from '../../geometry/picking';
import { MarkerAction } from '../../geometry/marker-data';
import { Loci, EmptyLoci } from 'mol-model/loci';
import { LocationIterator } from '../../util/location-iterator';
import { NullLocation } from 'mol-model/location';
import { createIdentityTransform } from '../../geometry/transform-data';
import { createRenderableState, updateRenderableState } from '../../geometry/geometry';
import { paramDefaultValues, NumberParam } from 'mol-view/parameter';
import { ValueCell } from 'mol-util';

export async function createVolumeSurface(ctx: RuntimeContext, volume: VolumeData, isoValue: VolumeIsoValue, mesh?: Mesh) {
    ctx.update({ message: 'Marching cubes...' });

    const surface = await computeMarchingCubesMesh({
        isoLevel: VolumeIsoValue.toAbsolute(isoValue).absoluteValue,
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
    isoValue: NumberParam('Iso Value', '', 2, -15, 15, 0.01),
}
export const DefaultIsosurfaceProps = paramDefaultValues(IsosurfaceParams)
export type IsosurfaceProps = typeof DefaultIsosurfaceProps

export function IsosurfaceVisual(): VolumeVisual<IsosurfaceProps> {
    let currentProps = DefaultIsosurfaceProps
    let renderObject: MeshRenderObject
    let currentVolume: VolumeData
    let mesh: Mesh

    async function create(ctx: RuntimeContext, volume: VolumeData, props: Partial<IsosurfaceProps> = {}) {
        currentProps = { ...DefaultIsosurfaceProps, ...props }

        mesh = await createVolumeSurface(ctx, volume,  VolumeIsoValue.relative(volume.dataStats, currentProps.isoValue))

        const locationIt = LocationIterator(1, 1, () => NullLocation)
        const transform = createIdentityTransform()

        const values = await Mesh.createValues(ctx, mesh, transform, locationIt, currentProps)
        const state = createRenderableState(currentProps)

        renderObject = createMeshRenderObject(values, state)
    }

    async function update(ctx: RuntimeContext, props: Partial<IsosurfaceProps> = {}) {
        const newProps = { ...currentProps, ...props }

        let createMesh = false

        if (newProps.isoValue !== currentProps.isoValue) createMesh = true

        if (createMesh) {
            mesh = await createVolumeSurface(ctx, currentVolume,  VolumeIsoValue.relative(currentVolume.dataStats, currentProps.isoValue), mesh)
            ValueCell.update(renderObject.values.drawCount, mesh.triangleCount * 3)
        }

        Mesh.updateValues(renderObject.values, newProps)
        updateRenderableState(renderObject.state, newProps)

        currentProps = newProps
    }

    return {
        get renderObject () { return renderObject },
        async createOrUpdate(ctx: RuntimeContext, props: Partial<IsosurfaceProps> = {}, volume?: VolumeData) {
            if (!volume && !currentVolume) {
                throw new Error('missing volume')
            } else if (volume && (!currentVolume || !renderObject)) {
                currentVolume = volume
                await create(ctx, volume, props)
            } else if (volume && volume !== currentVolume) {
                currentVolume = volume
                await create(ctx, volume, props)
            } else {
                await update(ctx, props)
            }

            currentProps = { ...DefaultIsosurfaceProps, ...props }
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