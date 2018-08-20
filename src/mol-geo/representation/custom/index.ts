/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Task } from 'mol-task'
import { RenderObject, createMeshRenderObject, MeshRenderObject } from 'mol-gl/render-object';
import { RepresentationProps, Representation } from '..';
import { PickingId } from '../../util/picking';
import { Loci, EmptyLoci, CustomLoci, isCustomLoci, isEveryLoci } from 'mol-model/loci';
import { MarkerAction, applyMarkerAction, createMarkers } from '../../util/marker-data';
import { createRenderableState, createMeshValues, createIdentityTransform, DefaultMeshProps } from '../util';
import { Mesh } from '../../shape/mesh';
import { getMeshData } from '../../util/mesh-data';
import { MeshValues } from 'mol-gl/renderable';
import { createValueColor } from '../../util/color-data';
import { Color } from 'mol-util/color';
import { CustomLocation } from 'mol-model/location';
import { ValueCell } from 'mol-util';

export interface CustomRepresentation<P extends RepresentationProps = {}> extends Representation<Mesh, P> { }

export const DefaultCustomProps = {
    ...DefaultMeshProps,
}
export type CustomProps = typeof DefaultCustomProps

export function CustomRepresentation<P extends CustomProps>(): CustomRepresentation<P> {
    const renderObjects: RenderObject[] = []
    let _renderObject: MeshRenderObject
    let _mesh: Mesh
    let _props: P

    function create(mesh: Mesh, props: Partial<P> = {}) {
        _props = Object.assign({}, DefaultCustomProps, _props, props)
        return Task.create('CustomRepresentation.create', async ctx => {
            renderObjects.length = 0
            _mesh = mesh

            const elementCount = mesh.triangleCount
            const instanceCount = 1

            const color = createValueColor(Color(0x7ec0ee))
            const marker = createMarkers(instanceCount * elementCount)
            const counts = { drawCount: mesh.triangleCount * 3, elementCount, instanceCount }

            const values: MeshValues = {
                ...getMeshData(mesh),
                ...createMeshValues(_props, counts),
                aTransform: createIdentityTransform(),
                ...color,
                ...marker,

                elements: mesh.indexBuffer,
            }
            const state = createRenderableState(_props)

            _renderObject = createMeshRenderObject(values, state)
            renderObjects.push(_renderObject)
        });
    }

    function update(props: Partial<P>) {
        return Task.create('CustomRepresentation.update', async ctx => {
            // TODO
        })
    }

    return {
        get renderObjects () { return renderObjects },
        get props () { return _props },
        create,
        update,
        getLoci(pickingId: PickingId) {
            const { objectId, elementId } = pickingId
            if (_renderObject.id === objectId) {
                return CustomLoci([ CustomLocation(_mesh, elementId) ])
            }
            return EmptyLoci
        },
        mark(loci: Loci, action: MarkerAction) {
            const { tMarker } = _renderObject.values
            let changed = false
            if (isEveryLoci(loci)) {
                if (applyMarkerAction(tMarker.ref.value.array, 0, _mesh.triangleCount, action)) changed = true
                changed = true
            } else if (isCustomLoci(loci)) {
                for (const l of loci.locations) {
                    if (l.data === _mesh) {
                        if (applyMarkerAction(tMarker.ref.value.array, 0, _mesh.triangleCount, action)) changed = true
                        // TODO
                        // const idx = l.key
                        // if (idx !== undefined) {
                        //     if (applyMarkerAction(tMarker.ref.value.array, idx, idx + 1, action)) changed = true
                        // }
                    }
                }
            }
            if (changed) {
                ValueCell.update(tMarker, tMarker.ref.value)
            }
        },
        destroy() {
            // TODO
        }
    }
}