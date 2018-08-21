/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Task } from 'mol-task'
import { RenderObject, createMeshRenderObject, MeshRenderObject } from 'mol-gl/render-object';
import { RepresentationProps, Representation } from '..';
import { PickingId } from '../../util/picking';
import { Loci, EmptyLoci, isEveryLoci } from 'mol-model/loci';
import { MarkerAction, applyMarkerAction, createMarkers } from '../../util/marker-data';
import { createRenderableState, createMeshValues, createIdentityTransform, DefaultMeshProps } from '../util';
import { getMeshData } from '../../util/mesh-data';
import { MeshValues } from 'mol-gl/renderable';
import { ValueCell } from 'mol-util';
import { ColorThemeProps } from 'mol-view/theme/color';
import { Shape } from 'mol-model/shape';
import { LocationIterator } from '../../util/location-iterator';
import { arrayMax } from 'mol-util/array';
import { createColors } from '../structure/visual/util/common';
import { OrderedSet, Interval } from 'mol-data/int';

export interface ShapeRepresentation<P extends RepresentationProps = {}> extends Representation<Shape, P> { }

export const DefaultShapeProps = {
    ...DefaultMeshProps,
    colorTheme: { name: 'shape-group' } as ColorThemeProps
}
export type ShapeProps = typeof DefaultShapeProps

export function ShapeRepresentation<P extends ShapeProps>(): ShapeRepresentation<P> {
    const renderObjects: RenderObject[] = []
    let _renderObject: MeshRenderObject
    let _shape: Shape
    let _props: P

    function create(shape: Shape, props: Partial<P> = {}) {
        _props = Object.assign({}, DefaultShapeProps, _props, props)
        _shape = shape

        return Task.create('ShapeRepresentation.create', async ctx => {
            renderObjects.length = 0

            const mesh = shape.mesh
            const locationIt = ShapeGroupIterator.fromShape(shape)
            const { groupCount, instanceCount } = locationIt

            const color = createColors(locationIt, _props.colorTheme)
            const marker = createMarkers(instanceCount * groupCount)
            const counts = { drawCount: mesh.triangleCount * 3, groupCount, instanceCount }

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
            console.log(_renderObject)
            renderObjects.push(_renderObject)
        });
    }

    function update(props: Partial<P>) {
        return Task.create('ShapeRepresentation.update', async ctx => {
            // TODO
        })
    }

    return {
        get renderObjects () { return renderObjects },
        get props () { return _props },
        create,
        update,
        getLoci(pickingId: PickingId) {
            const { objectId, groupId } = pickingId
            if (_renderObject.id === objectId) {
                return Shape.Loci([ { shape: _shape, ids: OrderedSet.ofSingleton(groupId) } ])
            }
            return EmptyLoci
        },
        mark(loci: Loci, action: MarkerAction) {
            const { tMarker } = _renderObject.values
            let changed = false
            if (isEveryLoci(loci)) {
                if (applyMarkerAction(tMarker.ref.value.array, 0, _shape.mesh.triangleCount, action)) changed = true
            } else if (Shape.isLoci(loci)) {
                for (const g of loci.groups) {
                    if (Interval.is(g.ids)) {
                        const start = Interval.start(g.ids)
                        const end = Interval.end(g.ids)
                        if (applyMarkerAction(tMarker.ref.value.array, start, end, action)) changed = true
                    } else {
                        for (let i = 0, _i = g.ids.length; i < _i; i++) {
                            const idx = g.ids[i];
                            if (applyMarkerAction(tMarker.ref.value.array, idx, idx + 1, action)) changed = true
                        }
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

export namespace ShapeGroupIterator {
    export function fromShape(shape: Shape): LocationIterator {
        const groupCount = arrayMax(shape.mesh.groupBuffer.ref.value) + 1
        const instanceCount = 1
        const location = Shape.Location(shape)
        const getLocation = (groupIndex: number) => {
            location.group = groupIndex
            return location
        }
        return LocationIterator(groupCount, instanceCount, getLocation)
    }
}