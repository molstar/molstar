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
import { getMeshData } from '../../util/mesh-data';
import { MeshValues } from 'mol-gl/renderable';
import { ValueCell } from 'mol-util';
import { ColorThemeProps } from 'mol-view/theme/color';
import { Shape } from 'mol-model/shape';
import { LocationIterator } from '../../util/location-iterator';
import { createColors } from '../structure/visual/util/common';
import { OrderedSet, Interval } from 'mol-data/int';
import { createIdentityTransform } from '../../util/transform-data';
import { DefaultMeshProps, createMeshValues, createRenderableState } from '../../geometry/geometry';

export interface ShapeRepresentation<P extends RepresentationProps = {}> extends Representation<Shape, P> { }

export const DefaultShapeProps = {
    ...DefaultMeshProps,

    colorTheme: { name: 'shape-group' } as ColorThemeProps
}
export type ShapeProps = typeof DefaultShapeProps

// TODO
// export type ShapeRepresentation = ShapeRepresentation<ShapeProps>

export function ShapeRepresentation<P extends ShapeProps>(): ShapeRepresentation<P> {
    const renderObjects: RenderObject[] = []
    let _renderObject: MeshRenderObject | undefined
    let _shape: Shape
    let _props: P

    function createOrUpdate(props: Partial<P> = {}, shape?: Shape) {
        _props = Object.assign({}, DefaultShapeProps, _props, props)
        if (shape) _shape = shape

        return Task.create('ShapeRepresentation.create', async ctx => {
            renderObjects.length = 0

            if (!_shape) return

            const mesh = _shape.mesh
            const locationIt = ShapeGroupIterator.fromShape(_shape)
            const { groupCount, instanceCount } = locationIt

            const transform = createIdentityTransform()
            const color = await createColors(ctx, locationIt, _props.colorTheme)
            const marker = createMarkers(instanceCount * groupCount)
            const counts = { drawCount: mesh.triangleCount * 3, groupCount, instanceCount }

            const values: MeshValues = {
                ...getMeshData(mesh),
                ...createMeshValues(_props, counts),
                ...transform,
                ...color,
                ...marker,

                elements: mesh.indexBuffer,
            }
            const state = createRenderableState(_props)

            _renderObject = createMeshRenderObject(values, state)
            renderObjects.push(_renderObject)
        });
    }

    return {
        label: 'Shape mesh',
        get renderObjects () { return renderObjects },
        get props () { return _props },
        createOrUpdate,
        getLoci(pickingId: PickingId) {
            const { objectId, groupId } = pickingId
            if (_renderObject && _renderObject.id === objectId) {
                return Shape.Loci([ { shape: _shape, ids: OrderedSet.ofSingleton(groupId) } ])
            }
            return EmptyLoci
        },
        mark(loci: Loci, action: MarkerAction) {
            if (!_renderObject) return false
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
            return changed
        },
        destroy() {
            // TODO
            renderObjects.length = 0
            _renderObject = undefined
        }
    }
}

export namespace ShapeGroupIterator {
    export function fromShape(shape: Shape): LocationIterator {
        const { groupCount } = shape
        const instanceCount = 1
        const location = Shape.Location(shape)
        const getLocation = (groupIndex: number) => {
            location.group = groupIndex
            return location
        }
        return LocationIterator(groupCount, instanceCount, getLocation)
    }
}