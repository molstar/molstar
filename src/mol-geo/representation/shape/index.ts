/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Task } from 'mol-task'
import { RenderObject, createMeshRenderObject, MeshRenderObject } from 'mol-gl/render-object';
import { RepresentationProps, Representation } from '..';
import { PickingId } from '../../geometry/picking';
import { Loci, EmptyLoci, isEveryLoci } from 'mol-model/loci';
import { MarkerAction, applyMarkerAction } from '../../geometry/marker-data';
import { ValueCell } from 'mol-util';
import { ColorThemeName, ColorThemeOptions } from 'mol-view/theme/color';
import { Shape } from 'mol-model/shape';
import { LocationIterator } from '../../util/location-iterator';
import { OrderedSet, Interval } from 'mol-data/int';
import { createIdentityTransform } from '../../geometry/transform-data';
import { createRenderableState } from '../../geometry/geometry';
import { Mesh } from '../../geometry/mesh/mesh';
import { paramDefaultValues, SelectParam } from 'mol-view/parameter';

export interface ShapeRepresentation<P extends RepresentationProps = {}> extends Representation<Shape, P> { }

export const ShapeParams = {
    ...Mesh.Params,
    colorTheme: SelectParam<ColorThemeName>('Color Theme', '', 'shape-group', ColorThemeOptions)
}
export const DefaultShapeProps = paramDefaultValues(ShapeParams)
export type ShapeProps = typeof DefaultShapeProps

// TODO
// export type ShapeRepresentation = ShapeRepresentation<ShapeProps>

export function ShapeRepresentation<P extends ShapeProps>(): ShapeRepresentation<P> {
    const renderObjects: RenderObject[] = []
    let _renderObject: MeshRenderObject | undefined
    let _shape: Shape
    let currentProps: P

    function createOrUpdate(props: Partial<P> = {}, shape?: Shape) {
        currentProps = Object.assign({}, DefaultShapeProps, currentProps, props)
        if (shape) _shape = shape

        return Task.create('ShapeRepresentation.create', async ctx => {
            renderObjects.length = 0

            if (!_shape) return

            const mesh = _shape.mesh
            const locationIt = ShapeGroupIterator.fromShape(_shape)
            const transform = createIdentityTransform()

            const values = await Mesh.createValues(ctx, mesh, transform, locationIt, currentProps)
            const state = createRenderableState(currentProps)

            _renderObject = createMeshRenderObject(values, state)
            renderObjects.push(_renderObject)
        });
    }

    return {
        label: 'Shape mesh',
        params: ShapeParams,
        get renderObjects () { return renderObjects },
        get props () { return currentProps },
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