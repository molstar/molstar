/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Task } from 'mol-task'
import { createRenderObject, GraphicsRenderObject } from 'mol-gl/render-object';
import { Representation } from '../representation';
import { Loci, EmptyLoci, isEveryLoci } from 'mol-model/loci';
import { ValueCell } from 'mol-util';
import { Shape } from 'mol-model/shape';
import { OrderedSet, Interval } from 'mol-data/int';
import { ParamDefinition as PD } from 'mol-util/param-definition';
import { createIdentityTransform } from 'mol-geo/geometry/transform-data';
import { PickingId } from 'mol-geo/geometry/picking';
import { MarkerAction, applyMarkerAction } from 'mol-geo/geometry/marker-data';
import { LocationIterator } from 'mol-geo/util/location-iterator';
import { createEmptyTheme, Theme } from 'mol-theme/theme';
import { Subject } from 'rxjs';
import { Geometry, GeometryUtils } from 'mol-geo/geometry/geometry';
import { ShapeGroupColorTheme } from 'mol-theme/color/shape-group';

export interface ShapeRepresentation<D, G extends Geometry, P extends Geometry.Params<G>> extends Representation<D, P> { }

export function ShapeRepresentation<D, G extends Geometry, P extends Geometry.Params<G>>(getShape: (data: D, props: PD.Values<P>, shape?: Shape<G>) => Shape<G>, geometryUtils: GeometryUtils<G>): ShapeRepresentation<D, G, P> {
    let version = 0
    const updated = new Subject<number>()
    const _state = Representation.createState()
    const renderObjects: GraphicsRenderObject[] = []
    let _renderObject: GraphicsRenderObject | undefined
    let _shape: Shape<G>
    let _theme = createEmptyTheme()
    let currentProps: PD.Values<P> = PD.getDefaultValues(geometryUtils.Params as P) // TODO avoid casting
    let currentParams: P
    let locationIt: LocationIterator

    function createOrUpdate(props: Partial<PD.Values<P>> = {}, data?: D) {
        currentProps = Object.assign(currentProps, props)
        const shape = data ? getShape(data, currentProps, _shape) : undefined

        return Task.create('ShapeRepresentation.create', async runtime => {
            if (!shape && !_shape) {
                console.error('no shape given')
                return
            } else if (shape && !_shape) {
                console.log('first shape')

            } else if (shape && _shape && shape.id === _shape.id) {
                console.log('same shape')

            } else if (shape && _shape && shape.id !== _shape.id) {
                console.log('new shape')

            } else {
                console.log('only props')

            }

            if (shape) _shape = shape
            renderObjects.length = 0
            locationIt = ShapeGroupIterator.fromShape(_shape)
            _theme.color = ShapeGroupColorTheme({ shape: _shape }, {})

            const transform = createIdentityTransform()
            const values = geometryUtils.createValues(_shape.geometry, transform, locationIt, _theme, currentProps)
            const state = geometryUtils.createRenderableState(currentProps)

            _renderObject = createRenderObject(_shape.geometry.kind, values, state)
            if (_renderObject) renderObjects.push(_renderObject)
            updated.next(version++)
        });
    }

    return {
        label: 'Shape geometry',
        get groupCount () { return locationIt ? locationIt.count : 0 },
        get renderObjects () { return renderObjects },
        get props () { return currentProps },
        get params () { return currentParams },
        get state() { return _state },
        get theme() { return _theme },
        updated,
        createOrUpdate,
        getLoci(pickingId: PickingId) {
            const { objectId, groupId } = pickingId
            if (_renderObject && _renderObject.id === objectId) {
                return Shape.Loci(_shape, [{ ids: OrderedSet.ofSingleton(groupId) }])
            }
            return EmptyLoci
        },
        mark(loci: Loci, action: MarkerAction) {
            if (!_renderObject) return false
            const { tMarker } = _renderObject.values
            let changed = false
            if (isEveryLoci(loci)) {
                if (applyMarkerAction(tMarker.ref.value.array, 0, _shape.groupCount, action)) changed = true
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
        setState(state: Partial<Representation.State>) {
            if (state.visible !== undefined) renderObjects.forEach(ro => ro.state.visible = !!state.visible)
            if (state.pickable !== undefined) renderObjects.forEach(ro => ro.state.pickable = !!state.pickable)
            // TODO state.transform

            Representation.updateState(_state, state)
        },
        setTheme(theme: Theme) {
            _theme = theme
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
        const instanceCount = 1
        const location = Shape.Location(shape)
        const getLocation = (groupIndex: number) => {
            location.group = groupIndex
            return location
        }
        return LocationIterator(shape.groupCount, instanceCount, getLocation)
    }
}