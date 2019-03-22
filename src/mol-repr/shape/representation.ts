/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Task, RuntimeContext } from 'mol-task'
import { createRenderObject, GraphicsRenderObject } from 'mol-gl/render-object';
import { Representation } from '../representation';
import { Loci, EmptyLoci, isEveryLoci } from 'mol-model/loci';
import { ValueCell } from 'mol-util';
import { Shape, ShapeGroup } from 'mol-model/shape';
import { OrderedSet, Interval } from 'mol-data/int';
import { ParamDefinition as PD } from 'mol-util/param-definition';
import { createTransform, TransformData } from 'mol-geo/geometry/transform-data';
import { PickingId } from 'mol-geo/geometry/picking';
import { MarkerAction, createMarkers } from 'mol-geo/geometry/marker-data';
import { LocationIterator } from 'mol-geo/util/location-iterator';
import { createEmptyTheme, Theme } from 'mol-theme/theme';
import { Subject } from 'rxjs';
import { Geometry, GeometryUtils } from 'mol-geo/geometry/geometry';
import { ShapeGroupColorTheme } from 'mol-theme/color/shape-group';
import { createColors } from 'mol-geo/geometry/color-data';
import { VisualUpdateState } from 'mol-repr/util';
import { Mat4 } from 'mol-math/linear-algebra';
import { Visual } from 'mol-repr/visual';
import { createSizes } from 'mol-geo/geometry/size-data';
import { ShapeGroupSizeTheme } from 'mol-theme/size/shape-group';

export interface ShapeRepresentation<D, G extends Geometry, P extends Geometry.Params<G>> extends Representation<D, P> { }

export type ShapeGetter<D, G extends Geometry, P extends Geometry.Params<G>> = (ctx: RuntimeContext, data: D, props: PD.Values<P>, shape?: Shape<G>) => Shape<G> | Promise<Shape<G>>

export function ShapeRepresentation<D, G extends Geometry, P extends Geometry.Params<G>>(getShape: ShapeGetter<D, G, P>, geometryUtils: GeometryUtils<G>): ShapeRepresentation<D, G, P> {
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

    const updateState = VisualUpdateState.create()

    function prepareUpdate(props: Partial<PD.Values<P>> = {}, shape?: Shape<G>) {
        VisualUpdateState.reset(updateState)

        if (!shape && !_shape) {
            // console.error('no shape given')
            return
        } else if (shape && !_shape) {
            // console.log('first shape')
            updateState.createNew = true
        } else if (shape && _shape && shape.id === _shape.id) {
            // console.log('same shape')
            // nothing to set
        } else if (shape && _shape && shape.id !== _shape.id) {
            // console.log('new shape')
            updateState.updateTransform = true
            updateState.createGeometry = true
        } else if (!shape) {
            // console.log('only props')
            // nothing to set
        } else {
            console.warn('unexpected state')
        }

        if (updateState.updateTransform) {
            updateState.updateColor = true
            updateState.updateSize = true
        }

        if (updateState.createGeometry) {
            updateState.updateColor = true
            updateState.updateSize = true
        }
    }

    function createOrUpdate(props: Partial<PD.Values<P>> = {}, data?: D) {
        return Task.create('ShapeRepresentation.create', async runtime => {
            const newProps = Object.assign(currentProps, props)
            const shape = data ? await getShape(runtime, data, newProps, _shape) : undefined

            prepareUpdate(props, shape)

            if (shape) {
                _shape = shape
                _theme.color = ShapeGroupColorTheme({ shape: _shape }, {})
                _theme.size = ShapeGroupSizeTheme({ shape: _shape }, {})
            }

            if (updateState.createNew) {
                renderObjects.length = 0 // clear list o renderObjects
                locationIt = ShapeGroupIterator.fromShape(_shape)
                const transform = createShapeTransform(_shape.transforms)
                const values = geometryUtils.createValues(_shape.geometry, transform, locationIt, _theme, newProps)
                const state = geometryUtils.createRenderableState(newProps)

                _renderObject = createRenderObject(_shape.geometry.kind, values, state)
                if (_renderObject) renderObjects.push(_renderObject) // add new renderObject to list
            } else {
                if (!_renderObject) {
                    throw new Error('expected renderObject to be available')
                }

                if (updateState.updateTransform) {
                    // console.log('update transform')
                    createShapeTransform(_shape.transforms, _renderObject.values)
                    locationIt = ShapeGroupIterator.fromShape(_shape)
                    const { instanceCount, groupCount } = locationIt
                    createMarkers(instanceCount * groupCount, _renderObject.values)
                }

                if (updateState.createGeometry) {
                    // console.log('update geometry')
                    ValueCell.update(_renderObject.values.drawCount, Geometry.getDrawCount(_shape.geometry))
                }

                if (updateState.updateTransform || updateState.createGeometry) {
                    // console.log('updateBoundingSphere')
                    geometryUtils.updateBoundingSphere(_renderObject.values, _shape.geometry)
                }

                if (updateState.updateColor) {
                    // console.log('update color')
                    createColors(locationIt, _theme.color, _renderObject.values)
                }

                if (updateState.updateSize) {
                    // not all geometries have size data, so check here
                    if ('uSize' in _renderObject.values) {
                        // console.log('update size')
                        createSizes(locationIt, _theme.size, _renderObject.values)
                    }
                }

                geometryUtils.updateValues(_renderObject.values, newProps)
                geometryUtils.updateRenderableState(_renderObject.state, newProps)
            }

            currentProps = newProps
            // increment version
            updated.next(version++)
        });
    }

    function lociApply(loci: Loci, apply: (interval: Interval) => boolean) {
        if (isEveryLoci(loci) || (Shape.isLoci(loci) && loci.shape === _shape)) {
            return apply(Interval.ofBounds(0, _shape.groupCount * _shape.transforms.length))
        } else {
            return eachShapeGroup(loci, _shape, apply)
        }
    }

    return {
        label: 'Shape geometry',
        get groupCount () { return locationIt ? locationIt.count : 0 },
        get props () { return currentProps },
        get params () { return currentParams },
        get state() { return _state },
        get theme() { return _theme },
        renderObjects,
        updated,
        createOrUpdate,
        getLoci(pickingId: PickingId) {
            const { objectId, groupId, instanceId } = pickingId
            if (_renderObject && _renderObject.id === objectId) {
                return ShapeGroup.Loci(_shape, [{ ids: OrderedSet.ofSingleton(groupId) }], instanceId)
            }
            return EmptyLoci
        },
        mark(loci: Loci, action: MarkerAction) {
            return Visual.mark(_renderObject, loci, action, lociApply)
        },
        setState(state: Partial<Representation.State>) {
            if (_renderObject) {
                if (state.visible !== undefined) Visual.setVisibility(_renderObject, state.visible)
                if (state.alphaFactor !== undefined) Visual.setAlphaFactor(_renderObject, state.alphaFactor)
                if (state.pickable !== undefined) Visual.setPickable(_renderObject, state.pickable)
                if (state.overpaint !== undefined) {
                    Visual.setOverpaint(_renderObject, state.overpaint, lociApply, true)
                }
                if (state.transform !== undefined) Visual.setTransform(_renderObject, state.transform)
            }

            Representation.updateState(_state, state)
        },
        setTheme(theme: Theme) {
            console.warn('The `ShapeRepresentation` theme is fixed to `ShapeGroupColorTheme` and `ShapeGroupSizeTheme`. Colors are taken from `Shape.getColor` and sizes from `Shape.getSize`')
        },
        destroy() {
            // TODO
            renderObjects.length = 0
            _renderObject = undefined
        }
    }
}

function createShapeTransform(transforms: Mat4[], transformData?: TransformData) {
    const transformArray = transformData && transformData.aTransform.ref.value.length >= transforms.length * 16 ? transformData.aTransform.ref.value : new Float32Array(transforms.length * 16)
    for (let i = 0, il = transforms.length; i < il; ++i) {
        Mat4.toArray(transforms[i], transformArray, i * 16)
    }
    return createTransform(transformArray, transforms.length, transformData)
}

function eachShapeGroup(loci: Loci, shape: Shape, apply: (interval: Interval) => boolean) {
    if (!ShapeGroup.isLoci(loci)) return false
    if (loci.shape !== shape) return false
    let changed = false
    const { groupCount } = shape
    const { instance, groups } = loci
    for (const g of groups) {
        if (Interval.is(g.ids)) {
            const start = instance * groupCount + Interval.start(g.ids)
            const end = instance * groupCount + Interval.end(g.ids)
            if (apply(Interval.ofBounds(start, end))) changed = true
        } else {
            for (let i = 0, _i = g.ids.length; i < _i; i++) {
                const idx = instance * groupCount + g.ids[i];
                if (apply(Interval.ofSingleton(idx))) changed = true
            }
        }
    }
    return changed
}

export namespace ShapeGroupIterator {
    export function fromShape(shape: Shape): LocationIterator {
        const instanceCount = shape.transforms.length
        const location = ShapeGroup.Location(shape)
        const getLocation = (groupIndex: number, instanceIndex: number) => {
            location.group = groupIndex
            location.instance = instanceIndex
            return location
        }
        return LocationIterator(shape.groupCount, instanceCount, getLocation)
    }
}