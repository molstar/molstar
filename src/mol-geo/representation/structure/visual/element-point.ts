/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ValueCell } from 'mol-util/value-cell'
import { createPointRenderObject, PointRenderObject } from 'mol-gl/render-object'
import { Unit } from 'mol-model/structure';
import { RuntimeContext } from 'mol-task'

import { UnitsVisual, DefaultStructureProps } from '..';
import { getElementLoci, StructureElementIterator } from './util/element';
import { createTransforms, createColors, createSizes } from './util/common';
import { deepEqual, defaults } from 'mol-util';
import { SortedArray } from 'mol-data/int';
import { RenderableState, PointValues } from 'mol-gl/renderable';
import { PickingId } from '../../../util/picking';
import { Loci, EmptyLoci } from 'mol-model/loci';
import { MarkerAction, createMarkers } from '../../../util/marker-data';
import { Vec3 } from 'mol-math/linear-algebra';
import { fillSerial } from 'mol-util/array';
import { SizeThemeProps } from 'mol-view/theme/size';
import { LocationIterator } from '../../../util/location-iterator';

export const DefaultElementPointProps = {
    ...DefaultStructureProps,
    sizeTheme: { name: 'physical' } as SizeThemeProps
}
export type ElementPointProps = Partial<typeof DefaultElementPointProps>

// TODO make async
export function createElementPointVertices(unit: Unit, vertices?: ValueCell<Float32Array>) {
    const elements = unit.elements
    const n = elements.length * 3
    const array = vertices && vertices.ref.value.length >= n ? vertices.ref.value : new Float32Array(n)

    const pos = unit.conformation.invariantPosition

    const p = Vec3.zero()
    for (let i = 0; i < n; i += 3) {
        pos(elements[i / 3], p)
        array[i] = p[0]
        array[i + 1] = p[1]
        array[i + 2] = p[2]
    }
    return vertices ? ValueCell.update(vertices, array) : ValueCell.create(array)
}

export function ElementPointVisual(): UnitsVisual<ElementPointProps> {
    let renderObject: PointRenderObject | undefined
    let currentProps = DefaultElementPointProps
    let currentGroup: Unit.SymmetryGroup
    let locationIt: LocationIterator

    let _units: ReadonlyArray<Unit>
    let _elements: SortedArray

    return {
        get renderObject () { return renderObject },
        async createOrUpdate(ctx: RuntimeContext, props: ElementPointProps = {}, group?: Unit.SymmetryGroup) {
            if (!group && !currentGroup) {
                throw new Error('missing group')
            } else if (group && !currentGroup) {
                currentProps = Object.assign({}, DefaultElementPointProps, props)
                currentGroup = group
                locationIt = StructureElementIterator.fromGroup(group)

                _units = group.units
                _elements = group.elements;

                const { colorTheme, sizeTheme } = currentProps
                const elementCount = _elements.length
                const instanceCount = group.units.length

                const vertices = createElementPointVertices(_units[0])
                const transform = createTransforms(group)
                // console.time('createColors point')
                const color = await createColors(ctx, locationIt, colorTheme)
                // console.timeEnd('createColors point')
                const size = createSizes(locationIt, sizeTheme)
                const marker = createMarkers(instanceCount * elementCount)

                const values: PointValues = {
                    aPosition: vertices,
                    aGroup: ValueCell.create(fillSerial(new Float32Array(elementCount))),
                    aInstance: ValueCell.create(fillSerial(new Float32Array(instanceCount))),
                    ...transform,
                    ...color,
                    ...marker,
                    ...size,

                    uAlpha: ValueCell.create(defaults(props.alpha, 1.0)),
                    uInstanceCount: ValueCell.create(instanceCount),
                    uGroupCount: ValueCell.create(group.elements.length),

                    drawCount: ValueCell.create(group.elements.length),
                    instanceCount: ValueCell.create(instanceCount),

                    dPointSizeAttenuation: ValueCell.create(true),
                    dUseFog: ValueCell.create(defaults(props.useFog, true)),
                }
                const state: RenderableState = {
                    depthMask: defaults(props.depthMask, true),
                    visible: defaults(props.visible, true)
                }

                renderObject = createPointRenderObject(values, state)
            } else if (renderObject) {
                if (group) currentGroup = group

                const newProps = { ...currentProps, ...props }

                let updateTransform = false
                let createVertices = false
                let updateColor = false
                let updateSize = false

                if (currentGroup.units.length !== locationIt.instanceCount) updateTransform = true
                if (!deepEqual(newProps.sizeTheme, currentProps.sizeTheme)) createVertices = true
                if (!deepEqual(newProps.colorTheme, currentProps.colorTheme)) updateColor = true
                if (!deepEqual(newProps.sizeTheme, currentProps.sizeTheme)) updateSize = true

                if (updateTransform) {
                    locationIt = StructureElementIterator.fromGroup(currentGroup)
                    const { instanceCount, groupCount } = locationIt
                    createTransforms(currentGroup, renderObject.values)
                    createMarkers(instanceCount * groupCount, renderObject.values)
                    ValueCell.update(renderObject.values.instanceCount, instanceCount)
                    ValueCell.update(renderObject.values.aInstance, fillSerial(new Float32Array(instanceCount))) // TODO reuse array
                    updateColor = true
                    updateSize = true
                }

                if (createVertices) {
                    createElementPointVertices(currentGroup.units[0], renderObject.values.aPosition)
                    ValueCell.update(renderObject.values.aGroup, fillSerial(new Float32Array(locationIt.groupCount))) // TODO reuse array
                    ValueCell.update(renderObject.values.drawCount, currentGroup.elements.length)
                    updateColor = true
                    updateSize = true
                }

                if (updateColor) {
                    await createColors(ctx, locationIt, newProps.colorTheme, renderObject.values)
                }

                if (updateSize) {
                    createSizes(locationIt, newProps.sizeTheme, renderObject.values)
                }

                currentProps = newProps
            }
        },
        getLoci(pickingId: PickingId) {
            return renderObject ? getElementLoci(pickingId, currentGroup, renderObject.id) : EmptyLoci
        },
        mark(loci: Loci, action: MarkerAction) {
            // TODO
            // markElement(loci, action, currentGroup, renderObject.values)
        },
        destroy() {
            // TODO
            renderObject = undefined
        }
    }
}
