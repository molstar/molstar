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
import { SizeThemeProps } from '../../../theme';
import { getElementLoci } from './util/element';
import { createTransforms, createColors, createSizes } from './util/common';
import { deepEqual, defaults } from 'mol-util';
import { SortedArray } from 'mol-data/int';
import { RenderableState, PointValues } from 'mol-gl/renderable';
import { PickingId } from '../../../util/picking';
import { Loci } from 'mol-model/loci';
import { MarkerAction, createMarkers } from '../../../util/marker-data';
import { Vec3 } from 'mol-math/linear-algebra';
import { fillSerial } from 'mol-util/array';
import { StructureElementIterator } from './util/location-iterator';

export const DefaultPointProps = {
    ...DefaultStructureProps,
    sizeTheme: { name: 'physical' } as SizeThemeProps
}
export type PointProps = Partial<typeof DefaultPointProps>

export function createPointVertices(unit: Unit) {
    const elements = unit.elements
    const elementCount = elements.length
    const vertices = new Float32Array(elementCount * 3)

    const pos = unit.conformation.invariantPosition

    const p = Vec3.zero()
    for (let i = 0; i < elementCount; i++) {
        const i3 = i * 3
        pos(elements[i], p)
        vertices[i3] = p[0]
        vertices[i3 + 1] = p[1]
        vertices[i3 + 2] = p[2]
    }
    return vertices
}

export default function PointVisual(): UnitsVisual<PointProps> {
    let renderObject: PointRenderObject
    let currentProps = DefaultPointProps
    let currentGroup: Unit.SymmetryGroup

    let _units: ReadonlyArray<Unit>
    let _elements: SortedArray

    return {
        get renderObject () { return renderObject },
        async create(ctx: RuntimeContext, group: Unit.SymmetryGroup, props: PointProps = {}) {
            currentProps = Object.assign({}, DefaultPointProps, props)
            currentGroup = group

            _units = group.units
            _elements = group.elements;

            const { colorTheme, sizeTheme } = currentProps
            const elementCount = _elements.length
            const instanceCount = group.units.length

            const locationIt = StructureElementIterator.fromGroup(group)

            const vertices = createPointVertices(_units[0])
            const transforms = createTransforms(group)
            const color = createColors(locationIt, colorTheme)
            const size = createSizes(locationIt, sizeTheme)
            const marker = createMarkers(instanceCount * elementCount)

            const values: PointValues = {
                aPosition: ValueCell.create(vertices),
                aElementId: ValueCell.create(fillSerial(new Float32Array(elementCount))),
                aTransform: transforms,
                aInstanceId: ValueCell.create(fillSerial(new Float32Array(instanceCount))),
                ...color,
                ...marker,
                ...size,

                uAlpha: ValueCell.create(defaults(props.alpha, 1.0)),
                uInstanceCount: ValueCell.create(instanceCount),
                uElementCount: ValueCell.create(group.elements.length),

                drawCount: ValueCell.create(vertices.length / 3),
                instanceCount: ValueCell.create(instanceCount),

                dPointSizeAttenuation: ValueCell.create(true),
                dUseFog: ValueCell.create(defaults(props.useFog, true)),
            }
            const state: RenderableState = {
                depthMask: defaults(props.depthMask, true),
                visible: defaults(props.visible, true)
            }

            renderObject = createPointRenderObject(values, state)
        },
        async update(ctx: RuntimeContext, props: PointProps) {
            if (!renderObject || !_units || !_elements) return false

            const newProps = { ...currentProps, ...props }
            if (deepEqual(currentProps, newProps)) {
                console.log('props identical, nothing to change')
                return true
            }

            if (!deepEqual(currentProps.colorTheme, newProps.colorTheme)) {
                console.log('colorTheme changed', currentProps.colorTheme, newProps.colorTheme)
            }

            if (!deepEqual(currentProps.sizeTheme, newProps.sizeTheme)) {
                console.log('sizeTheme changed', currentProps.sizeTheme, newProps.sizeTheme)
            }

            currentProps = newProps
            return false
        },
        getLoci(pickingId: PickingId) {
            return getElementLoci(pickingId, currentGroup, renderObject.id)
        },
        mark(loci: Loci, action: MarkerAction) {
            // TODO
            // markElement(loci, action, currentGroup, renderObject.values)
        },
        destroy() {
            // TODO
        }
    }
}
