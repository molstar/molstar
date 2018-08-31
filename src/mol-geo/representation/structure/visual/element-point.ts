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
    let renderObject: PointRenderObject | undefined
    let currentProps = DefaultPointProps
    let currentGroup: Unit.SymmetryGroup
    let locationIt: LocationIterator

    let _units: ReadonlyArray<Unit>
    let _elements: SortedArray

    return {
        get renderObject () { return renderObject },
        async createOrUpdate(ctx: RuntimeContext, props: PointProps = {}, group?: Unit.SymmetryGroup) {
            if (!group && !currentGroup) {
                throw new Error('missing group')
            } else if (group && !currentGroup) {
                currentProps = Object.assign({}, DefaultPointProps, props)
                currentGroup = group
                locationIt = StructureElementIterator.fromGroup(group)

                _units = group.units
                _elements = group.elements;

                const { colorTheme, sizeTheme } = currentProps
                const elementCount = _elements.length
                const instanceCount = group.units.length

                const vertices = createPointVertices(_units[0])
                const transform = createTransforms(group)
                const color = await createColors(ctx, locationIt, colorTheme)
                const size = createSizes(locationIt, sizeTheme)
                const marker = createMarkers(instanceCount * elementCount)

                const values: PointValues = {
                    aPosition: ValueCell.create(vertices),
                    aGroup: ValueCell.create(fillSerial(new Float32Array(elementCount))),
                    aInstance: ValueCell.create(fillSerial(new Float32Array(instanceCount))),
                    ...transform,
                    ...color,
                    ...marker,
                    ...size,

                    uAlpha: ValueCell.create(defaults(props.alpha, 1.0)),
                    uInstanceCount: ValueCell.create(instanceCount),
                    uGroupCount: ValueCell.create(group.elements.length),

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
            } else if (renderObject) {
                const newProps = { ...currentProps, ...props }

                if (!deepEqual(currentProps.colorTheme, newProps.colorTheme)) {
                    await createColors(ctx, locationIt, newProps.colorTheme, renderObject.values)
                }

                if (!deepEqual(currentProps.sizeTheme, newProps.sizeTheme)) {
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
