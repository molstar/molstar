/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ValueCell } from 'mol-util/value-cell'
import { createPointRenderObject, RenderObject, PointRenderObject } from 'mol-gl/render-object'
import { Unit, Element } from 'mol-model/structure';
import { Task } from 'mol-task'
import { fillSerial } from 'mol-gl/renderable/util';

import { UnitsRepresentation, DefaultStructureProps } from './index';
import VertexMap from '../../shape/vertex-map';
import { SizeTheme } from '../../theme';
import { createTransforms, createColors, createSizes, createFlags } from './utils';
import { deepEqual, defaults } from 'mol-util';
import { SortedArray } from 'mol-data/int';
import { RenderableState, PointValues } from 'mol-gl/renderable';
import { PickingId } from '../../util/picking';

export const DefaultPointProps = {
    ...DefaultStructureProps,
    sizeTheme: { name: 'vdw' } as SizeTheme
}
export type PointProps = Partial<typeof DefaultPointProps>

export function createPointVertices(unit: Unit) {
    const elements = unit.elements
    const elementCount = elements.length
    const vertices = new Float32Array(elementCount * 3)

    const { x, y, z } = unit.conformation
    const l = Element.Location()
    l.unit = unit

    for (let i = 0; i < elementCount; i++) {
        l.element = elements[i];
        const i3 = i * 3
        vertices[i3] = x(l.element)
        vertices[i3 + 1] = y(l.element)
        vertices[i3 + 2] = z(l.element)
    }
    return vertices
}

export default function Point(): UnitsRepresentation<PointProps> {
    const renderObjects: RenderObject[] = []
    let points: PointRenderObject
    let currentProps = DefaultPointProps
    let currentGroup: Unit.SymmetryGroup

    let _units: ReadonlyArray<Unit>
    let _elements: SortedArray

    return {
        renderObjects,
        create(group: Unit.SymmetryGroup, props: PointProps = {}) {
            currentProps = Object.assign({}, DefaultPointProps, props)

            return Task.create('Point.create', async ctx => {
                renderObjects.length = 0 // clear
                currentGroup = group

                _units = group.units
                _elements = group.elements;

                const { colorTheme, sizeTheme, hoverSelection } = currentProps
                const elementCount = _elements.length

                const vertexMap = VertexMap.create(
                    elementCount,
                    elementCount + 1,
                    fillSerial(new Uint32Array(elementCount)),
                    fillSerial(new Uint32Array(elementCount + 1))
                )

                await ctx.update('Computing point vertices');
                const vertices = createPointVertices(_units[0])

                await ctx.update('Computing point transforms');
                const transforms = createTransforms(group)

                await ctx.update('Computing point colors');
                const color = createColors(group, vertexMap, colorTheme)

                await ctx.update('Computing point sizes');
                const size = createSizes(group, vertexMap, sizeTheme)

                await ctx.update('Computing spacefill flags');
                const flag = createFlags(group, hoverSelection.instanceId, hoverSelection.elementId)

                const instanceCount = group.units.length

                const values: PointValues = {
                    aPosition: ValueCell.create(vertices),
                    aElementId: ValueCell.create(fillSerial(new Float32Array(elementCount))),
                    aTransform: transforms,
                    aInstanceId: ValueCell.create(fillSerial(new Float32Array(instanceCount))),
                    ...color,
                    ...flag,
                    ...size,

                    uAlpha: ValueCell.create(defaults(props.alpha, 1.0)),
                    uInstanceCount: ValueCell.create(instanceCount),
                    uElementCount: ValueCell.create(group.elements.length),

                    drawCount: ValueCell.create(vertices.length / 3),
                    instanceCount: ValueCell.create(instanceCount),

                    dPointSizeAttenuation: ValueCell.create(true)
                }
                const state: RenderableState = {
                    depthMask: defaults(props.depthMask, true),
                    visible: defaults(props.visible, true)
                }

                points = createPointRenderObject(values, state)
                renderObjects.push(points)
            })
        },
        update(props: PointProps) {
            return Task.create('Point.update', async ctx => {
                if (!points || !_units || !_elements) return false

                const newProps = { ...currentProps, ...props }
                if (deepEqual(currentProps, newProps)) {
                    console.log('props identical, nothing to change')
                    return true
                }

                // const elementCount = OrderedSet.size(_elementGroup.elements)
                // const unitCount = _units.length

                // const vertexMap = VertexMap.create(
                //     elementCount,
                //     elementCount + 1,
                //     fillSerial(new Uint32Array(elementCount)),
                //     fillSerial(new Uint32Array(elementCount + 1))
                // )

                if (!deepEqual(currentProps.colorTheme, newProps.colorTheme)) {
                    console.log('colorTheme changed', currentProps.colorTheme, newProps.colorTheme)
                    // await ctx.update('Computing point colors');
                    // const color = createColors(_units, _elementGroup, vertexMap, newProps.colorTheme)
                    // ValueCell.update(points.props.color, color)
                }

                if (!deepEqual(currentProps.sizeTheme, newProps.sizeTheme)) {
                    console.log('sizeTheme changed', currentProps.sizeTheme, newProps.sizeTheme)
                }

                currentProps = newProps
                return false
            })
        },
        getLoci(pickingId: PickingId) {
            const { objectId, instanceId, elementId } = pickingId
            if (points.id === objectId) {
                const unit = currentGroup.units[instanceId]
                const elements = SortedArray.ofSingleton(currentGroup.elements[elementId])
                return Element.Loci([{ unit, elements }])
            }
            return null
        }
    }
}
