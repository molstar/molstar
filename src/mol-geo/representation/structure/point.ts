/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ValueCell } from 'mol-util/value-cell'
import { createPointRenderObject, RenderObject, PointRenderObject } from 'mol-gl/scene'
import { Unit, Element } from 'mol-model/structure';
import { Task } from 'mol-task'
import { fillSerial } from 'mol-gl/renderable/util';

import { UnitsRepresentation } from './index';
import VertexMap from '../../shape/vertex-map';
import { ColorTheme, SizeTheme } from '../../theme';
import { createTransforms, createColors, createSizes } from './utils';
import { deepEqual } from 'mol-util';
import { SortedArray } from 'mol-data/int';

export const DefaultPointProps = {
    colorTheme: { name: 'instance-index' } as ColorTheme,
    sizeTheme: { name: 'vdw' } as SizeTheme,
    alpha: 1,
    visible: true,
    depthMask: true
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
    let curProps = DefaultPointProps

    let _units: ReadonlyArray<Unit>
    let _elements: SortedArray

    return {
        renderObjects,
        create(group: Unit.SymmetryGroup, props: PointProps = {}) {
            return Task.create('Point.create', async ctx => {
                renderObjects.length = 0 // clear
                curProps = { ...DefaultPointProps, ...props }

                _units = group.units
                _elements = group.elements;

                const { colorTheme, sizeTheme, alpha, visible, depthMask } = curProps
                const elementCount = _elements.length
                const unitCount = _units.length

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

                points = createPointRenderObject({
                    objectId: 0,
                    alpha,
                    visible,
                    depthMask,

                    position: ValueCell.create(vertices),
                    id: ValueCell.create(fillSerial(new Float32Array(elementCount))),
                    size: size,
                    color: color,
                    transform: ValueCell.create(transforms),

                    instanceCount: unitCount,
                    elementCount,
                    positionCount: vertices.length / 3,

                    usePointSizeAttenuation: true
                })
                renderObjects.push(points)
            })
        },
        update(props: PointProps) {
            return Task.create('Point.update', async ctx => {
                if (!points || !_units || !_elements) return false

                const newProps = { ...curProps, ...props }
                if (deepEqual(curProps, newProps)) {
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

                if (!deepEqual(curProps.colorTheme, newProps.colorTheme)) {
                    console.log('colorTheme changed', curProps.colorTheme, newProps.colorTheme)
                    // await ctx.update('Computing point colors');
                    // const color = createColors(_units, _elementGroup, vertexMap, newProps.colorTheme)
                    // ValueCell.update(points.props.color, color)
                }

                if (!deepEqual(curProps.sizeTheme, newProps.sizeTheme)) {
                    console.log('sizeTheme changed', curProps.sizeTheme, newProps.sizeTheme)
                }

                curProps = newProps
                return false
            })
        }
    }
}
