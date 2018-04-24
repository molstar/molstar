/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from 'mol-util/value-cell'

import { createPointRenderObject, RenderObject, PointRenderObject } from 'mol-gl/scene'

import { OrderedSet } from 'mol-data/int'
import { Unit, ElementGroup } from 'mol-model/structure';
import { RepresentationProps, UnitsRepresentation } from './index';
import { Task } from 'mol-task'
import { fillSerial } from 'mol-gl/renderable/util';

import VertexMap from '../../shape/vertex-map';
import { ColorTheme, SizeTheme } from '../../theme';
import { createTransforms, createColors, createSizes } from './utils';

export const DefaultPointProps = {
    colorTheme: { name: 'instance-index' } as ColorTheme,
    sizeTheme: { name: 'vdw' } as SizeTheme
}
export type PointProps = Partial<typeof DefaultPointProps>

export function createPointVertices(unit: Unit, elementGroup: ElementGroup) {
    const elementCount = OrderedSet.size(elementGroup.elements)
    const vertices = new Float32Array(elementCount * 3)
    const { x, y, z } = unit.model.conformation
    for (let i = 0; i < elementCount; i++) {
        const e = OrderedSet.getAt(elementGroup.elements, i)
        const i3 = i * 3
        vertices[i3] = x[e]
        vertices[i3 + 1] = y[e]
        vertices[i3 + 2] = z[e]
    }
    return vertices
}

export default function Point(): UnitsRepresentation<PointProps> {
    const renderObjects: RenderObject[] = []
    let points: PointRenderObject

    return {
        renderObjects,
        create(units: ReadonlyArray<Unit>, elementGroup: ElementGroup, props: PointProps = {}) {
            return Task.create('Point.create', async ctx => {
                renderObjects.length = 0 // clear

                const { colorTheme, sizeTheme } = { ...DefaultPointProps, ...props }
                const elementCount = OrderedSet.size(elementGroup.elements)
                const unitCount = units.length

                const vertexMap = VertexMap.create(
                    elementCount,
                    elementCount + 1,
                    fillSerial(new Uint32Array(elementCount)),
                    fillSerial(new Uint32Array(elementCount + 1))
                )

                await ctx.update('Computing point vertices');
                const vertices = createPointVertices(units[0], elementGroup)

                await ctx.update('Computing point transforms');
                const transforms = createTransforms(units)

                await ctx.update('Computing point colors');
                const color = createColors(units, elementGroup, vertexMap, colorTheme)

                await ctx.update('Computing point sizes');
                const size = createSizes(units, elementGroup, vertexMap, sizeTheme)

                points = createPointRenderObject({
                    objectId: 0,

                    position: ValueCell.create(vertices),
                    id: ValueCell.create(fillSerial(new Float32Array(unitCount))),
                    size,
                    color,
                    transform: ValueCell.create(transforms),

                    instanceCount: unitCount,
                    elementCount,
                    positionCount: vertices.length / 3,

                    usePointSizeAttenuation: true
                })
                renderObjects.push(points)
            })
        },
        update(props: RepresentationProps) {
            return Task.create('Point.update', async ctx => {
                if (!points) return false

                return false
            })
        }
    }
}
