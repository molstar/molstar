/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from 'mol-util/value-cell'

import { createRenderObject, RenderObject } from 'mol-gl/renderer'
import { createColorTexture } from 'mol-gl/util';
import { Mat4 } from 'mol-math/linear-algebra'
import { OrderedSet } from 'mol-data/int'
import { ChunkedArray } from 'mol-data/util';
import { Unit, ElementGroup } from 'mol-model/structure';
import { RepresentationProps, UnitsRepresentation } from './index';
import { Task } from 'mol-task'

export const DefaultPointProps = {

}
export type PointProps = Partial<typeof DefaultPointProps>

export default function Point(): UnitsRepresentation<PointProps> {
    const renderObjects: RenderObject[] = []
    const vertices = ChunkedArray.create(Float32Array, 3, 1024, 2048);

    return {
        renderObjects,
        create: (units: ReadonlyArray<Unit>, elementGroup: ElementGroup, props: PointProps = {}) => Task.create('Spacefill', async ctx => {
            // const l = Element.Location();

            const { x, y, z } = units[0].model.conformation
            const elementCount = OrderedSet.size(elementGroup.elements)
            for (let i = 0; i < elementCount; i++) {
                const e = OrderedSet.getAt(elementGroup.elements, i)
                ChunkedArray.add3(vertices, x[e], y[e], z[e])

                if (i % 10 === 0 && ctx.shouldUpdate) {
                    await ctx.update({ message: 'Point', current: i, max: elementCount });
                }
            }

            const unitsCount = units.length
            const transformArray = new Float32Array(unitsCount * 16)
            for (let i = 0; i < unitsCount; i++) {
                Mat4.toArray(units[i].operator.matrix, transformArray, i * 16)
            }

            const color = ValueCell.create(createColorTexture(unitsCount))
            color.ref.value.set([ 0, 0, 255 ])

            const points = createRenderObject('point', {
                position: ValueCell.create(ChunkedArray.compact(vertices, true) as Float32Array),
                color,
                transform: ValueCell.create(transformArray),

                instanceCount: unitsCount,
                positionCount: vertices.elementCount
            }, {})
            renderObjects.push(points)
        }),
        update: (props: RepresentationProps) => false
    }
}
