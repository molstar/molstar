/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from 'mol-util/value-cell'

import { createPointRenderObject, RenderObject } from 'mol-gl/scene'
import { Mat4 } from 'mol-math/linear-algebra'
import { OrderedSet } from 'mol-data/int'
import { ChunkedArray } from 'mol-data/util';
import { Unit, ElementGroup } from 'mol-model/structure';
import { RepresentationProps, UnitsRepresentation } from './index';
import { Task } from 'mol-task'
import { VdwRadius } from 'mol-model/structure/model/properties/atomic';
import { createUniformColor, createInstanceColor } from '../../color/data';
import { fillSerial } from 'mol-gl/renderable/util';
import { ColorScale } from '../../color/scale';
import { createUniformSize } from '../../size/data';
import { vdwSizeData } from '../../size/structure/vdw';

export const DefaultPointProps = {

}
export type PointProps = Partial<typeof DefaultPointProps>

export default function Point(): UnitsRepresentation<PointProps> {
    const renderObjects: RenderObject[] = []
    const vertices = ChunkedArray.create(Float32Array, 3, 1024, 2048);
    const sizes = ChunkedArray.create(Float32Array, 1, 1024, 2048);

    return {
        renderObjects,
        create: (units: ReadonlyArray<Unit>, elementGroup: ElementGroup, props: PointProps = {}) => Task.create('Spacefill', async ctx => {
            // const l = Element.Location();

            const { x, y, z } = units[0].model.conformation
            const { type_symbol } = units[0].model.hierarchy.atoms
            const elementCount = OrderedSet.size(elementGroup.elements)
            for (let i = 0; i < elementCount; i++) {
                const e = OrderedSet.getAt(elementGroup.elements, i)
                ChunkedArray.add3(vertices, x[e], y[e], z[e])
                ChunkedArray.add(sizes, VdwRadius(type_symbol.value(e)))

                if (i % 10000 === 0 && ctx.shouldUpdate) {
                    await ctx.update({ message: 'Point', current: i, max: elementCount });
                }
            }

            const unitCount = units.length
            const transformArray = new Float32Array(unitCount * 16)
            for (let i = 0; i < unitCount; i++) {
                Mat4.toArray(units[i].operator.matrix, transformArray, i * 16)
            }

            // const color = createUniformColor({ value: 0xFF4411 })
            const scale = ColorScale.create({ domain: [ 0, unitCount - 1 ] })
            const color = createInstanceColor({
                colorFn: scale.color,
                unitCount
            })

            // const size = createUniformSize({ value: 1 })
            const size = vdwSizeData({
                units,
                elementGroup,
                offsetData: {
                    primitiveCount: elementCount,
                    offsetCount: elementCount + 1,
                    offsets: fillSerial(new Uint32Array(elementCount + 1))
                }
            })
            console.log(size)

            const points = createPointRenderObject({
                objectId: 0,

                position: ValueCell.create(ChunkedArray.compact(vertices, true) as Float32Array),
                id: ValueCell.create(fillSerial(new Float32Array(unitCount))),
                size,
                color,
                transform: ValueCell.create(transformArray),

                instanceCount: unitCount,
                elementCount,
                positionCount: vertices.elementCount,

                usePointSizeAttenuation: true
            })
            renderObjects.push(points)
        }),
        update: (props: RepresentationProps) => false
    }
}
