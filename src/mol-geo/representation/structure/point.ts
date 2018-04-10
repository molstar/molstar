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
import { Element, Unit, ElementSet } from 'mol-model/structure';
import P from 'mol-model/structure/query/properties';
import { RepresentationProps, UnitRepresentation } from './index';
import { Task } from 'mol-task'


export default function Point(): UnitRepresentation {
    const renderObjects: RenderObject[] = []
    const vertices = ChunkedArray.create(Float32Array, 3, 1024, 2048);

    return {
        create: (units: ReadonlyArray<Unit>, elements: ElementSet, props: Partial<RepresentationProps> = {}) => Task.create('Spacefill', async ctx => {
            const l = Element.Location();

            const unitIds = ElementSet.unitIds(elements);
            for (let i = 0, _i = unitIds.length; i < _i; i++) {
                const unitId = unitIds[i];
                const unit = units[unitId];
                const elementGroup = ElementSet.unitGetByIndex(elements, i);
                const elementCount = OrderedSet.size(elementGroup.elements)
                l.unit = unit;

                for (let i = 0; i < elementCount; i++) {
                    l.element = OrderedSet.getAt(elementGroup.elements, i)
                    ChunkedArray.add3(vertices, P.atom.x(l), P.atom.y(l), P.atom.z(l))
                }

                if (i % 10 === 0 && ctx.shouldUpdate) {
                    await ctx.update({ message: 'Point', current: i, max: _i });
                }
            }

            const transformArray = new Float32Array(32)
            const m4 = Mat4.identity()
            Mat4.toArray(m4, transformArray, 0)

            const color = ValueCell.create(createColorTexture(1))
            color.ref.value.set([ 0, 0, 255 ])

            const points = createRenderObject('point', {
                position: ValueCell.create(ChunkedArray.compact(vertices, true) as Float32Array),
                color,
                transform: ValueCell.create(transformArray),

                instanceCount: transformArray.length / 16,
                positionCount: vertices.elementCount
            }, {})
            renderObjects.push(points)

            return renderObjects
        }),
        update: (props: RepresentationProps) => false
    }
}
