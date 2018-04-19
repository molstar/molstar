/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from 'mol-util/value-cell'

import { createRenderObject, RenderObject } from 'mol-gl/scene'
// import { createColorTexture } from 'mol-gl/util';
import { Vec3, Mat4 } from 'mol-math/linear-algebra'
import { OrderedSet } from 'mol-data/int'
import { Unit, ElementGroup } from 'mol-model/structure';
import { RepresentationProps, UnitsRepresentation } from './index';
import { Task } from 'mol-task'
import { MeshBuilder } from '../../shape/mesh-builder';
import { VdwRadius } from 'mol-model/structure/model/properties/atomic';
import { ElementColor, colorToArray, normalizedColorToArray, ColorScale } from '../../color';
import { Color } from 'mol-gl/renderable/mesh';
import { createColorTexture } from 'mol-gl/util';

export const DefaultSpacefillProps = {
    detail: 0
}
export type SpacefillProps = Partial<typeof DefaultSpacefillProps>

// function buildColorBuffer() {
//     if (props && props.color) {
//         colors = new Float32Array(icosahedron.vertices.length)
//         for (let i = 0, il = colors.length; i < il; i += 3) {
//             hexColorToArray(props.color, colors, i)
//         }
//     }
// }

export default function Spacefill(): UnitsRepresentation<SpacefillProps> {
    const renderObjects: RenderObject[] = []

    return {
        renderObjects,
        create: (units: ReadonlyArray<Unit>, elementGroup: ElementGroup, props: SpacefillProps = {}) => Task.create('Spacefill', async ctx => {
            const { detail } = { ...DefaultSpacefillProps, ...props }
            const meshBuilder = MeshBuilder.create()

            const v = Vec3.zero()
            const m = Mat4.identity()

            const { x, y, z } = units[0].model.conformation
            const { type_symbol } = units[0].model.hierarchy.atoms
            const elementCount = OrderedSet.size(elementGroup.elements)
            for (let i = 0; i < elementCount; i++) {
                const e = OrderedSet.getAt(elementGroup.elements, i)
                v[0] = x[e]
                v[1] = y[e]
                v[2] = z[e]
                Mat4.setTranslation(m, v)

                meshBuilder.setId(i)
                meshBuilder.addIcosahedron(m, {
                    radius: VdwRadius(type_symbol.value(e)),
                    detail
                })

                if (i % 10 === 0 && ctx.shouldUpdate) {
                    await ctx.update({ message: 'Spacefill', current: i, max: elementCount });
                }
            }

            const mesh = meshBuilder.getMesh()
            // console.log(mesh)
            if (!mesh.offsetBuffer.ref.value) return

            const unitsCount = units.length
            const transformArray = new Float32Array(unitsCount * 16)
            for (let i = 0; i < unitsCount; i++) {
                Mat4.toArray(units[i].operator.matrix, transformArray, i * 16)
            }

            // console.log({ unitsCount, elementCount })

            let colorType = 'element'
            let color: Color

            if (colorType === 'attribute') {
                const colors = new Float32Array(mesh.vertexCount * 3);
                const offsets = mesh.offsetBuffer.ref.value
                for (let i = 0, il = mesh.offsetCount - 1; i < il; ++i) {
                    const start = offsets[i]
                    const end = offsets[i + 1]
                    const e = OrderedSet.getAt(elementGroup.elements, i)
                    const hexColor = ElementColor(type_symbol.value(e))
                    for (let i = start, il = end; i < il; ++i) {
                        normalizedColorToArray(hexColor, colors, i * 3)
                    }
                }
                color = { type: 'attribute', value: ValueCell.create(colors) }
            } else if (colorType === 'instance') {
                const colors = createColorTexture(unitsCount)
                const scale = ColorScale.create({ domain: [ 0, unitsCount - 1 ] })
                for (let i = 0; i < unitsCount; i++) {
                    scale.colorToArray(i, colors, i * 3)
                }
                color = { type: 'instance', value: ValueCell.create(colors) }
            } else if (colorType === 'element') {
                const elementCount = mesh.offsetCount - 1
                const count = unitsCount * elementCount
                const colors = createColorTexture(count)
                const scale = ColorScale.create({ domain: [ 0, count - 1 ] })
                let colorOffset = 0
                for (let i = 0; i < unitsCount; i++) {
                    for (let j = 0, jl = elementCount; j < jl; ++j) {
                        const hexColor = scale.color(i * elementCount + j)
                        colorToArray(hexColor, colors, colorOffset)
                        colorOffset += 3
                    }
                }
                color = { type: 'element', value: ValueCell.create(colors) }
            }

            const spheres = createRenderObject('mesh', {
                position: mesh.vertexBuffer,
                normal: mesh.normalBuffer,
                color: color!,
                id: mesh.idBuffer,
                transform: ValueCell.create(transformArray),
                index: mesh.indexBuffer,

                instanceCount: unitsCount,
                indexCount: mesh.triangleCount,
                elementCount: mesh.offsetCount - 1,
                positionCount: mesh.vertexCount
            }, {})
            renderObjects.push(spheres)
        }),
        update: (props: RepresentationProps) => false
    }
}
