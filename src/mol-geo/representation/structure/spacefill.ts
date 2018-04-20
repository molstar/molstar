/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from 'mol-util/value-cell'

import { RenderObject, createMeshRenderObject } from 'mol-gl/scene'
// import { createColorTexture } from 'mol-gl/util';
import { Vec3, Mat4 } from 'mol-math/linear-algebra'
import { OrderedSet } from 'mol-data/int'
import { Unit, ElementGroup } from 'mol-model/structure';
import { RepresentationProps, UnitsRepresentation } from './index';
import { Task } from 'mol-task'
import { MeshBuilder } from '../../shape/mesh-builder';
import { VdwRadius } from 'mol-model/structure/model/properties/atomic';
import { elementSymbolColorData } from '../../color/structure/element';
import { ColorData } from '../../color';
import { createInstanceColor, createUniformColor, createElementInstanceColor } from '../../color/data';
import { ColorScale } from '../../color/scale';

export const DefaultSpacefillProps = {
    detail: 0
}
export type SpacefillProps = Partial<typeof DefaultSpacefillProps>

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

                if (i % 10000 === 0 && ctx.shouldUpdate) {
                    await ctx.update({ message: 'Spacefill', current: i, max: elementCount });
                }
            }

            const mesh = meshBuilder.getMesh()
            // console.log(mesh)
            if (!mesh.offsetBuffer.ref.value) return

            const unitCount = units.length
            const transformArray = new Float32Array(unitCount * 16)
            for (let i = 0; i < unitCount; i++) {
                Mat4.toArray(units[i].operator.matrix, transformArray, i * 16)
            }

            // console.log({ unitCount, elementCount })

            // const color = createUniformColor({ value: 0xFF4411 })

            // const color = elementSymbolColorData({ units, elementGroup, mesh })

            const scale = ColorScale.create({ domain: [ 0, unitCount - 1 ] })
            const color = createInstanceColor({
                colorFn: scale.color,
                unitCount
            })

            // const scale = ColorScale.create({ domain: [ 0, unitCount * elementCount - 1 ] })
            // const color = createElementInstanceColor({
            //     colorFn: (unitIdx, elementIdx) => scale.color(unitIdx * elementCount + elementIdx),
            //     unitCount,
            //     offsetCount: mesh.offsetCount,
            //     offsets: mesh.offsetBuffer as any
            // })

            const spheres = createMeshRenderObject({
                objectId: 0,

                position: mesh.vertexBuffer,
                normal: mesh.normalBuffer as ValueCell<Float32Array>,
                color: color as ColorData,
                id: mesh.idBuffer as ValueCell<Float32Array>,
                transform: ValueCell.create(transformArray),
                index: mesh.indexBuffer,

                instanceCount: unitCount,
                indexCount: mesh.triangleCount,
                elementCount: mesh.offsetCount - 1,
                positionCount: mesh.vertexCount
            })
            renderObjects.push(spheres)
        }),
        update: (props: RepresentationProps) => false
    }
}
