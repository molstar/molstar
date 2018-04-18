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
import { ElementColor } from '../../color';

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

                meshBuilder.addIcosahedron(m, {
                    radius: VdwRadius(type_symbol.value(e)),
                    color: ElementColor(type_symbol.value(e)),
                    detail
                })

                if (i % 10 === 0 && ctx.shouldUpdate) {
                    await ctx.update({ message: 'Spacefill', current: i, max: elementCount });
                }
            }

            const unitsCount = units.length
            const transformArray = new Float32Array(unitsCount * 16)
            for (let i = 0; i < unitsCount; i++) {
                Mat4.toArray(units[i].operator.matrix, transformArray, i * 16)
            }

            // const color = ValueCell.create(createColorTexture(unitsCount))
            // color.ref.value.set([ 0, 0, 255 ])

            const mesh = meshBuilder.getMesh()
            console.log(mesh)

            const spheres = createRenderObject('mesh', {
                position: mesh.vertexBuffer,
                normal: mesh.normalBuffer,
                color: { type: 'attribute', value: (mesh as any).colorBuffer },
                transform: ValueCell.create(transformArray),
                elements: mesh.indexBuffer,

                instanceCount: unitsCount,
                elementCount: mesh.triangleCount,
                positionCount: mesh.vertexCount
            }, {})
            renderObjects.push(spheres)
        }),
        update: (props: RepresentationProps) => false
    }
}
