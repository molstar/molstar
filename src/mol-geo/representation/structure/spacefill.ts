/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from 'mol-util/value-cell'

import { RenderObject, createMeshRenderObject, MeshRenderObject } from 'mol-gl/scene'
// import { createColorTexture } from 'mol-gl/util';
import { Vec3, Mat4 } from 'mol-math/linear-algebra'
import { OrderedSet } from 'mol-data/int'
import { Unit, ElementGroup } from 'mol-model/structure';
import { RepresentationProps, UnitsRepresentation } from './index';
import { Task } from 'mol-task'
import { MeshBuilder } from '../../shape/mesh-builder';
import { VdwRadius } from 'mol-model/structure/model/properties/atomic';
import { createTransforms, createColors } from './utils';
import { ColorTheme } from '../../theme';
import VertexMap from '../../shape/vertex-map';

export const DefaultSpacefillProps = {
    detail: 0,
    colorTheme: { name: 'instance-index' } as ColorTheme,
}
export type SpacefillProps = Partial<typeof DefaultSpacefillProps>

function createSpacefillMesh(unit: Unit, elementGroup: ElementGroup, detail: number) {
    return Task.create('Spacefill', async ctx => {
        const meshBuilder = MeshBuilder.create()

        const v = Vec3.zero()
        const m = Mat4.identity()

        const { x, y, z } = unit.model.atomSiteConformation
        const { type_symbol } = unit.model.hierarchy.atoms
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

        return meshBuilder.getMesh()
    })
}

export default function Spacefill(): UnitsRepresentation<SpacefillProps> {
    const renderObjects: RenderObject[] = []
    let spheres: MeshRenderObject

    return {
        renderObjects,
        create(units: ReadonlyArray<Unit>, elementGroup: ElementGroup, props: SpacefillProps = {}) {
            return Task.create('Spacefill.create', async ctx => {
                renderObjects.length = 0 // clear

                const { detail, colorTheme } = { ...DefaultSpacefillProps, ...props }

                await ctx.update('Computing spacefill mesh');
                const mesh = await ctx.runChild(createSpacefillMesh(units[0], elementGroup, detail))
                // console.log(mesh)

                const vertexMap = VertexMap.fromMesh(mesh)

                await ctx.update('Computing spacefill transforms');
                const transforms = createTransforms(units)

                await ctx.update('Computing spacefill colors');
                const color = createColors(units, elementGroup, vertexMap, colorTheme)

                spheres = createMeshRenderObject({
                    objectId: 0,

                    position: mesh.vertexBuffer,
                    normal: mesh.normalBuffer as ValueCell<Float32Array>,
                    color: color,
                    id: mesh.idBuffer as ValueCell<Float32Array>,
                    transform: ValueCell.create(transforms),
                    index: mesh.indexBuffer,

                    instanceCount: units.length,
                    indexCount: mesh.triangleCount,
                    elementCount: OrderedSet.size(elementGroup.elements),
                    positionCount: mesh.vertexCount
                })
                renderObjects.push(spheres)
            })
        },
        update(props: RepresentationProps) {
            return Task.create('Spacefill.update', async ctx => {
                if (!spheres) return false

                return false
            })
        }
    }
}
