/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ValueCell } from 'mol-util/value-cell'

import { RenderObject, createMeshRenderObject, MeshRenderObject } from 'mol-gl/scene'
// import { createColorTexture } from 'mol-gl/util';
import { Vec3, Mat4 } from 'mol-math/linear-algebra'
import { Unit, Element, Queries, StructureSymmetry } from 'mol-model/structure';
import { UnitsRepresentation } from './index';
import { Task } from 'mol-task'
import { MeshBuilder } from '../../shape/mesh-builder';
import { createTransforms, createColors } from './utils';
import { ColorTheme } from '../../theme';
import VertexMap from '../../shape/vertex-map';
import { icosahedronVertexCount } from '../../primitive/icosahedron';

export const DefaultSpacefillProps = {
    detail: 0,
    colorTheme: { name: 'instance-index' } as ColorTheme,
    alpha: 1,
    visible: true,
    doubleSided: false
}
export type SpacefillProps = Partial<typeof DefaultSpacefillProps>

function createSpacefillMesh(unit: Unit, detail: number) {
    return Task.create('Sphere mesh', async ctx => {
        const { elements } = unit;
        const elementCount = elements.length;
        const vertexCount = elementCount * icosahedronVertexCount(detail)
        const meshBuilder = MeshBuilder.create(vertexCount)

        let radius: Element.Property<number>
        if (Unit.isAtomic(unit)) {
            radius = Queries.props.atom.vdw_radius
        } else if (Unit.isSpheres(unit)) {
            radius = Queries.props.coarse_grained.sphere_radius
        } else {
            console.warn('Unsupported unit type')
            return meshBuilder.getMesh()
        }

        const v = Vec3.zero()
        const m = Mat4.identity()

        const { x, y, z } = unit.conformation
        const l = Element.Location()
        l.unit = unit

        for (let i = 0; i < elementCount; i++) {
            l.element = elements[i]
            v[0] = x(l.element)
            v[1] = y(l.element)
            v[2] = z(l.element)
            Mat4.setTranslation(m, v)

            meshBuilder.setId(i)
            meshBuilder.addIcosahedron(m, { radius: radius(l), detail })

            if (i % 10000 === 0 && ctx.shouldUpdate) {
                await ctx.update({ message: 'Sphere mesh', current: i, max: elementCount });
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
        create(group: StructureSymmetry.UnitGroup, props: SpacefillProps = {}) {
            return Task.create('Spacefill.create', async ctx => {
                renderObjects.length = 0 // clear

                const { detail, colorTheme, alpha, visible, doubleSided } = { ...DefaultSpacefillProps, ...props }

                await ctx.update('Computing spacefill mesh');
                const mesh = await ctx.runChild(createSpacefillMesh(group.units[0], detail))
                // console.log(mesh)

                const vertexMap = VertexMap.fromMesh(mesh)

                await ctx.update('Computing spacefill transforms');
                const transforms = createTransforms(group)

                await ctx.update('Computing spacefill colors');
                const color = createColors(group, vertexMap, colorTheme)

                spheres = createMeshRenderObject({
                    objectId: 0,
                    alpha,
                    visible,
                    doubleSided,

                    position: mesh.vertexBuffer,
                    normal: mesh.normalBuffer as ValueCell<Float32Array>,
                    color: color,
                    id: mesh.idBuffer as ValueCell<Float32Array>,
                    transform: ValueCell.create(transforms),
                    index: mesh.indexBuffer,

                    instanceCount: group.units.length,
                    indexCount: mesh.triangleCount,
                    elementCount: group.elements.length,
                    positionCount: mesh.vertexCount
                })
                renderObjects.push(spheres)
            })
        },
        update(props: SpacefillProps) {
            return Task.create('Spacefill.update', async ctx => {
                if (!spheres) return false

                return false
            })
        }
    }
}
