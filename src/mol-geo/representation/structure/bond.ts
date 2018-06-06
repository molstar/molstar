/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ValueCell } from 'mol-util/value-cell'

import { RenderObject, createMeshRenderObject, MeshRenderObject } from 'mol-gl/render-object'
import { Unit, Element } from 'mol-model/structure';
import { UnitsRepresentation, DefaultStructureProps } from './index';
import { Task } from 'mol-task'
import { createTransforms, createEmptyFlags } from './utils';
import { fillSerial } from 'mol-gl/renderable/util';
import { RenderableState, MeshValues } from 'mol-gl/renderable';
import { getMeshData } from '../../util/mesh-data';
import { Mesh } from '../../shape/mesh';
import { PickingId } from '../../util/picking';
import { MeshBuilder } from '../../shape/mesh-builder';
import { Vec3, Mat4 } from 'mol-math/linear-algebra';
import { createUniformColor } from '../../util/color-data';
import { defaults } from 'mol-util';
import { SortedArray } from 'mol-data/int';

function createBondMesh(unit: Unit, mesh?: Mesh) {
    return Task.create('Cylinder mesh', async ctx => {
        if (!Unit.isAtomic(unit)) return Mesh.createEmpty(mesh)

        const elements = unit.elements;
        const { count, offset, neighbor } = unit.bonds;
    
        if (!count) return Mesh.createEmpty(mesh)

        // TODO calculate vertextCount properly
        const vertexCount = 32 * count
        const meshBuilder = MeshBuilder.create(vertexCount, vertexCount / 2, mesh)

        const va = Vec3.zero()
        const vb = Vec3.zero()
        const vt = Vec3.zero()
        const m = Mat4.identity()

        const { x, y, z } = unit.conformation
        const l = Element.Location()
        l.unit = unit

        for (let j = 0; j < offset.length - 1; ++j) {
            const start = offset[j]
            const end = offset[j + 1]

            if (end <= start) continue

            const aI = elements[j]
            va[0] = x(aI)
            va[1] = y(aI)
            va[2] = z(aI)
            for (let _bI = start; _bI < end; ++_bI) {
                const bI = elements[neighbor[_bI]]
                if (bI > aI) continue

                vb[0] = x(bI)
                vb[1] = y(bI)
                vb[2] = z(bI)

                Vec3.scale(vt, Vec3.add(vt, va, vb), 0.5)
                Vec3.makeRotation(m, Vec3.create(0, 1, 0), Vec3.sub(vb, vb, va))
                Mat4.setTranslation(m, vt)

                meshBuilder.setId(j)
                meshBuilder.addCylinder(m, { radiusTop: 0.2, radiusBottom: 0.2 })
            }

            if (j % 10000 === 0 && ctx.shouldUpdate) {
                await ctx.update({ message: 'Cylinder mesh', current: j, max: count });
            }
        }

        return meshBuilder.getMesh()
    })
}

export const DefaultBondProps = {
    ...DefaultStructureProps,
    flipSided: false,
    flatShaded: false,
}
export type BondProps = Partial<typeof DefaultBondProps>

export default function Bond(): UnitsRepresentation<BondProps> {
    const renderObjects: RenderObject[] = []
    let cylinders: MeshRenderObject
    let currentProps: typeof DefaultBondProps
    let mesh: Mesh
    let currentGroup: Unit.SymmetryGroup
    // let vertexMap: VertexMap

    return {
        renderObjects,
        create(group: Unit.SymmetryGroup, props: BondProps = {}) {
            currentProps = Object.assign({}, DefaultBondProps, props)

            return Task.create('Bond.create', async ctx => {
                renderObjects.length = 0 // clear
                currentGroup = group

                mesh = await createBondMesh(group.units[0]).runAsChild(ctx, 'Computing bond mesh')
                // console.log(mesh)
                // vertexMap = VertexMap.fromMesh(mesh)

                await ctx.update('Computing bond transforms');
                const transforms = createTransforms(group)

                await ctx.update('Computing bond colors');
                const color = createUniformColor({ value: 0xFF0000 })

                await ctx.update('Computing bond flags');
                const flag = createEmptyFlags()

                const instanceCount = group.units.length

                const values: MeshValues = {
                    ...getMeshData(mesh),
                    aTransform: transforms,
                    aInstanceId: ValueCell.create(fillSerial(new Float32Array(instanceCount))),
                    ...color,
                    ...flag,

                    uAlpha: ValueCell.create(defaults(props.alpha, 1.0)),
                    uInstanceCount: ValueCell.create(instanceCount),
                    uElementCount: ValueCell.create(group.elements.length),

                    elements: mesh.indexBuffer,

                    drawCount: ValueCell.create(mesh.triangleCount * 3),
                    instanceCount: ValueCell.create(instanceCount),

                    dDoubleSided: ValueCell.create(defaults(props.doubleSided, true)),
                    dFlatShaded: ValueCell.create(defaults(props.flatShaded, false)),
                    dFlipSided: ValueCell.create(defaults(props.flipSided, false)),
                }
                const state: RenderableState = {
                    depthMask: defaults(props.depthMask, true),
                    visible: defaults(props.visible, true)
                }

                cylinders = createMeshRenderObject(values, state)
                renderObjects.push(cylinders)
            })
        },
        update(props: BondProps) {
            const newProps = Object.assign({}, currentProps, props)

            return Task.create('Bond.update', async ctx => {
                if (!cylinders) return false
                // TODO

                ValueCell.updateIfChanged(cylinders.values.uAlpha, newProps.alpha)
                ValueCell.updateIfChanged(cylinders.values.dDoubleSided, newProps.doubleSided)
                ValueCell.updateIfChanged(cylinders.values.dFlipSided, newProps.flipSided)
                ValueCell.updateIfChanged(cylinders.values.dFlatShaded, newProps.flatShaded)

                cylinders.state.visible = newProps.visible
                cylinders.state.depthMask = newProps.depthMask

                return true
            })
        },
        getLoci(pickingId: PickingId) {
            // const { objectId, instanceId, elementId } = pickingId
            // if (cylinders.id === objectId) {
            //     const unit = currentGroup.units[instanceId]
            //     const elements = SortedArray.ofSingleton(currentGroup.elements[elementId])
            //     return Element.Loci([{ unit, elements }])
            // }
            return null
        }
    }
}
