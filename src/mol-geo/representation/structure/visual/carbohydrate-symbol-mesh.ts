/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from 'mol-util/value-cell'

import { createMeshRenderObject, MeshRenderObject } from 'mol-gl/render-object'
import { Unit, Structure } from 'mol-model/structure';
import { DefaultStructureProps, StructureVisual } from '..';
import { RuntimeContext } from 'mol-task'
import { createIdentityTransform } from './util/common';
import { MeshValues } from 'mol-gl/renderable';
import { getMeshData } from '../../../util/mesh-data';
import { Mesh } from '../../../shape/mesh';
import { PickingId } from '../../../util/picking';
import { createMarkers, MarkerAction } from '../../../util/marker-data';
import { Loci, EmptyLoci } from 'mol-model/loci';
import { SizeTheme } from '../../../theme';
import { createMeshValues, updateMeshValues, updateRenderableState, createRenderableState, DefaultMeshProps } from '../../util';
import { MeshBuilder } from '../../../shape/mesh-builder';
import { Vec3, Mat4 } from 'mol-math/linear-algebra';
import { createUniformColor } from '../../../util/color-data';

async function createCarbohydrateSymbolMesh(ctx: RuntimeContext, structure: Structure, mesh?: Mesh) {
    const builder = MeshBuilder.create(256, 128, mesh)

    const t = Mat4.identity()
    const p = Vec3.zero()
    const { carbohydrates } = structure

    const linkParams = { radiusTop: 0.2, radiusBottom: 0.2 }

    for (let i = 0, il = carbohydrates.elements.length; i < il; ++i) {
        const c = carbohydrates.elements[i]
        Mat4.setTranslation(t, c.center)
        builder.addBox(t, { width: 2, height: 2, depth: 2 })
    }

    for (let i = 0, il = carbohydrates.links.length; i < il; ++i) {
        const l = carbohydrates.links[i]
        const centerA = carbohydrates.elements[l.carbohydrateIndexA].center
        const centerB = carbohydrates.elements[l.carbohydrateIndexB].center
        builder.addCylinder(centerA, centerB, 0.5, linkParams)
    }

    for (let i = 0, il = carbohydrates.terminalLinks.length; i < il; ++i) {
        const tl = carbohydrates.terminalLinks[i]
        const center = carbohydrates.elements[tl.carbohydrateIndex].center
        tl.elementUnit.conformation.position(tl.elementUnit.elements[tl.elementIndex], p)
        if (tl.fromCarbohydrate) {
            builder.addCylinder(center, p, 0.5, linkParams)
        } else {
            builder.addCylinder(p, center, 0.5, linkParams)
        }
    }

    return builder.getMesh()
}

export const DefaultCarbohydrateSymbolProps = {
    ...DefaultMeshProps,
    ...DefaultStructureProps,
    sizeTheme: { name: 'physical', factor: 1 } as SizeTheme,
    detail: 0,
    unitKinds: [ Unit.Kind.Atomic, Unit.Kind.Spheres ] as Unit.Kind[]
}
export type CarbohydrateSymbolProps = Partial<typeof DefaultCarbohydrateSymbolProps>

export function CarbohydrateSymbolVisual(): StructureVisual<CarbohydrateSymbolProps> {
    let renderObject: MeshRenderObject
    let currentProps: typeof DefaultCarbohydrateSymbolProps
    let mesh: Mesh
    let currentStructure: Structure

    return {
        get renderObject () { return renderObject },
        async create(ctx: RuntimeContext, structure: Structure, props: CarbohydrateSymbolProps = {}) {
            currentProps = Object.assign({}, DefaultCarbohydrateSymbolProps, props)
            currentStructure = structure

            const instanceCount = 1
            const elementCount = currentStructure.elementCount

            mesh = await createCarbohydrateSymbolMesh(ctx, currentStructure, mesh)
            // console.log(mesh)

            const transforms = createIdentityTransform()
            const color = createUniformColor({ value: 0x999911 }) // TODO
            const marker = createMarkers(instanceCount * elementCount)

            const counts = { drawCount: mesh.triangleCount * 3, elementCount, instanceCount }

            const values: MeshValues = {
                ...getMeshData(mesh),
                ...color,
                ...marker,
                aTransform: transforms,
                elements: mesh.indexBuffer,
                ...createMeshValues(currentProps, counts),
                aColor: ValueCell.create(new Float32Array(mesh.vertexCount * 3))
            }
            const state = createRenderableState(currentProps)

            renderObject = createMeshRenderObject(values, state)
        },
        async update(ctx: RuntimeContext, props: CarbohydrateSymbolProps) {
            const newProps = Object.assign({}, currentProps, props)

            if (!renderObject) return false

            updateMeshValues(renderObject.values, newProps)
            updateRenderableState(renderObject.state, newProps)

            return false
        },
        getLoci(pickingId: PickingId) {
            return EmptyLoci
        },
        mark(loci: Loci, action: MarkerAction) {
            // TODO
        },
        destroy() {
            // TODO
        }
    }
}
