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
import { getSaccharideShape, SaccharideShapes } from 'mol-model/structure/structure/carbohydrates/constants';

async function createCarbohydrateSymbolMesh(ctx: RuntimeContext, structure: Structure, mesh?: Mesh) {
    const builder = MeshBuilder.create(256, 128, mesh)

    const t = Mat4.identity()
    const sMat = Mat4.identity()
    const sVec = Vec3.zero()
    const p = Vec3.zero()
    const pd = Vec3.zero()
    const p1 = Vec3.zero()
    const p2 = Vec3.zero()
    const carbohydrates = structure.carbohydrates

    function centerAlign(center: Vec3, normal: Vec3, direction: Vec3) {
        Vec3.add(pd, center, direction)
        Mat4.targetTo(t, center, pd, normal)
        Mat4.setTranslation(t, center)
    }

    const side = 1.75 * 2 * 0.806; // 0.806 == Math.cos(Math.PI / 4)
    const radius = 1.75
    const coneParams = { radiusTop: 0.0, radiusBottom: radius, bottomCap: true }

    const linkParams = { radiusTop: 0.4, radiusBottom: 0.4 }

    for (let i = 0, il = carbohydrates.elements.length; i < il; ++i) {
        const c = carbohydrates.elements[i];
        if (!c.hasRing) continue;

        const cGeo = c.geometry!
        const shapeType = getSaccharideShape(c.component.type)
        switch (shapeType) {
            case SaccharideShapes.FilledSphere:
                builder.addSphere(cGeo.center, radius, 2)
                break;
            case SaccharideShapes.FilledCube:
                centerAlign(cGeo.center, cGeo.normal, cGeo.direction)
                builder.addBox(t, { width: side, height: side, depth: side })
                builder.addOctahedron(t)
                break;
            case SaccharideShapes.CrossedCube:
                // TODO split
                centerAlign(cGeo.center, cGeo.normal, cGeo.direction)
                builder.addBox(t, { width: side, height: side, depth: side })
                break;
            case SaccharideShapes.FilledCone:
                Vec3.scaleAndAdd(p1, cGeo.center, cGeo.direction, radius)
                Vec3.scaleAndSub(p2, cGeo.center, cGeo.direction, radius)
                builder.addCylinder(p1, p2, 1, coneParams)
                break
            case SaccharideShapes.DevidedCone:
                // TODO split
                Vec3.scaleAndAdd(p1, cGeo.center, cGeo.direction, radius)
                Vec3.scaleAndSub(p2, cGeo.center, cGeo.direction, radius)
                builder.addCylinder(p1, p2, 1, coneParams)
                break
            case SaccharideShapes.FlatBox:
                centerAlign(cGeo.center, cGeo.normal, cGeo.direction)
                builder.addBox(t, { width: side, height: side / 2, depth: side })
                break
            case SaccharideShapes.FilledStar:
                centerAlign(cGeo.center, cGeo.normal, cGeo.direction)
                builder.addStar(t, { outerRadius: side, innerRadius: side / 2, thickness: side / 2, pointCount: 5 })
                break
            case SaccharideShapes.FilledDiamond:
                centerAlign(cGeo.center, cGeo.normal, cGeo.direction)
                Mat4.fromScaling(sMat, Vec3.set(sVec, side * 1.4, side * 1.4, side * 1.4))
                Mat4.mul(t, t, sMat)
                builder.addOctahedron(t)
                break
            case SaccharideShapes.DividedDiamond:
                // TODO split
                centerAlign(cGeo.center, cGeo.normal, cGeo.direction)
                Mat4.fromScaling(sMat, Vec3.set(sVec, side * 1.4, side * 1.4, side * 1.4))
                Mat4.mul(t, t, sMat)
                builder.addOctahedron(t)
                break
            case SaccharideShapes.FlatDiamond:
            case SaccharideShapes.Pentagon:
                centerAlign(cGeo.center, cGeo.normal, cGeo.direction)
                builder.addBox(t, { width: side, height: side, depth: 0.5 })
                break
            case SaccharideShapes.FlatHexagon:
            default:
                centerAlign(cGeo.center, cGeo.normal, cGeo.direction)
                builder.addBox(t, { width: side, height: side, depth: 0.1 })
                break
        }
    }

    for (let i = 0, il = carbohydrates.links.length; i < il; ++i) {
        const l = carbohydrates.links[i]
        const centerA = carbohydrates.elements[l.carbohydrateIndexA].geometry!.center
        const centerB = carbohydrates.elements[l.carbohydrateIndexB].geometry!.center
        builder.addCylinder(centerA, centerB, 0.5, linkParams)
    }

    for (let i = 0, il = carbohydrates.terminalLinks.length; i < il; ++i) {
        const tl = carbohydrates.terminalLinks[i]
        const center = carbohydrates.elements[tl.carbohydrateIndex].geometry!.center
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
