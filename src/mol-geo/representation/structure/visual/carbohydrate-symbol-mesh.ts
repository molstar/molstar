/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from 'mol-util/value-cell'

import { createMeshRenderObject, MeshRenderObject } from 'mol-gl/render-object'
import { Unit, Structure, StructureElement } from 'mol-model/structure';
import { DefaultStructureProps, StructureVisual } from '..';
import { RuntimeContext } from 'mol-task'
import { createIdentityTransform, createColors } from './util/common';
import { MeshValues } from 'mol-gl/renderable';
import { getMeshData } from '../../../util/mesh-data';
import { Mesh } from '../../../shape/mesh';
import { PickingId } from '../../../util/picking';
import { createMarkers, MarkerAction, MarkerData, applyMarkerAction } from '../../../util/marker-data';
import { Loci, EmptyLoci, isEveryLoci } from 'mol-model/loci';
import { SizeTheme } from '../../../theme';
import { createMeshValues, updateMeshValues, updateRenderableState, createRenderableState, DefaultMeshProps } from '../../util';
import { MeshBuilder } from '../../../shape/mesh-builder';
import { Vec3, Mat4 } from 'mol-math/linear-algebra';
import { getSaccharideShape, SaccharideShapes } from 'mol-model/structure/structure/carbohydrates/constants';
import { deepEqual } from 'mol-util';
import { LocationIterator } from './util/location-iterator';
import { OrderedSet } from 'mol-data/int';

const t = Mat4.identity()
const sVec = Vec3.zero()
const pd = Vec3.zero()

async function createCarbohydrateSymbolMesh(ctx: RuntimeContext, structure: Structure, mesh?: Mesh) {
    const builder = MeshBuilder.create(256, 128, mesh)

    const carbohydrates = structure.carbohydrates

    const side = 1.75 * 2 * 0.806; // 0.806 == Math.cos(Math.PI / 4)
    const radius = 1.75

    for (let i = 0, il = carbohydrates.elements.length; i < il; ++i) {
        const c = carbohydrates.elements[i];
        const shapeType = getSaccharideShape(c.component.type)

        const { center, normal, direction } = c.geometry
        Vec3.add(pd, center, direction)
        Mat4.targetTo(t, center, pd, normal)
        Mat4.setTranslation(t, center)

        builder.setId(i)

        switch (shapeType) {
            case SaccharideShapes.FilledSphere:
                builder.addSphere(center, radius, 2)
                break;
            case SaccharideShapes.FilledCube:
                Mat4.scaleUniformly(t, t, side)
                builder.addBox(t)
                break;
            case SaccharideShapes.CrossedCube:
                // TODO split
                Mat4.scaleUniformly(t, t, side)
                builder.addBox(t)
                break;
            case SaccharideShapes.FilledCone:
                Mat4.scaleUniformly(t, t, side * 1.2)
                builder.addOctagonalPyramid(t)
                break
            case SaccharideShapes.DevidedCone:
                // TODO split
                Mat4.scaleUniformly(t, t, side * 1.2)
                builder.addOctagonalPyramid(t)
                break
            case SaccharideShapes.FlatBox:
                Mat4.mul(t, t, Mat4.rotZY90)
                Mat4.scale(t, t, Vec3.set(sVec, side, side, side / 2))
                builder.addBox(t)
                break
            case SaccharideShapes.FilledStar:
                Mat4.mul(t, t, Mat4.rotZY90)
                builder.addStar(t, { outerRadius: side, innerRadius: side / 2, thickness: side / 2, pointCount: 5 })
                break
            case SaccharideShapes.FilledDiamond:
                Mat4.mul(t, t, Mat4.rotZY90)
                Mat4.scale(t, t, Vec3.set(sVec, side * 1.4, side * 1.4, side * 1.4))
                builder.addOctahedron(t)
                break
            case SaccharideShapes.DividedDiamond:
                // TODO split
                Mat4.mul(t, t, Mat4.rotZY90)
                Mat4.scale(t, t, Vec3.set(sVec, side * 1.4, side * 1.4, side * 1.4))
                builder.addOctahedron(t)
                break
            case SaccharideShapes.FlatDiamond:
                Mat4.mul(t, t, Mat4.rotZY90)
                Mat4.scale(t, t, Vec3.set(sVec, side, side / 2, side / 2))
                builder.addDiamondPrism(t)
                break
            case SaccharideShapes.Pentagon:
                Mat4.mul(t, t, Mat4.rotZY90)
                Mat4.scale(t, t, Vec3.set(sVec, side, side, side / 2))
                builder.addPentagonalPrism(t)
                break
            case SaccharideShapes.FlatHexagon:
            default:
                Mat4.mul(t, t, Mat4.rotZYZ90)
                Mat4.scale(t, t, Vec3.set(sVec, side / 1.5, side , side / 2))
                builder.addHexagonalPrism(t)
                break
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

            const { colorTheme } = { ...DefaultCarbohydrateSymbolProps, ...props }
            const instanceCount = 1
            const elementCount = currentStructure.elementCount

            mesh = await createCarbohydrateSymbolMesh(ctx, currentStructure, mesh)

            const transforms = createIdentityTransform()
            const color = createColors(createCarbohydrateElementIterator(structure), colorTheme)
            const marker = createMarkers(instanceCount * elementCount)

            const counts = { drawCount: mesh.triangleCount * 3, elementCount, instanceCount }

            const values: MeshValues = {
                ...getMeshData(mesh),
                ...color,
                ...marker,
                aTransform: transforms,
                elements: mesh.indexBuffer,
                ...createMeshValues(currentProps, counts)
            }
            const state = createRenderableState(currentProps)

            renderObject = createMeshRenderObject(values, state)
        },
        async update(ctx: RuntimeContext, props: CarbohydrateSymbolProps) {
            const newProps = Object.assign({}, currentProps, props)

            if (!renderObject) return false

            let updateColor = false

            if (!deepEqual(newProps.colorTheme, currentProps.colorTheme)) {
                updateColor = true
            }

            if (updateColor) {
                createColors(createCarbohydrateElementIterator(currentStructure), newProps.colorTheme, renderObject.values)
            }

            updateMeshValues(renderObject.values, newProps)
            updateRenderableState(renderObject.state, newProps)

            currentProps = newProps
            return false
        },
        getLoci(pickingId: PickingId) {
            return getCarbohydrateLoci(pickingId, currentStructure, renderObject.id)
        },
        mark(loci: Loci, action: MarkerAction) {
            markCarbohydrate(loci, action, currentStructure, renderObject.values)
        },
        destroy() {
            // TODO
        }
    }
}

function createCarbohydrateElementIterator(structure: Structure): LocationIterator {
    const carbElements = structure.carbohydrates.elements
    const elementCount = carbElements.length
    const instanceCount = 1
    const location = StructureElement.create()
    const getLocation = (elementIndex: number, instanceIndex: number) => {
        const carb = carbElements[elementIndex]
        location.unit = carb.unit
        location.element = carb.anomericCarbon
        return location
    }
    return LocationIterator(elementCount, instanceCount, getLocation)
}

function getCarbohydrateLoci(pickingId: PickingId, structure: Structure, id: number) {
    const { objectId, elementId } = pickingId
    if (id === objectId) {
        const carb = structure.carbohydrates.elements[elementId]
        const { unit } = carb
        const index = OrderedSet.findPredecessorIndex(unit.elements, carb.anomericCarbon)
        const indices = OrderedSet.ofSingleton(index as StructureElement.UnitIndex)
        return StructureElement.Loci([{ unit, indices }])
    }
    return EmptyLoci
}

function markCarbohydrate(loci: Loci, action: MarkerAction, structure: Structure, values: MarkerData) {
    const tMarker = values.tMarker

    const { getElementIndex } = structure.carbohydrates
    const elementCount = structure.carbohydrates.elements.length

    let changed = false
    const array = tMarker.ref.value.array
    if (isEveryLoci(loci)) {
        if (applyMarkerAction(array, 0, elementCount, action)) {
            changed = true
        }
    } else if (StructureElement.isLoci(loci)) {
        for (const e of loci.elements) {
            OrderedSet.forEach(e.indices, index => {
                const idx = getElementIndex(e.unit, e.unit.elements[index])
                if (idx !== undefined) {
                    if (applyMarkerAction(array, idx, idx + 1, action) && !changed) {
                        changed = true
                    }
                }
            })
        }
    } else {
        return
    }
    if (changed) {
        ValueCell.update(tMarker, tMarker.ref.value)
    }
}