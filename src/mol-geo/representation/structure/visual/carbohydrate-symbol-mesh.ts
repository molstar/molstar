/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Unit, Structure, StructureElement } from 'mol-model/structure';
import { ComplexVisual } from '..';
import { RuntimeContext } from 'mol-task'
import { Mesh } from '../../../shape/mesh';
import { PickingId } from '../../../util/picking';
import { Loci, EmptyLoci } from 'mol-model/loci';
import { MeshBuilder } from '../../../shape/mesh-builder';
import { Vec3, Mat4 } from 'mol-math/linear-algebra';
import { getSaccharideShape, SaccharideShapes } from 'mol-model/structure/structure/carbohydrates/constants';
import { LocationIterator } from './util/location-iterator';
import { OrderedSet, Interval } from 'mol-data/int';
import { ComplexMeshVisual, DefaultComplexMeshProps } from '../complex-visual';
import { SizeThemeProps } from 'mol-view/theme/size';

const t = Mat4.identity()
const sVec = Vec3.zero()
const pd = Vec3.zero()

async function createCarbohydrateSymbolMesh(ctx: RuntimeContext, structure: Structure, props: CarbohydrateSymbolProps, mesh?: Mesh) {
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

        builder.setId(i * 2)

        switch (shapeType) {
            case SaccharideShapes.FilledSphere:
                builder.addSphere(center, radius, 2)
                break;
            case SaccharideShapes.FilledCube:
                Mat4.scaleUniformly(t, t, side)
                builder.addBox(t)
                break;
            case SaccharideShapes.CrossedCube:
                Mat4.scaleUniformly(t, t, side)
                builder.addPerforatedBox(t)
                Mat4.mul(t, t, Mat4.rotZ90X180)
                builder.setId(i * 2 + 1)
                builder.addPerforatedBox(t)
                break;
            case SaccharideShapes.FilledCone:
                Mat4.scaleUniformly(t, t, side * 1.2)
                builder.addOctagonalPyramid(t)
                break
            case SaccharideShapes.DevidedCone:
                Mat4.scaleUniformly(t, t, side * 1.2)
                builder.addPerforatedOctagonalPyramid(t)
                Mat4.mul(t, t, Mat4.rotZ90)
                builder.setId(i * 2 + 1)
                builder.addPerforatedOctagonalPyramid(t)
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
                Mat4.mul(t, t, Mat4.rotZY90)
                Mat4.scale(t, t, Vec3.set(sVec, side * 1.4, side * 1.4, side * 1.4))
                builder.addPerforatedOctahedron(t)
                Mat4.mul(t, t, Mat4.rotY90)
                builder.setId(i * 2 + 1)
                builder.addPerforatedOctahedron(t)
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
    ...DefaultComplexMeshProps,
    sizeTheme: { name: 'physical', factor: 1 } as SizeThemeProps,
    detail: 0,
    unitKinds: [ Unit.Kind.Atomic, Unit.Kind.Spheres ] as Unit.Kind[]
}
export type CarbohydrateSymbolProps = typeof DefaultCarbohydrateSymbolProps

export function CarbohydrateSymbolVisual(): ComplexVisual<CarbohydrateSymbolProps> {
    return ComplexMeshVisual<CarbohydrateSymbolProps>({
        defaultProps: DefaultCarbohydrateSymbolProps,
        createMesh: createCarbohydrateSymbolMesh,
        createLocationIterator: CarbohydrateElementIterator,
        getLoci: getCarbohydrateLoci,
        mark: markCarbohydrate,
        setUpdateState: () => {}
    })
}

function CarbohydrateElementIterator(structure: Structure): LocationIterator {
    const carbElements = structure.carbohydrates.elements
    const elementCount = carbElements.length * 2
    const instanceCount = 1
    const location = StructureElement.create()
    function getLocation (elementIndex: number, instanceIndex: number) {
        const carb = carbElements[Math.floor(elementIndex / 2)]
        location.unit = carb.unit
        location.element = carb.anomericCarbon
        return location
    }
    function isSecondary (elementIndex: number, instanceIndex: number) {
        return (elementIndex % 2) === 1
    }
    return LocationIterator(elementCount, instanceCount, getLocation, isSecondary)
}

function getCarbohydrateLoci(pickingId: PickingId, structure: Structure, id: number) {
    const { objectId, elementId } = pickingId
    if (id === objectId) {
        const carb = structure.carbohydrates.elements[Math.floor(elementId / 2)]
        const { unit } = carb
        const index = OrderedSet.findPredecessorIndex(unit.elements, carb.anomericCarbon)
        const indices = OrderedSet.ofSingleton(index as StructureElement.UnitIndex)
        return StructureElement.Loci([{ unit, indices }])
    }
    return EmptyLoci
}

function markCarbohydrate(loci: Loci, structure: Structure, apply: (interval: Interval) => boolean) {
    const { getElementIndex } = structure.carbohydrates

    let changed = false
    if (StructureElement.isLoci(loci)) {
        for (const e of loci.elements) {
            OrderedSet.forEach(e.indices, index => {
                const idx = getElementIndex(e.unit, e.unit.elements[index])
                if (idx !== undefined) {
                    if (apply(Interval.ofBounds(idx * 2, idx * 2 + 2))) changed = true
                }
            })
        }
    }
    return changed
}