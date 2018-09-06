/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Unit, Structure, StructureElement } from 'mol-model/structure';
import { ComplexVisual, MeshUpdateState } from '..';
import { RuntimeContext } from 'mol-task'
import { Mesh } from '../../../mesh/mesh';
import { PickingId } from '../../../util/picking';
import { Loci, EmptyLoci } from 'mol-model/loci';
import { MeshBuilder } from '../../../mesh/mesh-builder';
import { Vec3, Mat4 } from 'mol-math/linear-algebra';
import { getSaccharideShape, SaccharideShapes } from 'mol-model/structure/structure/carbohydrates/constants';
import { LocationIterator } from '../../../util/location-iterator';
import { OrderedSet, Interval } from 'mol-data/int';
import { ComplexMeshVisual, DefaultComplexMeshProps } from '../complex-visual';
import { SizeThemeProps, SizeTheme } from 'mol-view/theme/size';
import { addSphere } from '../../../mesh/builder/sphere';
import { Box, PerforatedBox } from '../../../primitive/box';
import { OctagonalPyramid, PerforatedOctagonalPyramid } from '../../../primitive/pyramid';
import { Star } from '../../../primitive/star';
import { Octahedron, PerforatedOctahedron } from '../../../primitive/octahedron';
import { DiamondPrism, PentagonalPrism, HexagonalPrism } from '../../../primitive/prism';

const t = Mat4.identity()
const sVec = Vec3.zero()
const pd = Vec3.zero()

const sideFactor = 1.75 * 2 * 0.806; // 0.806 == Math.cos(Math.PI / 4)
const radiusFactor = 1.75

const box = Box()
const perforatedBox = PerforatedBox()
const octagonalPyramid = OctagonalPyramid()
const perforatedOctagonalPyramid = PerforatedOctagonalPyramid()
const star = Star({ outerRadius: 1, innerRadius: 0.5, thickness: 0.5, pointCount: 5 })
const octahedron = Octahedron()
const perforatedOctahedron = PerforatedOctahedron()
const diamondPrism = DiamondPrism()
const pentagonalPrism = PentagonalPrism()
const hexagonalPrism = HexagonalPrism()

async function createCarbohydrateSymbolMesh(ctx: RuntimeContext, structure: Structure, props: CarbohydrateSymbolProps, mesh?: Mesh) {
    const builder = MeshBuilder.create(256, 128, mesh)

    const sizeTheme = SizeTheme(props.sizeTheme)
    const { detail } = props

    const carbohydrates = structure.carbohydrates
    const l = StructureElement.create()

    for (let i = 0, il = carbohydrates.elements.length; i < il; ++i) {
        const c = carbohydrates.elements[i];
        const shapeType = getSaccharideShape(c.component.type)

        l.unit = c.unit
        l.element = c.unit.elements[c.anomericCarbon]
        const size = sizeTheme.size(l)
        const radius = size * radiusFactor
        const side = size * sideFactor

        const { center, normal, direction } = c.geometry
        Vec3.add(pd, center, direction)
        Mat4.targetTo(t, center, pd, normal)
        Mat4.setTranslation(t, center)

        builder.setGroup(i * 2)

        switch (shapeType) {
            case SaccharideShapes.FilledSphere:
                addSphere(builder, center, radius, detail)
                break;
            case SaccharideShapes.FilledCube:
                Mat4.scaleUniformly(t, t, side)
                builder.add(t, box)
                break;
            case SaccharideShapes.CrossedCube:
                Mat4.scaleUniformly(t, t, side)
                builder.add(t, perforatedBox)
                Mat4.mul(t, t, Mat4.rotZ90X180)
                builder.setGroup(i * 2 + 1)
                builder.add(t, perforatedBox)
                break;
            case SaccharideShapes.FilledCone:
                Mat4.scaleUniformly(t, t, side * 1.2)
                builder.add(t, octagonalPyramid)
                break
            case SaccharideShapes.DevidedCone:
                Mat4.scaleUniformly(t, t, side * 1.2)
                builder.add(t, perforatedOctagonalPyramid)
                Mat4.mul(t, t, Mat4.rotZ90)
                builder.setGroup(i * 2 + 1)
                builder.add(t, perforatedOctagonalPyramid)
                break
            case SaccharideShapes.FlatBox:
                Mat4.mul(t, t, Mat4.rotZY90)
                Mat4.scale(t, t, Vec3.set(sVec, side, side, side / 2))
                builder.add(t, box)
                break
            case SaccharideShapes.FilledStar:
                Mat4.scaleUniformly(t, t, side)
                Mat4.mul(t, t, Mat4.rotZY90)
                builder.add(t, star)
                break
            case SaccharideShapes.FilledDiamond:
                Mat4.mul(t, t, Mat4.rotZY90)
                Mat4.scale(t, t, Vec3.set(sVec, side * 1.4, side * 1.4, side * 1.4))
                builder.add(t, octahedron)
                break
            case SaccharideShapes.DividedDiamond:
                Mat4.mul(t, t, Mat4.rotZY90)
                Mat4.scale(t, t, Vec3.set(sVec, side * 1.4, side * 1.4, side * 1.4))
                builder.add(t, perforatedOctahedron)
                Mat4.mul(t, t, Mat4.rotY90)
                builder.setGroup(i * 2 + 1)
                builder.add(t, perforatedOctahedron)
                break
            case SaccharideShapes.FlatDiamond:
                Mat4.mul(t, t, Mat4.rotZY90)
                Mat4.scale(t, t, Vec3.set(sVec, side, side / 2, side / 2))
                builder.add(t, diamondPrism)
                break
            case SaccharideShapes.Pentagon:
                Mat4.mul(t, t, Mat4.rotZY90)
                Mat4.scale(t, t, Vec3.set(sVec, side, side, side / 2))
                builder.add(t, pentagonalPrism)
                break
            case SaccharideShapes.FlatHexagon:
            default:
                Mat4.mul(t, t, Mat4.rotZYZ90)
                Mat4.scale(t, t, Vec3.set(sVec, side / 1.5, side , side / 2))
                builder.add(t, hexagonalPrism)
                break
        }
    }

    return builder.getMesh()
}

export const DefaultCarbohydrateSymbolProps = {
    ...DefaultComplexMeshProps,
    sizeTheme: { name: 'uniform', value: 1, factor: 1 } as SizeThemeProps,
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
        setUpdateState: (state: MeshUpdateState, newProps: CarbohydrateSymbolProps, currentProps: CarbohydrateSymbolProps) => {
            state.createMesh = newProps.detail !== currentProps.detail
        }
    })
}

function CarbohydrateElementIterator(structure: Structure): LocationIterator {
    const carbElements = structure.carbohydrates.elements
    const groupCount = carbElements.length * 2
    const instanceCount = 1
    const location = StructureElement.create()
    function getLocation (groupIndex: number, instanceIndex: number) {
        const carb = carbElements[Math.floor(groupIndex / 2)]
        location.unit = carb.unit
        location.element = carb.anomericCarbon
        return location
    }
    function isSecondary (elementIndex: number, instanceIndex: number) {
        return (elementIndex % 2) === 1
    }
    return LocationIterator(groupCount, instanceCount, getLocation, true, isSecondary)
}

function getCarbohydrateLoci(pickingId: PickingId, structure: Structure, id: number) {
    const { objectId, groupId } = pickingId
    if (id === objectId) {
        const carb = structure.carbohydrates.elements[Math.floor(groupId / 2)]
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