/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Mat4, Vec3 } from 'mol-math/linear-algebra';
import { Box, PerforatedBox } from 'mol-geo/primitive/box';
import { OctagonalPyramid, PerforatedOctagonalPyramid } from 'mol-geo/primitive/pyramid';
import { Star } from 'mol-geo/primitive/star';
import { Octahedron, PerforatedOctahedron } from 'mol-geo/primitive/octahedron';
import { DiamondPrism, PentagonalPrism, HexagonalPrism } from 'mol-geo/primitive/prism';
import { Structure, StructureElement } from 'mol-model/structure';
import { Mesh } from 'mol-geo/geometry/mesh/mesh';
import { MeshBuilder } from 'mol-geo/geometry/mesh/mesh-builder';
import { getSaccharideShape, SaccharideShapes } from 'mol-model/structure/structure/carbohydrates/constants';
import { addSphere } from 'mol-geo/geometry/mesh/builder/sphere';
import { ComplexMeshParams, ComplexMeshVisual } from '../complex-visual';
import { ParamDefinition as PD } from 'mol-util/param-definition';
import { ComplexVisual } from '../representation';
import { VisualUpdateState } from '../../util';
import { LocationIterator } from 'mol-geo/util/location-iterator';
import { PickingId } from 'mol-geo/geometry/picking';
import { OrderedSet, Interval } from 'mol-data/int';
import { EmptyLoci, Loci } from 'mol-model/loci';
import { VisualContext } from 'mol-repr/visual';
import { Theme } from 'mol-theme/theme';
import { getAltResidueLoci } from './util/common';

const t = Mat4.identity()
const sVec = Vec3.zero()
const pd = Vec3.zero()

const SideFactor = 2 * 0.806; // 0.806 == Math.cos(Math.PI / 4)

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

function createCarbohydrateSymbolMesh(ctx: VisualContext, structure: Structure, theme: Theme, props: PD.Values<CarbohydrateSymbolParams>, mesh?: Mesh) {
    const builderState = MeshBuilder.createState(256, 128, mesh)

    const { detail, sizeFactor } = props

    const carbohydrates = structure.carbohydrates
    const n = carbohydrates.elements.length
    const l = StructureElement.create()

    for (let i = 0; i < n; ++i) {
        const c = carbohydrates.elements[i];
        const shapeType = getSaccharideShape(c.component.type)

        l.unit = c.unit
        l.element = c.unit.elements[c.anomericCarbon]
        const size = theme.size.size(l)
        const radius = size * sizeFactor
        const side = size * sizeFactor * SideFactor

        const { center, normal, direction } = c.geometry
        Vec3.add(pd, center, direction)
        Mat4.targetTo(t, center, pd, normal)
        Mat4.setTranslation(t, center)

        builderState.currentGroup = i * 2

        switch (shapeType) {
            case SaccharideShapes.FilledSphere:
                addSphere(builderState, center, radius, detail)
                break;
            case SaccharideShapes.FilledCube:
                Mat4.scaleUniformly(t, t, side)
                MeshBuilder.addPrimitive(builderState, t, box)
                break;
            case SaccharideShapes.CrossedCube:
                Mat4.scaleUniformly(t, t, side)
                MeshBuilder.addPrimitive(builderState, t, perforatedBox)
                Mat4.mul(t, t, Mat4.rotZ90X180)
                builderState.currentGroup += 1
                MeshBuilder.addPrimitive(builderState, t, perforatedBox)
                break;
            case SaccharideShapes.FilledCone:
                Mat4.scaleUniformly(t, t, side * 1.2)
                MeshBuilder.addPrimitive(builderState, t, octagonalPyramid)
                break
            case SaccharideShapes.DevidedCone:
                Mat4.scaleUniformly(t, t, side * 1.2)
                MeshBuilder.addPrimitive(builderState, t, perforatedOctagonalPyramid)
                Mat4.mul(t, t, Mat4.rotZ90)
                builderState.currentGroup += 1
                MeshBuilder.addPrimitive(builderState, t, perforatedOctagonalPyramid)
                break
            case SaccharideShapes.FlatBox:
                Mat4.mul(t, t, Mat4.rotZY90)
                Mat4.scale(t, t, Vec3.set(sVec, side, side, side / 2))
                MeshBuilder.addPrimitive(builderState, t, box)
                break
            case SaccharideShapes.FilledStar:
                Mat4.scaleUniformly(t, t, side)
                Mat4.mul(t, t, Mat4.rotZY90)
                MeshBuilder.addPrimitive(builderState, t, star)
                break
            case SaccharideShapes.FilledDiamond:
                Mat4.mul(t, t, Mat4.rotZY90)
                Mat4.scale(t, t, Vec3.set(sVec, side * 1.4, side * 1.4, side * 1.4))
                MeshBuilder.addPrimitive(builderState, t, octahedron)
                break
            case SaccharideShapes.DividedDiamond:
                Mat4.mul(t, t, Mat4.rotZY90)
                Mat4.scale(t, t, Vec3.set(sVec, side * 1.4, side * 1.4, side * 1.4))
                MeshBuilder.addPrimitive(builderState, t, perforatedOctahedron)
                Mat4.mul(t, t, Mat4.rotY90)
                builderState.currentGroup += 1
                MeshBuilder.addPrimitive(builderState, t, perforatedOctahedron)
                break
            case SaccharideShapes.FlatDiamond:
                Mat4.mul(t, t, Mat4.rotZY90)
                Mat4.scale(t, t, Vec3.set(sVec, side, side / 2, side / 2))
                MeshBuilder.addPrimitive(builderState, t, diamondPrism)
                break
            case SaccharideShapes.Pentagon:
                Mat4.mul(t, t, Mat4.rotZY90)
                Mat4.scale(t, t, Vec3.set(sVec, side, side, side / 2))
                MeshBuilder.addPrimitive(builderState, t, pentagonalPrism)
                break
            case SaccharideShapes.FlatHexagon:
            default:
                Mat4.mul(t, t, Mat4.rotZYZ90)
                Mat4.scale(t, t, Vec3.set(sVec, side / 1.5, side , side / 2))
                MeshBuilder.addPrimitive(builderState, t, hexagonalPrism)
                break
        }
    }

    return MeshBuilder.getMesh(builderState)
}

export const CarbohydrateSymbolParams = {
    ...ComplexMeshParams,
    detail: PD.Numeric(0, { min: 0, max: 3, step: 1 }),
    sizeFactor: PD.Numeric(1.75, { min: 0, max: 10, step: 0.01 }),
}
export type CarbohydrateSymbolParams = typeof CarbohydrateSymbolParams

export function CarbohydrateSymbolVisual(materialId: number): ComplexVisual<CarbohydrateSymbolParams> {
    return ComplexMeshVisual<CarbohydrateSymbolParams>({
        defaultProps: PD.getDefaultValues(CarbohydrateSymbolParams),
        createGeometry: createCarbohydrateSymbolMesh,
        createLocationIterator: CarbohydrateElementIterator,
        getLoci: getCarbohydrateLoci,
        eachLocation: eachCarbohydrate,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<CarbohydrateSymbolParams>, currentProps: PD.Values<CarbohydrateSymbolParams>) => {
            state.createGeometry = (
                newProps.sizeFactor !== currentProps.sizeFactor ||
                newProps.detail !== currentProps.detail
            )
        }
    }, materialId)
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

/** Return a Loci for the elements of the whole residue of a carbohydrate. */
function getCarbohydrateLoci(pickingId: PickingId, structure: Structure, id: number) {
    const { objectId, groupId } = pickingId
    if (id === objectId) {
        const carb = structure.carbohydrates.elements[Math.floor(groupId / 2)]
        return getAltResidueLoci(structure, carb.unit, carb.anomericCarbon)
    }
    return EmptyLoci
}

/** For each carbohydrate (usually a monosaccharide) when all its residue's elements are in a loci. */
function eachCarbohydrate(loci: Loci, structure: Structure, apply: (interval: Interval) => boolean) {
    const { getElementIndex, getAnomericCarbons } = structure.carbohydrates
    let changed = false
    if (!StructureElement.isLoci(loci)) return false
    if (!Structure.areEquivalent(loci.structure, structure)) return false
    for (const e of loci.elements) {
        // TODO make more efficient by handling/grouping `e.indices` by residue index
        // TODO only call apply when the full alt-residue of the unit is part of `e`
        OrderedSet.forEach(e.indices, v => {
            const { model, elements } = e.unit
            const { index } = model.atomicHierarchy.residueAtomSegments
            const rI = index[elements[v]]
            const eIndices = getAnomericCarbons(e.unit, rI)
            for (let i = 0, il = eIndices.length; i < il; ++i) {
                const eI = eIndices[i]
                if (!OrderedSet.has(e.indices, OrderedSet.indexOf(elements, eI))) continue
                const idx = getElementIndex(e.unit, eI)
                if (idx !== undefined) {
                    if (apply(Interval.ofBounds(idx * 2, idx * 2 + 2))) changed = true
                }
            }
        })
    }
    return changed
}