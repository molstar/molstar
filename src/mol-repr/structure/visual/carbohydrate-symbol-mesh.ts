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
import { VisualContext } from 'mol-repr/representation';
import { Theme } from 'mol-theme/theme';
import { getResidueLoci } from './util/common';

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

async function createCarbohydrateSymbolMesh(ctx: VisualContext, structure: Structure, theme: Theme, props: PD.Values<CarbohydrateSymbolParams>, mesh?: Mesh) {
    const builder = MeshBuilder.create(256, 128, mesh)

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

        if (i % 10000 === 0 && ctx.runtime.shouldUpdate) {
            await ctx.runtime.update({ message: 'Carbohydrate symbols', current: i, max: n });
        }
    }

    return builder.getMesh()
}

export const CarbohydrateSymbolParams = {
    ...ComplexMeshParams,
    detail: PD.Numeric(0, { min: 0, max: 3, step: 1 }),
    sizeFactor: PD.Numeric(1.75, { min: 0, max: 10, step: 0.01 }),
}
export type CarbohydrateSymbolParams = typeof CarbohydrateSymbolParams

export function CarbohydrateSymbolVisual(): ComplexVisual<CarbohydrateSymbolParams> {
    return ComplexMeshVisual<CarbohydrateSymbolParams>({
        defaultProps: PD.getDefaultValues(CarbohydrateSymbolParams),
        createGeometry: createCarbohydrateSymbolMesh,
        createLocationIterator: CarbohydrateElementIterator,
        getLoci: getCarbohydrateLoci,
        mark: markCarbohydrate,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<CarbohydrateSymbolParams>, currentProps: PD.Values<CarbohydrateSymbolParams>) => {
            state.createGeometry = (
                newProps.sizeFactor !== currentProps.sizeFactor ||
                newProps.detail !== currentProps.detail
            )
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

/** Return a Loci for the elements of the whole residue of a carbohydrate. */
function getCarbohydrateLoci(pickingId: PickingId, structure: Structure, id: number) {
    const { objectId, groupId } = pickingId
    if (id === objectId) {
        const carb = structure.carbohydrates.elements[Math.floor(groupId / 2)]
        return getResidueLoci(structure, carb.unit, carb.anomericCarbon)
    }
    return EmptyLoci
}

/** Mark a carbohydrate (usually a monosaccharide) when all its residue's elements are in a loci. */
function markCarbohydrate(loci: Loci, structure: Structure, apply: (interval: Interval) => boolean) {
    const { getElementIndex, getAnomericCarbon } = structure.carbohydrates
    let changed = false
    if (StructureElement.isLoci(loci)) {
        for (const e of loci.elements) {
            OrderedSet.forEach(e.indices, v => {
                const { model, elements } = e.unit
                const { index, offsets } = model.atomicHierarchy.residueAtomSegments        
                const rI = index[elements[v]]
                const unitIndexMin = OrderedSet.findPredecessorIndex(elements, offsets[rI])
                const unitIndexMax = OrderedSet.findPredecessorIndex(elements, offsets[rI + 1] - 1)
                const unitIndexInterval = Interval.ofRange(unitIndexMin, unitIndexMax)
                if(!OrderedSet.isSubset(e.indices, unitIndexInterval)) return
                const eI = getAnomericCarbon(e.unit, rI)
                if (eI !== undefined) {
                    const idx = getElementIndex(e.unit, eI)
                    if (idx !== undefined) {
                        if (apply(Interval.ofBounds(idx * 2, idx * 2 + 2))) changed = true
                    }
                }
            })
        }
    }
    return changed
}