/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Vec3 } from '../../../../mol-math/linear-algebra';
import { BondType } from '../../../../mol-model/structure/model/types';
import { Unit, StructureElement, Structure, Bond } from '../../../../mol-model/structure';
import { ParamDefinition as PD } from '../../../../mol-util/param-definition';
import { Mesh } from '../../../../mol-geo/geometry/mesh/mesh';
import { MeshBuilder } from '../../../../mol-geo/geometry/mesh/mesh-builder';
import { CylinderProps } from '../../../../mol-geo/primitive/cylinder';
import { addFixedCountDashedCylinder, addCylinder, addDoubleCylinder } from '../../../../mol-geo/geometry/mesh/builder/cylinder';
import { LocationIterator } from '../../../../mol-geo/util/location-iterator';
import { VisualContext } from '../../../visual';
import { StructureGroup } from '../../units-visual';

export const BondCylinderParams = {
    bondScale: PD.Numeric(0.4, { min: 0, max: 1, step: 0.1 }),
    bondSpacing: PD.Numeric(1, { min: 0, max: 2, step: 0.01 }),
    bondCap: PD.Boolean(false),
    radialSegments: PD.Numeric(16, { min: 2, max: 56, step: 2 }),
    includeTypes: PD.MultiSelect(Object.keys(BondType.Names) as BondType.Names[], PD.objectToOptions(BondType.Names)),
    excludeTypes: PD.MultiSelect([] as BondType.Names[], PD.objectToOptions(BondType.Names)),
}
export const DefaultBondCylinderProps = PD.getDefaultValues(BondCylinderParams)
export type BondCylinderProps = typeof DefaultBondCylinderProps

export function ignoreBondType(include: BondType.Flag, exclude: BondType.Flag, f: BondType.Flag) {
    return !BondType.is(include, f) || BondType.is(exclude, f)
}

const tmpShiftV12 = Vec3.zero()
const tmpShiftV13 = Vec3.zero()

/** Calculate 'shift' direction that is perpendiculat to v1 - v2 and goes through v3 */
export function calculateShiftDir (out: Vec3, v1: Vec3, v2: Vec3, v3: Vec3 | null) {
    Vec3.normalize(tmpShiftV12, Vec3.sub(tmpShiftV12, v1, v2))
    if (v3 !== null) {
        Vec3.sub(tmpShiftV13, v1, v3)
    } else {
        Vec3.copy(tmpShiftV13, v1)  // no reference point, use v1
    }
    Vec3.normalize(tmpShiftV13, tmpShiftV13)

    // ensure v13 and v12 are not colinear
    let dp = Vec3.dot(tmpShiftV12, tmpShiftV13)
    if (1 - Math.abs(dp) < 1e-5) {
        Vec3.set(tmpShiftV13, 1, 0, 0)
        dp = Vec3.dot(tmpShiftV12, tmpShiftV13)
        if (1 - Math.abs(dp) < 1e-5) {
            Vec3.set(tmpShiftV13, 0, 1, 0)
            dp = Vec3.dot(tmpShiftV12, tmpShiftV13)
        }
    }

    Vec3.setMagnitude(tmpShiftV12, tmpShiftV12, dp)
    Vec3.sub(tmpShiftV13, tmpShiftV13, tmpShiftV12)
    return Vec3.normalize(out, tmpShiftV13)
}

export interface BondCylinderMeshBuilderProps {
    bondCount: number
    referencePosition(edgeIndex: number): Vec3 | null
    position(posA: Vec3, posB: Vec3, edgeIndex: number): void
    order(edgeIndex: number): number
    flags(edgeIndex: number): BondType
    radius(edgeIndex: number): number
    ignore(edgeIndex: number): boolean
}

/**
 * Each edge is included twice to allow for coloring/picking
 * the half closer to the first vertex, i.e. vertex a.
 */
export function createBondCylinderMesh(ctx: VisualContext, bondBuilder: BondCylinderMeshBuilderProps, props: BondCylinderProps, mesh?: Mesh) {
    const { bondCount, referencePosition, position, order, flags, radius, ignore } = bondBuilder

    if (!bondCount) return Mesh.createEmpty(mesh)

    const { bondScale, bondSpacing, radialSegments, bondCap } = props

    const vertexCountEstimate = radialSegments * 2 * bondCount * 2
    const builderState = MeshBuilder.createState(vertexCountEstimate, vertexCountEstimate / 4, mesh)

    const va = Vec3.zero()
    const vb = Vec3.zero()
    const vShift = Vec3.zero()
    const cylinderProps: CylinderProps = {
        radiusTop: 1,
        radiusBottom: 1,
        radialSegments,
        topCap: bondCap,
        bottomCap: bondCap
    }

    for (let edgeIndex = 0, _eI = bondCount; edgeIndex < _eI; ++edgeIndex) {
        if (ignore(edgeIndex)) continue

        position(va, vb, edgeIndex)

        const linkRadius = radius(edgeIndex)
        const o = order(edgeIndex)
        const f = flags(edgeIndex)
        builderState.currentGroup = edgeIndex

        if (BondType.is(f, BondType.Flag.MetallicCoordination) || BondType.is(f, BondType.Flag.HydrogenBond)) {
            // show metall coordinations and hydrogen bonds with dashed cylinders
            cylinderProps.radiusTop = cylinderProps.radiusBottom = linkRadius / 3
            cylinderProps.topCap = cylinderProps.bottomCap = true
            addFixedCountDashedCylinder(builderState, va, vb, 0.5, 7, cylinderProps)
        } else if (o === 2 || o === 3) {
            // show bonds with order 2 or 3 using 2 or 3 parallel cylinders
            const multiRadius = linkRadius * (bondScale / (0.5 * o))
            const absOffset = (linkRadius - multiRadius) * bondSpacing

            calculateShiftDir(vShift, va, vb, referencePosition(edgeIndex))
            Vec3.setMagnitude(vShift, vShift, absOffset)

            cylinderProps.radiusTop = cylinderProps.radiusBottom = multiRadius
            cylinderProps.topCap = cylinderProps.bottomCap = bondCap

            if (o === 3) addCylinder(builderState, va, vb, 0.5, cylinderProps)
            addDoubleCylinder(builderState, va, vb, 0.5, vShift, cylinderProps)
        } else {
            cylinderProps.radiusTop = cylinderProps.radiusBottom = linkRadius
            cylinderProps.topCap = cylinderProps.bottomCap = bondCap
            addCylinder(builderState, va, vb, 0.5, cylinderProps)
        }
    }

    return MeshBuilder.getMesh(builderState)
}

export namespace BondIterator {
    export function fromGroup(structureGroup: StructureGroup): LocationIterator {
        const { group } = structureGroup
        const unit = group.units[0]
        const groupCount = Unit.isAtomic(unit) ? unit.bonds.edgeCount * 2 : 0
        const instanceCount = group.units.length
        const location = StructureElement.Location.create()
        const getLocation = (groupIndex: number, instanceIndex: number) => {
            const unit = group.units[instanceIndex]
            location.unit = unit
            location.element = unit.elements[(unit as Unit.Atomic).bonds.a[groupIndex]]
            return location
        }
        return LocationIterator(groupCount, instanceCount, getLocation)
    }

    export function fromStructure(structure: Structure): LocationIterator {
        const groupCount = structure.interUnitBonds.edgeCount
        const instanceCount = 1
        const location = Bond.Location()
        const getLocation = (groupIndex: number) => {
            const bond = structure.interUnitBonds.edges[groupIndex]
            location.aUnit = bond.unitA
            location.aIndex = bond.indexA as StructureElement.UnitIndex
            location.bUnit = bond.unitB
            location.bIndex = bond.indexB as StructureElement.UnitIndex
            return location
        }
        return LocationIterator(groupCount, instanceCount, getLocation, true)
    }
}