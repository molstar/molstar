/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Vec3 } from 'mol-math/linear-algebra';
import { LinkType } from 'mol-model/structure/model/types';
import { Unit, StructureElement, Structure, Link } from 'mol-model/structure';
import { ParamDefinition as PD } from 'mol-util/param-definition';
import { Mesh } from 'mol-geo/geometry/mesh/mesh';
import { MeshBuilder } from 'mol-geo/geometry/mesh/mesh-builder';
import { CylinderProps } from 'mol-geo/primitive/cylinder';
import { addFixedCountDashedCylinder, addCylinder, addDoubleCylinder } from 'mol-geo/geometry/mesh/builder/cylinder';
import { LocationIterator } from 'mol-geo/util/location-iterator';
import { VisualContext } from 'mol-repr/representation';

export const LinkCylinderParams = {
    linkScale: PD.Numeric(0.4, { min: 0, max: 1, step: 0.1 }),
    linkSpacing: PD.Numeric(1, { min: 0, max: 2, step: 0.01 }),
    radialSegments: PD.Numeric(16, { min: 3, max: 56, step: 1 }),
}
export const DefaultLinkCylinderProps = PD.getDefaultValues(LinkCylinderParams)
export type LinkCylinderProps = typeof DefaultLinkCylinderProps

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

export interface LinkCylinderMeshBuilderProps {
    linkCount: number
    referencePosition(edgeIndex: number): Vec3 | null
    position(posA: Vec3, posB: Vec3, edgeIndex: number): void
    order(edgeIndex: number): number
    flags(edgeIndex: number): LinkType
    radius(edgeIndex: number): number
}

/**
 * Each edge is included twice to allow for coloring/picking
 * the half closer to the first vertex, i.e. vertex a.
 */
export function createLinkCylinderMesh(ctx: VisualContext, linkBuilder: LinkCylinderMeshBuilderProps, props: LinkCylinderProps, mesh?: Mesh) {
    const { linkCount, referencePosition, position, order, flags, radius } = linkBuilder

    if (!linkCount) return Mesh.createEmpty(mesh)

    const { linkScale, linkSpacing, radialSegments } = props

    const vertexCountEstimate = radialSegments * 2 * linkCount * 2
    const builderState = MeshBuilder.createState(vertexCountEstimate, vertexCountEstimate / 4, mesh)

    const va = Vec3.zero()
    const vb = Vec3.zero()
    const vShift = Vec3.zero()
    const cylinderProps: CylinderProps = { radiusTop: 1, radiusBottom: 1, radialSegments }

    for (let edgeIndex = 0, _eI = linkCount; edgeIndex < _eI; ++edgeIndex) {
        position(va, vb, edgeIndex)

        const linkRadius = radius(edgeIndex)
        const o = order(edgeIndex)
        const f = flags(edgeIndex)
        builderState.currentGroup = edgeIndex

        if (LinkType.is(f, LinkType.Flag.MetallicCoordination) || LinkType.is(f, LinkType.Flag.Hydrogen)) {
            // show metall coordinations and hydrogen bonds with dashed cylinders
            cylinderProps.radiusTop = cylinderProps.radiusBottom = linkRadius / 3
            addFixedCountDashedCylinder(builderState, va, vb, 0.5, 7, cylinderProps)
        } else if (o === 2 || o === 3) {
            // show bonds with order 2 or 3 using 2 or 3 parallel cylinders
            const multiRadius = linkRadius * (linkScale / (0.5 * o))
            const absOffset = (linkRadius - multiRadius) * linkSpacing

            calculateShiftDir(vShift, va, vb, referencePosition(edgeIndex))
            Vec3.setMagnitude(vShift, vShift, absOffset)

            cylinderProps.radiusTop = cylinderProps.radiusBottom = multiRadius

            if (o === 3) addCylinder(builderState, va, vb, 0.5, cylinderProps)
            addDoubleCylinder(builderState, va, vb, 0.5, vShift, cylinderProps)
        } else {
            cylinderProps.radiusTop = cylinderProps.radiusBottom = linkRadius
            addCylinder(builderState, va, vb, 0.5, cylinderProps)
        }
    }

    return MeshBuilder.getMesh(builderState)
}

export namespace LinkIterator {
    export function fromGroup(group: Unit.SymmetryGroup): LocationIterator {
        const unit = group.units[0]
        const groupCount = Unit.isAtomic(unit) ? unit.links.edgeCount * 2 : 0
        const instanceCount = group.units.length
        const location = StructureElement.create()
        const getLocation = (groupIndex: number, instanceIndex: number) => {
            const unit = group.units[instanceIndex]
            location.unit = unit
            location.element = unit.elements[(unit as Unit.Atomic).links.a[groupIndex]]
            return location
        }
        return LocationIterator(groupCount, instanceCount, getLocation)
    }

    export function fromStructure(structure: Structure): LocationIterator {
        const groupCount = structure.links.bondCount
        const instanceCount = 1
        const location = Link.Location()
        const getLocation = (groupIndex: number) => {
            const bond = structure.links.bonds[groupIndex]
            location.aUnit = bond.unitA
            location.aIndex = bond.indexA as StructureElement.UnitIndex
            location.bUnit = bond.unitB
            location.bIndex = bond.indexB as StructureElement.UnitIndex
            return location
        }
        return LocationIterator(groupCount, instanceCount, getLocation, true)
    }
}