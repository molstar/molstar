/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Vec3 } from 'mol-math/linear-algebra';
import { RuntimeContext } from 'mol-task';
import { Mesh } from '../../../../shape/mesh';
import { MeshBuilder } from '../../../../shape/mesh-builder';

export const DefaultLinkCylinderProps = {
    linkScale: 0.4,
    linkSpacing: 1,
    linkRadius: 0.25,
    radialSegments: 16
}
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
}

/**
 * Each edge is included twice to allow for coloring/picking
 * the half closer to the first vertex, i.e. vertex a.
 */
export async function createLinkCylinderMesh(ctx: RuntimeContext, linkBuilder: LinkCylinderMeshBuilderProps, props: LinkCylinderProps, mesh?: Mesh) {
    const { linkCount, referencePosition, position, order } = linkBuilder

    if (!linkCount) return Mesh.createEmpty(mesh)

    // approximate vertextCount, exact calculation would need to take link orders into account
    const vertexCount = 32 * linkCount
    const meshBuilder = MeshBuilder.create(vertexCount, vertexCount / 2, mesh)

    const va = Vec3.zero()
    const vb = Vec3.zero()
    const vShift = Vec3.zero()

    const { linkScale, linkSpacing, linkRadius, radialSegments } = props

    const cylinderParams = {
        height: 1,
        radiusTop: linkRadius,
        radiusBottom: linkRadius,
        radialSegments
    }

    for (let edgeIndex = 0, _eI = linkCount; edgeIndex < _eI; ++edgeIndex) {
        position(va, vb, edgeIndex)

        const o = order(edgeIndex)
        meshBuilder.setId(edgeIndex)

        if (o === 2 || o === 3) {
            const multiRadius = linkRadius * (linkScale / (0.5 * o))
            const absOffset = (linkRadius - multiRadius) * linkSpacing

            calculateShiftDir(vShift, va, vb, referencePosition(edgeIndex))
            Vec3.setMagnitude(vShift, vShift, absOffset)

            cylinderParams.radiusTop = cylinderParams.radiusBottom = multiRadius

            if (o === 3) meshBuilder.addCylinder(va, vb, 0.5, cylinderParams)
            meshBuilder.addDoubleCylinder(va, vb, 0.5, vShift, cylinderParams)
        } else {
            cylinderParams.radiusTop = cylinderParams.radiusBottom = linkRadius
            meshBuilder.addCylinder(va, vb, 0.5, cylinderParams)
        }

        if (edgeIndex % 10000 === 0 && ctx.shouldUpdate) {
            await ctx.update({ message: 'Cylinder mesh', current: edgeIndex, max: linkCount });
        }
    }

    return meshBuilder.getMesh()
}