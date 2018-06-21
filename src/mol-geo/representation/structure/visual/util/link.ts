/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Vec3, Mat4 } from 'mol-math/linear-algebra';
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
    Vec3.sub(tmpShiftV12, v1, v2)

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

export interface LinkCylinderMeshBuilder {
    linkCount: number
    eachLink(ctx: RuntimeContext, cb: (edgeIndex: number, aI: number, bI: number) => void): Promise<void>
    // assumes aI < bI
    getRefPos(aI: number, bI: number): Vec3 | null
    setPositions(posA: Vec3, posB: Vec3, edgeIndex: number): void
    getOrder(edgeIndex: number): number
}

/**
 * Each edge is included twice to allow for coloring/picking
 * the half closer to the first vertex, i.e. vertex a.
 */
export async function createLinkCylinderMesh(ctx: RuntimeContext, linkBuilder: LinkCylinderMeshBuilder, props: LinkCylinderProps, mesh?: Mesh) {
    const { linkCount, eachLink, getRefPos, setPositions, getOrder } = linkBuilder

    if (!linkCount) return Mesh.createEmpty(mesh)

    // approximate vertextCount, exact calculation would need to take link orders into account
    const vertexCount = 32 * linkCount
    const meshBuilder = MeshBuilder.create(vertexCount, vertexCount / 2, mesh)

    const va = Vec3.zero()
    const vb = Vec3.zero()
    const vd = Vec3.zero()
    const vc = Vec3.zero()
    const m = Mat4.identity()
    const mt = Mat4.identity()

    const vShift = Vec3.zero()
    const vCenter = Vec3.zero()

    const { linkScale, linkSpacing, linkRadius, radialSegments } = props

    const cylinderParams = {
        height: 1,
        radiusTop: linkRadius,
        radiusBottom: linkRadius,
        radialSegments,
        openEnded: true
    }

    // for (let edgeIndex = 0, _eI = edgeCount * 2; edgeIndex < _eI; ++edgeIndex) {
    await eachLink(ctx, (edgeIndex, aI, bI) => {
        // const aI = a[edgeIndex], bI = b[edgeIndex];

        setPositions(va, vb, edgeIndex)
        const d = Vec3.distance(va, vb)

        Vec3.sub(vd, vb, va)
        Vec3.scale(vd, Vec3.normalize(vd, vd), d / 4)
        Vec3.add(vc, va, vd)
        // ensure both edge halfs are pointing in the the same direction so the triangles align
        if (aI > bI) Vec3.scale(vd, vd, -1)
        Vec3.makeRotation(m, Vec3.create(0, 1, 0), vd)

        const order = getOrder(edgeIndex)
        meshBuilder.setId(edgeIndex)
        cylinderParams.height = d / 2

        if (order === 2 || order === 3) {
            const multiRadius = linkRadius * (linkScale / (0.5 * order))
            const absOffset = (linkRadius - multiRadius) * linkSpacing

            if (aI < bI) {
                calculateShiftDir(vShift, va, vb, getRefPos(aI, bI))
            } else {
                calculateShiftDir(vShift, vb, va, getRefPos(bI, aI))
            }
            Vec3.setMagnitude(vShift, vShift, absOffset)

            cylinderParams.radiusTop = multiRadius
            cylinderParams.radiusBottom = multiRadius

            if (order === 3) {
                Mat4.fromTranslation(mt, vc)
                Mat4.mul(mt, mt, m)
                meshBuilder.addCylinder(mt, cylinderParams)
            }

            Vec3.add(vCenter, vc, vShift)
            Mat4.fromTranslation(mt, vCenter)
            Mat4.mul(mt, mt, m)
            meshBuilder.addCylinder(mt, cylinderParams)

            Vec3.sub(vCenter, vc, vShift)
            Mat4.fromTranslation(mt, vCenter)
            Mat4.mul(mt, mt, m)
            meshBuilder.addCylinder(mt, cylinderParams)
        } else {
            cylinderParams.radiusTop = linkRadius
            cylinderParams.radiusBottom = linkRadius

            Mat4.setTranslation(m, vc)
            meshBuilder.addCylinder(m, cylinderParams)
        }
    })

    return meshBuilder.getMesh()
}