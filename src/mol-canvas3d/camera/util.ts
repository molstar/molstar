/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Mat4, Vec3, Vec4, EPSILON } from '../../mol-math/linear-algebra'
import { Camera } from '../camera'
import { Sphere3D } from '../../mol-math/geometry'
import { Canvas3DProps } from '../canvas3d'

export { Viewport }

type Viewport = {
    x: number
    y: number
    width: number
    height: number
}

function Viewport() {
    return Viewport.zero()
}

namespace Viewport {
    export function zero(): Viewport {
        return { x: 0, y: 0, width: 0, height: 0 }
    }
    export function create(x: number, y: number, width: number, height: number): Viewport {
        return { x, y, width, height }
    }
    export function clone(viewport: Viewport): Viewport {
        return { ...viewport }
    }
    export function copy(target: Viewport, source: Viewport): Viewport {
        return Object.assign(target, source)
    }
    export function set(viewport: Viewport, x: number, y: number, width: number, height: number): Viewport {
        viewport.x = x
        viewport.y = y
        viewport.width = width
        viewport.height = height
        return viewport
    }

    export function toVec4(v4: Vec4, viewport: Viewport): Vec4 {
        v4[0] = viewport.x
        v4[1] = viewport.y
        v4[2] = viewport.width
        v4[3] = viewport.height
        return v4
    }

    export function equals(a: Viewport, b: Viewport) {
        return a.x === b.x && a.y === b.y && a.width === b.width && a.height === b.height
    }
}

//

const tmpVec3 = Vec3()

/** Modifies the direction & up vectors in place, both are normalized */
export function cameraLookAt(position: Vec3, up: Vec3, direction: Vec3, target: Vec3) {
    Vec3.sub(tmpVec3, target, position)
    Vec3.normalize(tmpVec3, tmpVec3)

    if (!Vec3.isZero(tmpVec3)) {
        // change direction vector to look at target
        const d = Vec3.dot(tmpVec3, up)
        if (Math.abs(d - 1) < EPSILON) { // parallel
            Vec3.scale(up, direction, -1)
        } else if (Math.abs(d + 1) < EPSILON) { // anti parallel
            Vec3.copy(up, direction)
        }
        Vec3.copy(direction, tmpVec3)

        // normalize up vector
        Vec3.cross(tmpVec3, direction, up)
        Vec3.normalize(tmpVec3, tmpVec3)
        Vec3.cross(up, tmpVec3, direction)
        Vec3.normalize(up, up)
    }
}

export function cameraSetClipping(state: Camera.Snapshot, boundingSphere: Sphere3D, p: Canvas3DProps) {
    const cDist = Vec3.distance(state.position, state.target)
    const bRadius = Math.max(10, boundingSphere.radius)

    const nearFactor = (50 - p.clip[0]) / 50
    const farFactor = -(50 - p.clip[1]) / 50
    let near = cDist - (bRadius * nearFactor)
    let far = cDist + (bRadius * farFactor)

    const fogNearFactor = (50 - p.fog[0]) / 50
    const fogFarFactor = -(50 - p.fog[1]) / 50
    let fogNear = cDist - (bRadius * fogNearFactor)
    let fogFar = cDist + (bRadius * fogFarFactor)

    if (state.mode === 'perspective') {
        // set at least to 5 to avoid slow sphere impostor rendering
        near = Math.max(5, p.cameraClipDistance, near)
        far = Math.max(5, far)
        fogNear = Math.max(5, fogNear)
        fogFar = Math.max(5, fogFar)
    } else if (state.mode === 'orthographic') {
        if (p.cameraClipDistance > 0) {
            near = Math.max(p.cameraClipDistance, near)
        }
    }

    state.near = near;
    state.far = far;
    state.fogNear = fogNear;
    state.fogFar = fogFar;

    // if (near !== currentNear || far !== currentFar || fogNear !== currentFogNear || fogFar !== currentFogFar) {
    //     camera.setState({ near, far, fogNear, fogFar })
    //     currentNear = near, currentFar = far, currentFogNear = fogNear, currentFogFar = fogFar
    // }
}

const NEAR_RANGE = 0
const FAR_RANGE = 1

const tmpVec4 = Vec4()

/** Transform point into 2D window coordinates. */
export function cameraProject (out: Vec4, point: Vec3, viewport: Viewport, projectionView: Mat4) {
    const { x: vX, y: vY, width: vWidth, height: vHeight } = viewport

    // clip space -> NDC -> window coordinates, implicit 1.0 for w component
    Vec4.set(tmpVec4, point[0], point[1], point[2], 1.0)

    // transform into clip space
    Vec4.transformMat4(tmpVec4, tmpVec4, projectionView)

    // transform into NDC
    const w = tmpVec4[3]
    if (w !== 0) {
        tmpVec4[0] /= w
        tmpVec4[1] /= w
        tmpVec4[2] /= w
    }

    // transform into window coordinates, set fourth component is (1/clip.w) as in gl_FragCoord.w
    out[0] = vX + vWidth / 2 * tmpVec4[0] + (0 + vWidth / 2)
    out[1] = vY + vHeight / 2 * tmpVec4[1] + (0 + vHeight / 2)
    out[2] = (FAR_RANGE - NEAR_RANGE) / 2 * tmpVec4[2] + (FAR_RANGE + NEAR_RANGE) / 2
    out[3] = w === 0 ? 0 : 1 / w
    return out
}

/**
 * Transform point from screen space to 3D coordinates.
 * The point must have x and y set to 2D window coordinates and z between 0 (near) and 1 (far).
 */
export function cameraUnproject (out: Vec3, point: Vec3, viewport: Viewport, inverseProjectionView: Mat4) {
    const { x: vX, y: vY, width: vWidth, height: vHeight } = viewport

    const x = point[0] - vX
    const y = (vHeight - point[1] - 1) - vY
    const z = point[2]

    out[0] = (2 * x) / vWidth - 1
    out[1] = (2 * y) / vHeight - 1
    out[2] = 2 * z - 1
    return Vec3.transformMat4(out, out, inverseProjectionView)
}