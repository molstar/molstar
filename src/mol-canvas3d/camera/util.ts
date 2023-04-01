/**
 * Copyright (c) 2018-2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Mat4, Vec3, Vec4 } from '../../mol-math/linear-algebra';

export { Viewport };

type Viewport = {
    x: number
    y: number
    width: number
    height: number
}

function Viewport() {
    return Viewport.zero();
}

namespace Viewport {
    export function zero(): Viewport {
        return { x: 0, y: 0, width: 0, height: 0 };
    }
    export function create(x: number, y: number, width: number, height: number): Viewport {
        return { x, y, width, height };
    }
    export function clone(viewport: Viewport): Viewport {
        return { ...viewport };
    }
    export function copy(target: Viewport, source: Viewport): Viewport {
        return Object.assign(target, source);
    }
    export function set(viewport: Viewport, x: number, y: number, width: number, height: number): Viewport {
        viewport.x = x;
        viewport.y = y;
        viewport.width = width;
        viewport.height = height;
        return viewport;
    }

    export function toVec4(v4: Vec4, viewport: Viewport): Vec4 {
        v4[0] = viewport.x;
        v4[1] = viewport.y;
        v4[2] = viewport.width;
        v4[3] = viewport.height;
        return v4;
    }

    export function equals(a: Viewport, b: Viewport) {
        return a.x === b.x && a.y === b.y && a.width === b.width && a.height === b.height;
    }
}

//

const tmpVec4 = Vec4();

/** Transform point into 2D window coordinates. */
export function cameraProject(out: Vec4, point: Vec3, viewport: Viewport, projectionView: Mat4) {
    const { x, y, width, height } = viewport;

    // clip space -> NDC -> window coordinates, implicit 1.0 for w component
    Vec4.set(tmpVec4, point[0], point[1], point[2], 1.0);

    // transform into clip space
    Vec4.transformMat4(tmpVec4, tmpVec4, projectionView);

    // transform into NDC
    const w = tmpVec4[3];
    if (w !== 0) {
        tmpVec4[0] /= w;
        tmpVec4[1] /= w;
        tmpVec4[2] /= w;
    }

    // transform into window coordinates, set fourth component to 1 / clip.w as in gl_FragCoord.w
    out[0] = (tmpVec4[0] + 1) * width * 0.5 + x;
    out[1] = (tmpVec4[1] + 1) * height * 0.5 + y;
    out[2] = (tmpVec4[2] + 1) * 0.5;
    out[3] = w === 0 ? 0 : 1 / w;
    return out;
}

/**
 * Transform point from screen space to 3D coordinates.
 * The point must have `x` and `y` set to 2D window coordinates
 * and `z` between 0 (near) and 1 (far); the optional `w` is not used.
 */
export function cameraUnproject(out: Vec3, point: Vec3 | Vec4, viewport: Viewport, inverseProjectionView: Mat4) {
    const { x, y, width, height } = viewport;

    const px = point[0] - x;
    const py = point[1] - y;
    const pz = point[2];

    out[0] = (2 * px) / width - 1;
    out[1] = (2 * py) / height - 1;
    out[2] = 2 * pz - 1;
    return Vec3.transformMat4(out, out, inverseProjectionView);
}