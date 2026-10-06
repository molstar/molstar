/**
 * Copyright (c) 2021-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Gianluca Tomasello <giagitom@gmail.com>
 */

import { EPSILON } from '../mol-math/linear-algebra/3d/common';
import { Mat4 } from '../mol-math/linear-algebra/3d/mat4';
import { Quat } from '../mol-math/linear-algebra/3d/quat';
import { Vec3 } from '../mol-math/linear-algebra/3d/vec3';
import { Vec4 } from '../mol-math/linear-algebra/3d/vec4';
import { Plane3D } from '../mol-math/geometry/primitives/plane3d';
import { Sphere3D } from '../mol-math/geometry/primitives/sphere3d';
import { degToRad } from '../mol-math/misc';
import { ParamDefinition as PD } from './param-definition';
import { stringToWords } from './string';

export interface Clip {
    variant: Clip.Variant,
    objects: Clip.Objects
}

export function Clip() {

}

export namespace Clip {
    /** Clip object types */
    export const Type = {
        none: 0, // to switch clipping off
        plane: 1,
        sphere: 2,
        cube: 3,
        cylinder: 4,
        infiniteCone: 5,
    };

    export type Variant = 'instance' | 'pixel'

    export type Objects = {
        count: number
        type: number[]
        invert: boolean[]
        position: number[]
        rotation: number[]
        scale: number[]
        /** Transform point by this before testing */
        transform: number[]
    }

    export const Params = {
        variant: PD.Select('pixel', PD.arrayToOptions<Variant>(['instance', 'pixel'])),
        objects: PD.ObjectList({
            type: PD.Select('plane', PD.objectToOptions(Type, t => stringToWords(t))),
            invert: PD.Boolean(false),
            position: PD.Vec3(Vec3()),
            rotation: PD.Group({
                axis: PD.Vec3(Vec3.create(1, 0, 0)),
                angle: PD.Numeric(0, { min: -180, max: 180, step: 1 }, { description: 'Angle in Degrees' }),
            }, { isExpanded: true }),
            scale: PD.Vec3(Vec3.create(1, 1, 1)),
            transform: PD.Mat4(Mat4.identity()),
        }, o => stringToWords(o.type))
    };
    export type Params = typeof Params
    export type Props = PD.Values<Params>

    function createClipObjects(count: number) {
        return {
            count: 0,
            type: (new Array(count)).fill(1),
            invert: (new Array(count)).fill(false),
            position: (new Array(count * 3)).fill(0),
            rotation: (new Array(count * 4)).fill(0),
            scale: (new Array(count * 3)).fill(1),
            transform: (new Array(count * 16)).fill(0),
        };
    }

    const qA = Quat();
    const qB = Quat();
    const vA = Vec3();
    const vB = Vec3();
    const mA = Mat4();
    const mB = Mat4();

    export function getClip(props: Props, clip?: Clip): Clip {
        const count = props.objects.length;
        const { type, invert, position, rotation, scale, transform } = clip?.objects || createClipObjects(count);
        for (let i = 0; i < count; ++i) {
            const p = props.objects[i];
            type[i] = Type[p.type];
            invert[i] = p.invert;
            Vec3.toArray(p.position, position, i * 3);
            Vec3.normalize(vA, p.rotation.axis);
            Quat.toArray(Quat.setAxisAngle(qA, vA, degToRad(p.rotation.angle)), rotation, i * 4);
            Vec3.toArray(p.scale, scale, i * 3);
            Mat4.toArray(p.transform, transform, i * 16);
        }
        return {
            variant: props.variant,
            objects: { count, type, invert, position, rotation, scale, transform }
        };
    }

    const v4A = Vec4();
    const pA = Plane3D();

    /**
     * Plane of the clip object at `index`, in the space of the points
     * given to the clip test (i.e. with `transform` undone) and with
     * the clipped half-space on its positive side
     */
    export function getPlane(out: Plane3D, objects: Objects, index: number): Plane3D {
        const { invert, position, rotation, transform } = objects;
        Vec3.transformQuat(vA, Vec3.unitY, Quat.fromArray(qA, rotation, index * 4));
        Vec4.set(v4A, vA[0], vA[1], vA[2], -Vec3.dot(vA, Vec3.fromArray(vB, position, index * 3)));
        if (invert[index]) Vec4.scale(v4A, v4A, -1);
        Mat4.fromArray(mA, transform, index * 16);
        Vec4.transformMat4(v4A, v4A, Mat4.transpose(mA, mA));
        return Plane3D.setUnnormalized(out, v4A[0], v4A[1], v4A[2], v4A[3]);
    }

    /**
     * Test if the surface of the clip object at `index` can intersect `sphere`,
     * exact for planes, conservative for spheres, cubes and cylinders
     */
    export function canIntersectSphere(objects: Objects, index: number, sphere: Sphere3D): boolean {
        const { type, position, scale, transform } = objects;
        const t = type[index];
        if (t === Type.plane) {
            return Math.abs(Plane3D.distanceToPoint(getPlane(pA, objects, index), sphere.center)) <= sphere.radius;
        }
        if (t !== Type.sphere && t !== Type.cube && t !== Type.cylinder) return true;

        Mat4.fromArray(mA, transform, index * 16);
        Vec3.transformMat4(vA, sphere.center, mA);
        const radius = sphere.radius * Mat4.getMaxScaleOnAxis(mA);
        const distance = Vec3.distance(vA, Vec3.fromArray(vB, position, index * 3));
        const sx = scale[index * 3], sy = scale[index * 3 + 1], sz = scale[index * 3 + 2];
        const inner = (t === Type.cylinder ? Math.min(sx, sy) : Math.min(sx, sy, sz)) / 2;
        const outer = (t === Type.sphere ? Math.max(sx, sy, sz) : t === Type.cube ? Math.hypot(sx, sy, sz) : Math.hypot(sx, sy)) / 2;
        return distance <= radius + outer && distance + radius >= inner;
    }

    export function areEqual(cA: Clip, cB: Clip) {
        if (cA.variant !== cB.variant) return false;
        if (cA.objects.count !== cB.objects.count) return false;

        const oA = cA.objects, oB = cB.objects;
        for (let i = 0, il = oA.count; i < il; ++i) {
            if (oA.invert[i] !== oB.invert[i]) return false;
            if (oA.type[i] !== oB.type[i]) return false;

            Vec3.fromArray(vA, oA.position, i * 3);
            Vec3.fromArray(vB, oB.position, i * 3);
            if (!Vec3.equals(vA, vB)) return false;

            Vec3.fromArray(vA, oA.scale, i * 3);
            Vec3.fromArray(vB, oB.scale, i * 3);
            if (!Vec3.equals(vA, vB)) return false;

            Quat.fromArray(qA, oA.rotation, i * 4);
            Quat.fromArray(qB, oB.rotation, i * 4);
            if (!Quat.equals(qA, qB)) return false;

            Mat4.fromArray(mA, oA.transform, i * 16);
            Mat4.fromArray(mB, oB.transform, i * 16);
            if (!Mat4.areEqual(mA, mB, EPSILON)) return false;
        }
        return true;
    }
}