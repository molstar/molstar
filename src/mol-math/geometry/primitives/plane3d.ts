/**
 * Copyright (c) 2022-2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 *
 * This code has been modified from https://github.com/mrdoob/three.js/,
 * copyright (c) 2010-2022 three.js authors. MIT License
 */

import { NumberArray } from '../../../mol-util/type-helpers';
import { Vec3 } from '../../linear-algebra/3d/vec3';
import { Sphere3D } from './sphere3d';

interface Plane3D { normal: Vec3, constant: number }

function Plane3D() {
    return Plane3D.create(Vec3.create(1, 0, 0), 0);
}

namespace Plane3D {
    export function create(normal: Vec3, constant: number): Plane3D { return { normal, constant }; }

    export function copy(out: Plane3D, p: Plane3D): Plane3D {
        Vec3.copy(out.normal, p.normal);
        out.constant = p.constant;
        return out;
    }

    export function clone(p: Plane3D): Plane3D {
        return copy(Plane3D(), p);
    }

    export function normalize(out: Plane3D, p: Plane3D): Plane3D {
        // Note: will lead to a divide by zero if the plane is invalid.
        const inverseNormalLength = 1.0 / Vec3.magnitude(p.normal);
        Vec3.scale(out.normal, p.normal, inverseNormalLength);
        out.constant = p.constant * inverseNormalLength;
        return out;
    }

    export function negate(out: Plane3D, p: Plane3D): Plane3D {
        Vec3.negate(out.normal, p.normal);
        out.constant = -p.constant;
        return out;
    }

    export function toArray<T extends NumberArray>(p: Plane3D, out: T, offset: number) {
        Vec3.toArray(p.normal, out, offset);
        out[offset + 3] = p.constant;
        return out;
    }

    export function fromArray(out: Plane3D, array: NumberArray, offset: number) {
        Vec3.fromArray(out.normal, array, offset);
        out.constant = array[offset + 3];
        return out;
    }

    export function fromNormalAndCoplanarPoint(out: Plane3D, normal: Vec3, point: Vec3) {
        Vec3.copy(out.normal, normal);
        out.constant = -Vec3.dot(out.normal, point);
        return out;
    }

    export function fromCoplanarPoints(out: Plane3D, a: Vec3, b: Vec3, c: Vec3) {
        const normal = Vec3.triangleNormal(Vec3(), a, b, c);
        fromNormalAndCoplanarPoint(out, normal, a);
        return out;
    }

    const unnormTmpV = Vec3();
    export function setUnnormalized(out: Plane3D, nx: number, ny: number, nz: number, constant: number) {
        Vec3.set(unnormTmpV, nx, ny, nz);
        const inverseNormalLength = 1.0 / Vec3.magnitude(unnormTmpV);
        Vec3.scale(out.normal, unnormTmpV, inverseNormalLength);
        out.constant = constant * inverseNormalLength;
        return out;
    }

    export function distanceToPoint(plane: Plane3D, point: Vec3) {
        return Vec3.dot(plane.normal, point) + plane.constant;
    }

    export function distanceToSpher3D(plane: Plane3D, sphere: Sphere3D) {
        return distanceToPoint(plane, sphere.center) - sphere.radius;
    }

    export function projectPoint(out: Vec3, plane: Plane3D, point: Vec3) {
        return Vec3.scaleAndAdd(out, out, plane.normal, -distanceToPoint(plane, point));
    }
}

export { Plane3D };
