/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Gianluca Tomasello <giagitom@gmail.com>
 */

import { Mat4 } from '../../linear-algebra/3d/mat4';
import { Vec3 } from '../../linear-algebra/3d/vec3';
import { Box3D } from './box3d';
import { Sphere3D } from './sphere3d';

interface Ray3D { origin: Vec3, direction: Vec3 }

function Ray3D() {
    return Ray3D.create(Vec3.create(0, 0, 0), Vec3.create(1, 0, 0));
}

namespace Ray3D {
    export function create(origin: Vec3, direction: Vec3): Ray3D { return { origin, direction }; }

    export function copy(out: Ray3D, r: Ray3D): Ray3D {
        Vec3.copy(out.origin, r.origin);
        Vec3.copy(out.direction, r.direction);
        return out;
    }

    export function clone(r: Ray3D): Ray3D {
        return copy(Ray3D(), r);
    }

    export function targetTo(out: Ray3D, ray: Ray3D, target: Vec3): Ray3D {
        Vec3.copy(out.origin, ray.origin);
        Vec3.normalize(out.direction, Vec3.sub(out.direction, target, ray.origin));
        return out;
    }

    /** Transform ray with a Mat4 */
    export function transform(out: Ray3D, ray: Ray3D, m: Mat4): Ray3D {
        Vec3.transformMat4(out.origin, ray.origin, m);
        Vec3.transformDirection(out.direction, ray.direction, m);
        return out;
    }

    //

    const tmpIR = Vec3();
    function _intersectSphere3D(ray: Ray3D, sphere: Sphere3D): number {
        const { center, radius } = sphere;
        const { origin, direction } = ray;

        const oc = Vec3.sub(tmpIR, origin, center);
        const a = Vec3.dot(direction, direction);
        const b = 2.0 * Vec3.dot(oc, direction);
        const c = Vec3.dot(oc, oc) - radius * radius;
        const discriminant = b * b - 4 * a * c;

        if (discriminant < 0) return -1; // no intersection

        const t = (-b - Math.sqrt(discriminant)) / (2.0 * a);
        if (t < 0) return -1; // behind the ray

        return t;
    }

    export function intersectSphere3D(out: Vec3, ray: Ray3D, sphere: Sphere3D): boolean {
        const t = _intersectSphere3D(ray, sphere);
        if (t < 0) return false;

        Vec3.scaleAndAdd(out, ray.origin, ray.direction, t);
        return true;
    }

    export function isIntersectingSphere3D(ray: Ray3D, sphere: Sphere3D): boolean {
        return _intersectSphere3D(ray, sphere) >= 0;
    }

    export function isInsideSphere3D(ray: Ray3D, sphere: Sphere3D): boolean {
        return Vec3.distance(ray.origin, sphere.center) < sphere.radius;
    }

    //

    function _intersectBox3D(ray: Ray3D, box: Box3D): number {
        const { origin, direction } = ray;
        const [minX, minY, minZ] = box.min;
        const [maxX, maxY, maxZ] = box.max;
        const [x, y, z] = origin;
        const invDirX = 1.0 / direction[0];
        const invDirY = 1.0 / direction[1];
        const invDirZ = 1.0 / direction[2];
        let tmin, tmax, tymin, tymax, tzmin, tzmax;
        if (invDirX >= 0) {
            tmin = (minX - x) * invDirX;
            tmax = (maxX - x) * invDirX;
        } else {
            tmin = (maxX - x) * invDirX;
            tmax = (minX - x) * invDirX;
        }
        if (invDirY >= 0) {
            tymin = (minY - y) * invDirY;
            tymax = (maxY - y) * invDirY;
        } else {
            tymin = (maxY - y) * invDirY;
            tymax = (minY - y) * invDirY;
        }
        if ((tmin > tymax) || (tymin > tmax)) return -1;
        if (tymin > tmin) tmin = tymin;
        if (tymax < tmax) tmax = tymax;
        if (invDirZ >= 0) {
            tzmin = (minZ - z) * invDirZ;
            tzmax = (maxZ - z) * invDirZ;
        } else {
            tzmin = (maxZ - z) * invDirZ;
            tzmax = (minZ - z) * invDirZ;
        }
        if ((tmin > tzmax) || (tzmin > tmax)) return -1;
        if (tzmin > tmin) tmin = tzmin;
        if (tzmax < tmax) tmax = tzmax;
        return tmin >= 0 ? tmin : -1;
    }

    export function intersectBox3D(out: Vec3, ray: Ray3D, box: Box3D): boolean {
        const t = _intersectBox3D(ray, box);
        if (t < 0) return false;

        Vec3.scaleAndAdd(out, ray.origin, ray.direction, t);
        return true;
    }

    export function isIntersectingBox3D(ray: Ray3D, box: Box3D): boolean {
        return _intersectBox3D(ray, box) >= 0;
    }
}

export { Ray3D };
