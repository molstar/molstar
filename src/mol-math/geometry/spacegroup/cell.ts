/**
 * Copyright (c) 2019-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Mat4, Vec3 } from '../../linear-algebra';
import { radToDeg } from '../../misc';

export { Cell };

interface Cell {
    readonly size: Vec3
    readonly anglesInRadians: Vec3
    readonly volume: number,
    /** Number of symmetry operators of the spacegroup the cell belongs to, 1 for a plain box */
    readonly order: number,
    /** Transfrom cartesian -> fractional coordinates within the cell */
    readonly toFractional: Mat4,
    /** Transfrom fractional coordinates within the cell -> cartesian */
    readonly fromFractional: Mat4
}

function Cell() {
    return Cell.empty();
}

namespace Cell {
    export function create(size: Vec3, anglesInRadians: Vec3, order = 1): Cell {
        const alpha = anglesInRadians[0];
        const beta = anglesInRadians[1];
        const gamma = anglesInRadians[2];

        const ca = Math.cos(alpha), cb = Math.cos(beta), cg = Math.cos(gamma);
        const volume = size[0] * size[1] * size[2] * Math.sqrt(1.0 - ca * ca - cb * cb - cg * cg + 2.0 * ca * cb * cg);

        const xScale = size[0], yScale = size[1], zScale = size[2];

        const z1 = cb;
        const z2 = (ca - cb * cg) / Math.sin(gamma);
        const z3 = Math.sqrt(1.0 - z1 * z1 - z2 * z2);

        const x = [xScale, 0.0, 0.0];
        const y = [cg * yScale, Math.sin(gamma) * yScale, 0.0];
        const z = [z1 * zScale, z2 * zScale, z3 * zScale];

        const fromFractional = Mat4.ofRows([
            [x[0], y[0], z[0], 0],
            [0, y[1], z[1], 0],
            [0, 0, z[2], 0],
            [0, 0, 0, 1.0]
        ]);
        const toFractional = Mat4.invert(Mat4.zero(), fromFractional);

        return { size, anglesInRadians, volume, order, toFractional, fromFractional };
    }

    export function empty(): Cell {
        return create(Vec3(), Vec3());
    }

    export function fromBasis(x: Vec3, y: Vec3, z: Vec3) {
        const a = Vec3.magnitude(x);
        const b = Vec3.magnitude(y);
        const c = Vec3.magnitude(z);

        const alpha = Math.acos(Vec3.dot(y, z) / (b * c));
        const beta = Math.acos(Vec3.dot(x, z) / (a * c));
        const gamma = Math.acos(Vec3.dot(x, y) / (a * b));

        if (a <= 0 || b <= 0 || c <= 0 || alpha >= Math.PI || beta >= Math.PI || gamma >= Math.PI) {
            return empty();
        } else {
            return create(Vec3.create(a, b, c), Vec3.create(alpha, beta, gamma));
        }
    }

    /** Human readable dimensions and angles, prefixed by `title` */
    export function getLabel(cell: Cell, title = 'Unit Cell') {
        const { size, anglesInRadians } = cell;
        const a = size[0].toFixed(2);
        const b = size[1].toFixed(2);
        const c = size[2].toFixed(2);
        const alpha = radToDeg(anglesInRadians[0]).toFixed(2);
        const beta = radToDeg(anglesInRadians[1]).toFixed(2);
        const gamma = radToDeg(anglesInRadians[2]).toFixed(2);
        return [
            title,
            `${a}\u00D7${b}\u00D7${c} \u212B`,
            `\u03b1=${alpha}\u00B0 \u03b2=${beta}\u00B0 \u03b3=${gamma}\u00B0`
        ].join(' | ');
    }
}