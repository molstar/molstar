/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Vec3, Mat4 } from '../../linear-algebra'
import { SpacegroupName, TransformData, GroupData, getSpacegroupIndex, OperatorData, SpacegroupNames } from './tables'
import { SymmetryOperator } from '../../geometry';

interface SpacegroupCell {
    /** Zero based spacegroup number */
    readonly index: number,
    readonly size: Vec3,
    readonly volume: number,
    readonly anglesInRadians: Vec3,
    /** Transfrom cartesian -> fractional coordinates within the cell */
    readonly toFractional: Mat4,
    /** Transfrom fractional coordinates within the cell -> cartesian */
    readonly fromFractional: Mat4
}

interface Spacegroup {
    readonly name: string,
    readonly cell: SpacegroupCell,
    readonly operators: ReadonlyArray<Mat4>
}

namespace SpacegroupCell {
    /** Create a 'P 1' with cellsize [1, 1, 1] */
    export const Zero: SpacegroupCell = create('P 1', Vec3.create(1, 1, 1), Vec3.create(Math.PI / 2, Math.PI / 2, Math.PI / 2));

    /** True if 'P 1' with cellsize [1, 1, 1] */
    export function isZero(cell: SpacegroupCell) {
        return cell.index === 0 && cell.size[0] === 1 && cell.size[1] === 1 && cell.size[1] === 1;
    }

    /** Returns Zero cell if the spacegroup does not exist */
    export function create(nameOrNumber: number | string | SpacegroupName, size: Vec3, anglesInRadians: Vec3): SpacegroupCell {
        const index = getSpacegroupIndex(nameOrNumber);
        if (index < 0) {
            console.warn(`Unknown spacegroup '${nameOrNumber}', returning a 'P 1' with cellsize [1, 1, 1]`)
            return Zero;
        }

        const volume = size[0] * size[1] * size[2]

        const alpha = anglesInRadians[0];
        const beta = anglesInRadians[1];
        const gamma = anglesInRadians[2];

        const xScale = size[0], yScale = size[1], zScale = size[2];

        const z1 = Math.cos(beta);
        const z2 = (Math.cos(alpha) - Math.cos(beta) * Math.cos(gamma)) / Math.sin(gamma);
        const z3 = Math.sqrt(1.0 - z1 * z1 - z2 * z2);

        const x = [xScale, 0.0, 0.0];
        const y = [Math.cos(gamma) * yScale, Math.sin(gamma) * yScale, 0.0];
        const z = [z1 * zScale, z2 * zScale, z3 * zScale];

        const fromFractional = Mat4.ofRows([
            [x[0], y[0], z[0], 0],
            [0, y[1], z[1], 0],
            [0, 0, z[2], 0],
            [0, 0, 0, 1.0]
        ]);
        const toFractional = Mat4.invert(Mat4.zero(), fromFractional)!;

        return { index, size, volume, anglesInRadians, toFractional, fromFractional };
    }
}


namespace Spacegroup {
    /** P1 with [1, 1, 1] cell */
    export const ZeroP1 = create(SpacegroupCell.Zero);

    export function create(cell: SpacegroupCell): Spacegroup {
        const operators = GroupData[cell.index].map(i => getOperatorMatrix(OperatorData[i]));
        return { name: SpacegroupNames[cell.index], cell, operators };
    }

    const _ijkVec = Vec3();
    const _tempMat = Mat4();
    export function setOperatorMatrix(spacegroup: Spacegroup, index: number, i: number, j: number, k: number, target: Mat4) {
        Vec3.set(_ijkVec, i, j, k);

        Mat4.fromTranslation(_tempMat, _ijkVec);
        return Mat4.mul(
            target,
            Mat4.mul(
                target,
                Mat4.mul(target, spacegroup.cell.fromFractional, _tempMat),
                spacegroup.operators[index]
            ),
            spacegroup.cell.toFractional
        );
    }

    export function getSymmetryOperator(spacegroup: Spacegroup, index: number, i: number, j: number, k: number): SymmetryOperator {
        const operator = setOperatorMatrix(spacegroup, index, i, j, k, Mat4.zero());
        return SymmetryOperator.create(`${index + 1}_${5 + i}${5 + j}${5 + k}`, operator, { id: '', operList: [] }, '', Vec3.create(i, j, k), index);
    }

    const _translationRef = Vec3()
    const _translationRefSymop = Vec3()
    const _translationSymop = Vec3()
    export function setOperatorMatrixRef(spacegroup: Spacegroup, index: number, i: number, j: number, k: number, ref: Vec3, target: Mat4) {
        Vec3.set(_ijkVec, i, j, k);
        Vec3.floor(_translationRef, ref)

        Mat4.copy(target, spacegroup.operators[index])

        Vec3.floor(_translationRefSymop, Vec3.transformMat4(_translationRefSymop, ref, target))

        Mat4.getTranslation(_translationSymop, target)
        Vec3.sub(_translationSymop, _translationSymop, _translationRefSymop)
        Vec3.add(_translationSymop, _translationSymop, _translationRef)
        Vec3.add(_translationSymop, _translationSymop, _ijkVec)

        Mat4.setTranslation(target, _translationSymop)
        Mat4.mul(target, spacegroup.cell.fromFractional, target)
        Mat4.mul(target, target, spacegroup.cell.toFractional)
        return target
    }

    /**
     * Get Symmetry operator for transformation around the given
     * reference point `ref` in fractional coordinates
     */
    export function getSymmetryOperatorRef(spacegroup: Spacegroup, index: number, i: number, j: number, k: number, ref: Vec3) {
        const operator = setOperatorMatrixRef(spacegroup, index, i, j, k, ref, Mat4.zero());
        return SymmetryOperator.create(`${index + 1}_${5 + i}${5 + j}${5 + k}`, operator, { id: '', operList: [] }, '', Vec3.create(i, j, k), index);
    }

    function getOperatorMatrix(ids: number[]) {
        const r1 = TransformData[ids[0]];
        const r2 = TransformData[ids[1]];
        const r3 = TransformData[ids[2]];
        return Mat4.ofRows([r1, r2, r3, [0, 0, 0, 1]]);
    }
}

export { Spacegroup, SpacegroupCell }