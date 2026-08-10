/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Vec3, Mat4 } from '../../linear-algebra';
import { findSpacegroupEntry, getSpacegroupNumber, operatorsForEntry, orderForEntry } from './tables';
import { SymmetryOperator } from '../../geometry/symmetry-operator';
import { Cell } from './cell';

interface SpacegroupCell extends Cell {
    /** Hermann-Mauguin spacegroup name (uniquely identifies the spacegroup) */
    readonly name: string
}

interface Spacegroup {
    /** Hermann-Mauguin spacegroup name */
    readonly name: string,
    /**
     * Spacegroup number - the CCP4-style number (e.g. 1146) when the
     * resolved setting has one, else its plain ITA number as a fallback
     * (not unique across the several non-CCP4-numbered settings that can
     * share one ITA number - see `SpacegroupEntry` in `./tables`).
     */
    readonly num: number,
    readonly cell: SpacegroupCell,
    /**
     * ITA change-of-basis coordinate expression relative to the canonical
     * (identity `'x,y,z'`) description, in the same xyz-style `basisop`
     * notation as `SpacegroupEntry.basisop` (e.g. `'x,y,z'`, `'z,x,y'`,
     * `'x-1/4,y-1/4,z-1/4'`, `'y,x,-z'`, `'x+y,-x+y,z'` - see
     * `./common`/`./notation`). Used by `createWithSetting` callers such as
     * MOL2 CRYSIN as well as by `Spacegroup.create`.
     */
    readonly basisop: string,
    readonly operators: ReadonlyArray<Mat4>
}

namespace SpacegroupCell {
    /** Create a 'P 1' with cellsize [1, 1, 1] */
    export const Zero: SpacegroupCell = create('P 1', Vec3.create(1, 1, 1), Vec3.create(Math.PI / 2, Math.PI / 2, Math.PI / 2));

    /** True if the cell is a `SpacegroupCell` */
    export function is(cell: Cell): cell is SpacegroupCell {
        return typeof (cell as SpacegroupCell).name === 'string';
    }

    /** True if 'P 1' with cellsize [1, 1, 1] */
    export function isZero(cell?: SpacegroupCell) {
        if (!cell) return true;
        return cell.name === 'P 1' && cell.size[0] === 1 && cell.size[1] === 1 && cell.size[2] === 1;
    }

    /** Returns Zero cell if the spacegroup does not exist */
    export function create(nameOrNumber: number | string, size: Vec3, anglesInRadians: Vec3): SpacegroupCell {
        const entry = findSpacegroupEntry(nameOrNumber);
        if (!entry) {
            console.warn(`Unknown spacegroup '${nameOrNumber}', returning a 'P 1' with cellsize [1, 1, 1]`);
            return Zero;
        }

        const cell = Cell.create(size, anglesInRadians, orderForEntry(entry));

        return { ...cell, name: entry.names[0] };
    }

    export function getLabel(cell: SpacegroupCell) {
        return Cell.getLabel(cell, `Unit Cell <b>${cell.name}</b> #${getSpacegroupNumber(cell.name)}`);
    }
}

namespace Spacegroup {
    /** P1 with [1, 1, 1] cell */
    export const ZeroP1 = create(SpacegroupCell.Zero);

    export function create(cell: SpacegroupCell): Spacegroup {
        const entry = findSpacegroupEntry(cell.name);
        if (!entry) throw new Error(`missing spacegroup entry for '${cell.name}'`);
        const num = entry.ccp4Number !== 0 ? entry.ccp4Number : entry.itaNumber;
        const operators = operatorsForEntry(entry);
        return { name: cell.name, num, cell, basisop: entry.basisop, operators };
    }

    /**
     * Assembles a `Spacegroup` from an already-computed operator list for a
     * non-default axis setting/description, given its xyz-style `basisop`
     * change-of-basis coordinate expression (e.g. `'y,x,-z'` or
     * `'x+y,-x+y,z'` - see the `Spacegroup.basisop` field). `num` is derived
     * from `cell.name` exactly as in `create`. This is intentionally
     * agnostic of any particular file format's own setting-numbering scheme
     * - callers that need to resolve such a scheme (e.g. MOL2 CRYSIN) into
     * `name` + `basisop` + `operators` should do so themselves and pass the
     * result here.
     */
    export function createWithSetting(cell: SpacegroupCell, name: string, basisop: string, operators: ReadonlyArray<Mat4>): Spacegroup {
        const entry = findSpacegroupEntry(cell.name);
        const num = entry ? (entry.ccp4Number !== 0 ? entry.ccp4Number : entry.itaNumber) : -1;
        return { name, num, cell: { ...cell, order: operators.length }, basisop, operators };
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

    export function getSymmetryOperator(spacegroup: Spacegroup, spgrOp: number, i: number, j: number, k: number): SymmetryOperator {
        const operator = setOperatorMatrix(spacegroup, spgrOp, i, j, k, Mat4.zero());
        return SymmetryOperator.create(SymmetryOperator.getSymmetryOperatorName(spgrOp, i, j, k), operator, { hkl: Vec3.create(i, j, k), spgrOp });
    }

    const _translationRef = Vec3();
    const _translationRefSymop = Vec3();
    const _translationRefOffset = Vec3();
    const _translationSymop = Vec3();

    /**
     * Get Symmetry operator for transformation around the given
     * reference point `ref` in fractional coordinates
     */
    export function getSymmetryOperatorRef(spacegroup: Spacegroup, spgrOp: number, i: number, j: number, k: number, ref: Vec3) {

        const operator = Mat4.zero();

        Vec3.set(_ijkVec, i, j, k);
        Vec3.floor(_translationRef, ref);

        Mat4.copy(operator, spacegroup.operators[spgrOp]);

        Vec3.floor(_translationRefSymop, Vec3.transformMat4(_translationRefSymop, ref, operator));

        Mat4.getTranslation(_translationSymop, operator);
        Vec3.sub(_translationSymop, _translationSymop, _translationRefSymop);
        Vec3.add(_translationSymop, _translationSymop, _translationRef);
        Vec3.add(_translationSymop, _translationSymop, _ijkVec);

        Mat4.setTranslation(operator, _translationSymop);
        Mat4.mul(operator, spacegroup.cell.fromFractional, operator);
        Mat4.mul(operator, operator, spacegroup.cell.toFractional);

        Vec3.sub(_translationRefOffset, _translationRefSymop, _translationRef);

        const _i = i - _translationRefOffset[0];
        const _j = j - _translationRefOffset[1];
        const _k = k - _translationRefOffset[2];

        // const operator = setOperatorMatrixRef(spacegroup, spgrOp, i, j, k, ref, Mat4.zero());
        return SymmetryOperator.create(SymmetryOperator.getSymmetryOperatorName(spgrOp, _i, _j, _k), operator, { hkl: Vec3.create(_i, _j, _k), spgrOp });
    }

    export function getOperatorXyz(op: Mat4) {
        return [
            formatElement(getRotation(op[0], op[4], op[8]), getShift(op[12])),
            formatElement(getRotation(op[1], op[5], op[9]), getShift(op[13])),
            formatElement(getRotation(op[2], op[6], op[10]), getShift(op[14]))
        ].join(',');
    }

    function getRotation(x: number, y: number, z: number) {
        const r: string[] = [];
        if (x > 0) r.push('+X');
        else if (x < 0) r.push('-X');
        if (y > 0) r.push('+Y');
        else if (y < 0) r.push('-Y');
        if (z > 0) r.push('+Z');
        else if (z < 0) r.push('-Z');

        if (r.length === 1) {
            return r[0].charAt(0) === '+' ? r[0].substr(1) : r[0];
        }
        if (r.length === 2) {
            const s0 = r[0].charAt(0);
            const s1 = r[1].charAt(0);
            if (s0 === '+') return `${r[0].substr(1)}${r[1]}`;
            if (s1 === '+') return `${r[1].substr(1)}${r[0]}`;
        }
        throw new Error(`unknown rotation '${r}', ${x} ${y} ${z}`);
    }

    function getShift(s: number) {
        switch (s) {
            case 1 / 2: return '1/2';
            case 1 / 4: return '1/4';
            case 3 / 4: return '3/4';
            case 1 / 3: return '1/3';
            case 2 / 3: return '2/3';
            case 1 / 6: return '1/6';
            case 5 / 6: return '5/6';
        }
        return '';
    }

    function formatElement(rotation: string, shift: string) {
        if (shift === '') return rotation;
        if (rotation.length > 2) return `${rotation}+${shift}`;
        return rotation.charAt(0) === '-' ? `${shift}${rotation}` : `${shift}+${rotation}`;
    }
}

export { Spacegroup, SpacegroupCell };