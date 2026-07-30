/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Mat4, ReadonlyMat4 } from '../../linear-algebra/3d/mat4';
import { Vec3 } from '../../linear-algebra/3d/vec3';
import { closeGroup, seitz } from './common';

const v3 = Vec3.create;

const HallLatticeVectors: { [k: string]: Vec3[] } = {
    P: [v3(0, 0, 0)],
    A: [v3(0, 0, 0), v3(0, 0.5, 0.5)],
    B: [v3(0, 0, 0), v3(0.5, 0, 0.5)],
    C: [v3(0, 0, 0), v3(0.5, 0.5, 0)],
    I: [v3(0, 0, 0), v3(0.5, 0.5, 0.5)],
    R: [v3(0, 0, 0), v3(2 / 3, 1 / 3, 1 / 3), v3(1 / 3, 2 / 3, 2 / 3)],
    F: [v3(0, 0, 0), v3(0, 0.5, 0.5), v3(0.5, 0, 0.5), v3(0.5, 0.5, 0)],
};

/** Lattice-centering translation vectors for a lattice letter (P/A/B/C/I/R/F). */
export function centringVectors(letter: string): Vec3[] {
    const v = HallLatticeVectors[letter.toUpperCase()];
    if (!v) throw new Error(`unknown lattice symbol '${letter}'`);
    return v;
}

/**
 * Rotation matrices (as rows) about the c axis for the Hall n-fold orders, plus
 * the diagonal axes
 *
 * - ' (2-fold along [1 -1 0])
 * - " (2-fold along [1 1 0]) and
 * - * (3-fold along [1 1 1]).
 */
function hallRotation(axis: number | string): Vec3[] {
    switch (axis) {
        case 1: return [v3(1, 0, 0), v3(0, 1, 0), v3(0, 0, 1)];
        case 2: return [v3(-1, 0, 0), v3(0, -1, 0), v3(0, 0, 1)];
        case 3: return [v3(0, -1, 0), v3(1, -1, 0), v3(0, 0, 1)];
        case 4: return [v3(0, -1, 0), v3(1, 0, 0), v3(0, 0, 1)];
        case 6: return [v3(1, -1, 0), v3(1, 0, 0), v3(0, 0, 1)];
        case '\'': return [v3(0, -1, 0), v3(-1, 0, 0), v3(0, 0, -1)];
        case '"': return [v3(0, 1, 0), v3(1, 0, 0), v3(0, 0, -1)];
        case '*': return [v3(0, 0, 1), v3(1, 0, 0), v3(0, 1, 0)];
        default: throw new Error(`incorrect Hall axis '${axis}'`);
    }
}

const HallTranslations: { [k: string]: Vec3 } = {
    a: v3(0.5, 0, 0), b: v3(0, 0.5, 0), c: v3(0, 0, 0.5),
    n: v3(0.5, 0.5, 0.5), u: v3(0.25, 0, 0), v: v3(0, 0.25, 0),
    w: v3(0, 0, 0.25), d: v3(0.25, 0.25, 0.25),
};

function negateRotation(r: Vec3[]): Vec3[] {
    return [
        v3(-r[0][0], -r[0][1], -r[0][2]),
        v3(-r[1][0], -r[1][1], -r[1][2]),
        v3(-r[2][0], -r[2][1], -r[2][2])
    ];
}

/**
 * Reorients a c-axis rotation onto the x or y axis by simultaneously permuting
 * its rows and columns.
 */
function reorient(r: Vec3[], i: number, j: number, k: number): Vec3[] {
    return [
        v3(r[i][i], r[i][j], r[i][k]),
        v3(r[j][i], r[j][j], r[j][k]),
        v3(r[k][i], r[k][j], r[k][k]),
    ];
}

function hallGenerator(sym: string, position: number, previousOrder: number): { op: Mat4, order: number } {
    let p = 0;
    const negated = sym[p] === '-';
    if (negated) p++;
    const orderChar = sym[p];
    if (orderChar < '1' || orderChar === '5' || orderChar > '6') throw new Error(`invalid Hall rotation order in '${sym}'`);
    const order = orderChar.charCodeAt(0) - 48;
    p++;

    let screw = 0;
    let principalAxis = '';
    let diagonalAxis = '';
    const tran = v3(0, 0, 0);
    for (; p < sym.length; p++) {
        const ch = sym[p];
        if (ch >= '1' && ch <= '5') {
            screw = ch.charCodeAt(0) - 48;
        } else if (ch === '\'' || ch === '"' || ch === '*') {
            diagonalAxis = ch;
        } else if (ch === 'x' || ch === 'y' || ch === 'z') {
            principalAxis = ch;
        } else {
            const t = HallTranslations[ch];
            if (!t) throw new Error(`unknown Hall translation '${ch}' in '${sym}'`);
            tran[0] += t[0]; tran[1] += t[1]; tran[2] += t[2];
        }
    }

    if (!principalAxis && !diagonalAxis) {
        if (position === 1) principalAxis = 'z';
        else if (position === 2 && order === 2) {
            if (previousOrder === 2 || previousOrder === 4) principalAxis = 'x';
            else if (previousOrder === 3 || previousOrder === 6) diagonalAxis = '\'';
        } else if (position === 3 && order === 3) {
            diagonalAxis = '*';
        } else if (order !== 1) {
            throw new Error(`missing axis in Hall symbol part '${sym}'`);
        }
    }

    let rot = hallRotation(diagonalAxis ? diagonalAxis : order);
    if (negated) rot = negateRotation(rot);
    if (principalAxis === 'x') rot = reorient(rot, 2, 0, 1);
    else if (principalAxis === 'y') rot = reorient(rot, 1, 2, 0);

    if (screw) {
        const axisIndex = principalAxis ? principalAxis.charCodeAt(0) - 'x'.charCodeAt(0) : 2;
        tran[axisIndex] += screw / order;
    }

    return { op: seitz(rot, tran), order };
}

const Identity3: Vec3[] = [v3(1, 0, 0), v3(0, 1, 0), v3(0, 0, 1)];

function parseHallOriginShift(inner: string): Mat4 {
    if (inner.indexOf(',') >= 0) throw new Error(`unsupported long-form Hall change-of-basis '(${inner})'`);
    const v = inner.trim().split(/\s+/).map(x => parseInt(x, 10));
    const shift = v.map(n => (((n % 12) + 12) % 12) / 12);
    return seitz(Identity3, shift);
}

const OperatorsFromHallCache = new Map<string, ReadonlyArray<ReadonlyMat4>>();

/**
 * Generates the full set of symmetry operators (as fractional `Mat4`) for a
 * Hall space-group symbol, e.g. `"-P 2ac 2ab"` or `"P 31 2 (0 0 4)"`. The
 * operators are derived purely from the symbol - no per-space-group operator
 * table is consulted. Results are memoized per symbol; callers must treat the
 * returned array (and its matrices) as read-only.
 */
export function operatorsFromHall(hall: string): ReadonlyArray<ReadonlyMat4> {
    const cached = OperatorsFromHallCache.get(hall);
    if (cached) return cached;

    let s = hall.trim();
    const centrosymmetric = s[0] === '-';
    if (centrosymmetric) s = s.slice(1).trim();

    const latticeLetter = s[0];
    const centering = centringVectors(latticeLetter);
    let rest = s.slice(1).trim();

    let originShift: string | undefined;
    const open = rest.indexOf('(');
    if (open >= 0) {
        const close = rest.indexOf(')', open);
        if (close < 0) throw new Error(`missing ')' in Hall symbol '${hall}'`);
        originShift = rest.slice(open + 1, close).trim();
        rest = rest.slice(0, open).trim();
    }

    const generators: Mat4[] = [];
    for (const cv of centering) {
        if (cv[0] !== 0 || cv[1] !== 0 || cv[2] !== 0) generators.push(seitz(Identity3, cv));
    }

    const parts = rest.length ? rest.split(/[\s_]+/) : [];
    let previousOrder = 0;
    let position = 0;
    for (const part of parts) {
        position++;
        if (part === '1') continue;
        const { op, order } = hallGenerator(part, position, previousOrder);
        generators.push(op);
        previousOrder = order;
    }

    if (centrosymmetric) generators.push(seitz([v3(-1, 0, 0), v3(0, -1, 0), v3(0, 0, -1)], v3(0, 0, 0)));

    if (originShift !== undefined) {
        const cob = parseHallOriginShift(originShift);
        const cobInv = Mat4.invert(Mat4(), cob);
        if (!cobInv) throw new Error('non-invertible Hall change-of-basis');
        for (let i = 0; i < generators.length; i++) {
            const tmp = Mat4.mul(Mat4(), cob, generators[i]);
            generators[i] = Mat4.mul(Mat4(), tmp, cobInv);
        }
    }

    const operators = closeGroup(generators);
    OperatorsFromHallCache.set(hall, operators);
    return operators;
}

/**
 * `true` if a `symbol Hall` string carries a long-form (comma-separated)
 * change-of-basis suffix, e.g. `'P 2y (z,x,y)'` - which `operatorsFromHall`
 * cannot parse (only the short numeric-shift form, e.g.
 * `'P 31 2 (0 0 4)'`, is supported).
 */
export function hasLongFormChangeOfBasis(hall: string): boolean {
    const open = hall.indexOf('(');
    if (open < 0) return false;
    const close = hall.indexOf(')', open);
    return close > open && hall.slice(open + 1, close).includes(',');
}
