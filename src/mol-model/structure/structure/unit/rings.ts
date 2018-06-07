/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import computeRings from './rings/compute'
import Unit from '../unit';

interface UnitRings {
    /** Each ring is specified as an array of indices in Unit.elements. */
    readonly all: ReadonlyArray<ReadonlyArray<number>>,
    readonly byFingerprint: Map<string, ReadonlyArray<number>>
}

namespace UnitRings {
    export function getRingFingerprint(unit: Unit.Atomic, ring: ArrayLike<number>) {
        const { elements } = unit;
        const { type_symbol } = unit.model.atomicHierarchy.atoms;

        const symbols: string[] = [];
        for (let i = 0, _i = ring.length; i < _i; i++) symbols[symbols.length] = type_symbol.value(elements[ring[i]]) as String as string;
        return getFingerprint(symbols);
    }

    export function create(unit: Unit.Atomic): UnitRings {
        const rings = computeRings(unit);
        const byFingerprint = new Map<string, number[]>();

        let idx = 0;
        for (const r of rings) {
            const fp = getRingFingerprint(unit, r);
            if (byFingerprint.has(fp)) byFingerprint.get(fp)!.push(idx);
            else byFingerprint.set(fp, [idx]);
            idx++;
        }

        return { all: rings, byFingerprint };
    }
}

export { UnitRings }

function getFingerprint(elements: string[]) {
    const len = elements.length;
    const reversed: string[] = new Array(len);

    for (let i = 0; i < len; i++) reversed[i] = elements[len - i - 1];

    const rotNormal = getMinimalRotation(elements);
    const rotReversed = getMinimalRotation(reversed);

    let isNormalSmaller = false;

    for (let i = 0; i < len; i++) {
        const u = elements[(i + rotNormal) % len], v = reversed[(i + rotReversed) % len];
        if (u !== v) {
            isNormalSmaller = u < v;
            break;
        }
    }

    if (isNormalSmaller) return buildFinderprint(elements, rotNormal);
    return buildFinderprint(reversed, rotReversed);
}

function getMinimalRotation(elements: string[]) {
    // adapted from http://en.wikipedia.org/wiki/Lexicographically_minimal_string_rotation

    const len = elements.length;
    const f = new Int32Array(len * 2);
    for (let i = 0; i < f.length; i++) f[i] = -1;

    let u = '', v = '', k = 0;

    for (let j = 1; j < f.length; j++) {
        let i = f[j - k - 1];
        while (i !== -1) {
            u = elements[j % len]; v = elements[(k + i + 1) % len];
            if (u === v) break;
            if (u < v) k = j - i - 1;
            i = f[i];
        }

        if (i === -1) {
            u = elements[j % len]; v = elements[(k + i + 1) % len];
            if (u !== v) {
                if (u < v) k = j;
                f[j - k] = -1;
            } else f[j - k] = i + 1;
        } else f[j - k] = i + 1;
    }

    return k;
}

function buildFinderprint(elements: string[], offset: number) {
    const len = elements.length;
    const ret: string[] = [];
    let i;
    for (i = 0; i < len - 1; i++) {
        ret.push(elements[(i + offset) % len]);
        ret.push('-');
    }
    ret.push(elements[(i + offset) % len]);
    return ret.join('');
}