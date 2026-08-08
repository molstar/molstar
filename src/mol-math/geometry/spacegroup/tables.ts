/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Mat4 } from '../../linear-algebra';
import { coordinateExpressionToOperator, transformOperators } from './common';
import { centringVectors, hasLongFormChangeOfBasis, operatorsFromHall } from './notation';
import { RawSpacegroupData } from './syminfo';

/**
 * The full spacegroup reference table: one entry per `syminfo.lib` setting
 * (540 total - all CCP4-numbered settings AND non-numbered axis/cell/origin
 * alternatives), built from `./syminfo`'s code-generated (compact tuple)
 * `RawSpacegroupData` (see `src/cli/syminfo` for how that's produced from
 * `data/sym/syminfo.lib`).
 */
export interface SpacegroupEntry {
    /** International Tables for Crystallography (ITA) spacegroup number (1-230), always present. */
    readonly itaNumber: number;
    /**
     * CCP4-style spacegroup number (mol*'s historical numbering, e.g. 1146,
     * 1003) - `0` for settings CCP4 doesn't assign a distinct number to
     * (most of `syminfo.lib`'s non-default axis/cell/origin descriptions).
     */
    readonly ccp4Number: number;
    /** Accepted Hermann-Mauguin name(s); index 0 is the canonical name. */
    readonly names: readonly string[];
    /**
     * Hall symbol for this exact setting. For entries with a non-zero
     * `ccp4Number` this is an independently-crafted symbol with no
     * change-of-basis suffix (so `operatorsFromHall` can parse it directly,
     * byte-identical to the pre-`syminfo.lib` table); for the rest it is
     * `syminfo.lib`'s own `symbol Hall` value verbatim, which may carry a
     * `(...)` comma-form change-of-basis suffix that `operatorsFromHall`
     * cannot parse - see `operatorsForEntry`, which handles both cases
     * transparently.
     */
    readonly hall: string;
    /**
     * ITA change-of-basis coordinate expression (`syminfo.lib`'s `basisop`)
     * relative to this ITA number's canonical (identity `'x,y,z'`) setting,
     * e.g. `'z,x,y'`, `'x-1/4,y-1/4,z-1/4'`, `'-y+z,x+z,-x+y+z'`.
     */
    readonly basisop: string;
}

export const SpacegroupData: readonly SpacegroupEntry[] = RawSpacegroupData.map(([itaNumber, ccp4Number, names, hall, basisop]) => ({ itaNumber, ccp4Number, names, hall, basisop }));

const EntriesByItaNumber: ReadonlyMap<number, readonly SpacegroupEntry[]> = (function () {
    const map = new Map<number, SpacegroupEntry[]>();
    for (const entry of SpacegroupData) {
        const list = map.get(entry.itaNumber);
        if (list) list.push(entry);
        else map.set(entry.itaNumber, [entry]);
    }
    return map;
}());

/** Maps a Hermann-Mauguin name (canonical or alias) to its `SpacegroupEntry`. */
export const SpacegroupEntryByName: ReadonlyMap<string, SpacegroupEntry> = (function () {
    const map = new Map<string, SpacegroupEntry>();
    // Two passes so a name shared between a CCP4-numbered setting and a
    // non-numbered alternate (this happens for the ~14 origin-choice-1/2
    // spacegroups, whose base Hermann-Mauguin name doesn't distinguish origin
    // choice) always resolves to the numbered one, exactly matching this
    // table's pre-`syminfo.lib` behavior.
    for (const entry of SpacegroupData) {
        if (entry.ccp4Number === 0) continue;
        for (const name of entry.names) map.set(name, entry);
    }
    for (const entry of SpacegroupData) {
        if (entry.ccp4Number !== 0) continue;
        for (const name of entry.names) if (!map.has(name)) map.set(name, entry);
    }
    return map;
}());

/** Maps a CCP4-style spacegroup number to its `SpacegroupEntry` (only the ~268 numbered settings). */
export const SpacegroupEntryByNumber: ReadonlyMap<number, SpacegroupEntry> = (function () {
    const map = new Map<number, SpacegroupEntry>();
    for (const entry of SpacegroupData) {
        if (entry.ccp4Number !== 0) map.set(entry.ccp4Number, entry);
    }
    return map;
}());

/** Resolves a Hermann-Mauguin name or CCP4-style number to its `SpacegroupEntry`, or `undefined` if unknown. */
export function findSpacegroupEntry(nameOrNumber: number | string): SpacegroupEntry | undefined {
    return typeof nameOrNumber === 'number'
        ? SpacegroupEntryByNumber.get(nameOrNumber)
        : SpacegroupEntryByName.get(nameOrNumber);
}

/** Maps a Hermann-Mauguin name (canonical or alias) to its spacegroup number. */
export const SpacegroupNameToNumberMap: { [name: string]: number } = (function () {
    const map: { [name: string]: number } = Object.create(null);
    for (const entry of SpacegroupData) {
        if (entry.ccp4Number === 0) continue;
        for (const name of entry.names) map[name] = entry.ccp4Number;
    }
    return map;
}());

/** Maps a CCP4-style spacegroup number to its canonical Hermann-Mauguin name. */
export const SpacegroupName: { [num: number]: string } = (function () {
    const names: { [num: number]: string } = Object.create(null);
    for (const [num, entry] of SpacegroupEntryByNumber) names[num] = entry.names[0];
    return names;
}());

/**
 * Returns the spacegroup number for a name or number - `-1` if unknown.
 * Numeric input is only resolved against the ~268 CCP4-style numbers (exact
 * match, same as before this table gained the extra `syminfo.lib` settings).
 * A NAME may resolve to a setting without its own CCP4-style number, in
 * which case its ITA number is returned instead (not unique across the
 * several such alternate settings that can share one ITA number, but a
 * strictly more useful result than "unknown").
 */
export function getSpacegroupNumber(nameOrNumber: number | string): number {
    const entry = findSpacegroupEntry(nameOrNumber);
    return entry ? (entry.ccp4Number !== 0 ? entry.ccp4Number : entry.itaNumber) : -1;
}

/**
 * Returns the Hall symbol for a CCP4-style spacegroup number (e.g. 1146), or
 * `undefined` if the number is unknown. Together with `operatorsFromHall`
 * this generates the full operator set from the symbol alone - but only for
 * the ~268 CCP4-numbered settings; use `operatorsForEntry` for any
 * `SpacegroupEntry`, numbered or not.
 */
export function getHallSymbol(spacegroupNumber: number): string | undefined {
    return SpacegroupEntryByNumber.get(spacegroupNumber)?.hall;
}

/**
 * `[first ITA number, point group order]` pairs, ascending - each pair covers
 * the ITA numbers up to (excluding) the next pair's first number.
 */
const PointGroupOrderRanges: readonly (readonly [number, number])[] = [
    [1, 1], [2, 2], [10, 4], [47, 8], // triclinic, monoclinic, orthorhombic
    [75, 4], [83, 8], [123, 16], // tetragonal
    [143, 3], [147, 6], [162, 12], // trigonal
    [168, 6], [175, 12], [191, 24], // hexagonal
    [195, 12], [200, 24], [221, 48], // cubic
];

const PointGroupOrderByItaNumber = (function () {
    const orders = new Uint8Array(231);
    for (let i = 0; i < PointGroupOrderRanges.length; i++) {
        const [start, order] = PointGroupOrderRanges[i];
        const end = i + 1 < PointGroupOrderRanges.length ? PointGroupOrderRanges[i + 1][0] : orders.length;
        for (let n = start; n < end; n++) orders[n] = order;
    }
    return orders;
}());

/**
 * Order of the point group (crystal class) an ITA spacegroup number belongs
 * to, i.e. the number of distinct rotation parts among its operators.
 */
export function pointGroupOrder(itaNumber: number): number {
    return PointGroupOrderByItaNumber[itaNumber] ?? 0;
}

/**
 * Number of symmetry operators of an entry (the multiplicity of the general
 * position) - its point group order times its lattice centering count. Both
 * are table lookups, so unlike `operatorsForEntry(entry).length` this builds
 * no matrices.
 */
export function orderForEntry(entry: SpacegroupEntry): number {
    // the Hall symbol's lattice letter, unlike the Hermann-Mauguin name's, is
    // 'P' for a rhombohedral-axes description of an R spacegroup
    const hall = entry.hall.trim();
    const latticeLetter = hall[0] === '-' ? hall.slice(1).trim()[0] : hall[0];
    return pointGroupOrder(entry.itaNumber) * centringVectors(latticeLetter).length;
}

const OperatorsByEntryCache = new Map<SpacegroupEntry, ReadonlyArray<Mat4>>();

/**
 * Resolves the full operator set for any `SpacegroupEntry`, numbered or not.
 * For entries with a directly-parseable Hall symbol (true for all
 * CCP4-numbered entries, since their Hall symbols are independently crafted
 * to avoid a change-of-basis suffix) this is just `operatorsFromHall(entry.hall)`.
 * Otherwise (most non-numbered `syminfo.lib` settings), operators are derived
 * by finding a sibling setting (same `itaNumber`) with a resolvable Hall
 * symbol, "undoing" that sibling's own `basisop` to get back to the common
 * canonical frame, then re-applying this entry's own `basisop`. Memoized per
 * entry (entries are stable singleton objects from `SpacegroupData`).
 */
export function operatorsForEntry(entry: SpacegroupEntry): ReadonlyArray<Mat4> {
    const cached = OperatorsByEntryCache.get(entry);
    if (cached) return cached;

    let operators: ReadonlyArray<Mat4>;
    if (!hasLongFormChangeOfBasis(entry.hall)) {
        operators = operatorsFromHall(entry.hall);
    } else {
        const siblings = EntriesByItaNumber.get(entry.itaNumber) ?? [];
        const anchor = siblings.find(e => e !== entry && !hasLongFormChangeOfBasis(e.hall));
        if (!anchor) throw new Error(`operatorsForEntry: no resolvable Hall symbol among ITA ${entry.itaNumber}'s settings`);
        const anchorOps = operatorsFromHall(anchor.hall) as Mat4[];
        const anchorP = coordinateExpressionToOperator(anchor.basisop);
        // Undo the anchor's own basisop to get back to the common frame that
        // every basisop in this ITA-number group is expressed relative to.
        const canonicalOps = transformOperators(anchorOps, anchorP);

        const P = coordinateExpressionToOperator(entry.basisop);
        const Pinv = Mat4.invert(Mat4(), P);
        if (!Pinv) throw new Error(`operatorsForEntry: non-invertible basisop '${entry.basisop}'`);
        operators = transformOperators(canonicalOps, Pinv);
    }
    OperatorsByEntryCache.set(entry, operators);
    return operators;
}
