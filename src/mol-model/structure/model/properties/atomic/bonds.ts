/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { AminoAcidNames, BaseNames } from '../../types';

/**
 * Map of intra component bond orders in aminoacids and nucleotides assuming standard IUPAC naming.
 * The key is constructed as `${compId}|${atomId1}|${atomId2}` with `atomId1 < atomId2`.
 */
const IntraBondOrderTable = new Map([
    ['HIS|CD2|CG', 2],
    ['HIS|CE1|ND1', 2],
    ['ARG|CZ|NH2', 2],
    ['PHE|CE1|CZ', 2],
    ['PHE|CD2|CE2', 2],
    ['PHE|CD1|CG', 2],
    ['TRP|CD1|CG', 2],
    ['TRP|CD2|CE2', 2],
    ['TRP|CE3|CZ3', 2],
    ['TRP|CH2|CZ2', 2],
    ['ASN|CG|OD1', 2],
    ['GLN|CD|OE1', 2],
    ['TYR|CD1|CG', 2],
    ['TYR|CD2|CE2', 2],
    ['TYR|CE1|CZ', 2],
    ['ASP|CG|OD1', 2],
    ['GLU|CD|OE1', 2],

    ['G|C8|N7', 2],
    ['G|C4|C5', 2],
    ['G|C2|N3', 2],
    ['G|C6|O6', 2],
    ['C|C4|N3', 2],
    ['C|C5|C6', 2],
    ['C|C2|O2', 2],
    ['A|C2|N3', 2],
    ['A|C6|N1', 2],
    ['A|C4|C5', 2],
    ['A|C8|N7', 2],
    ['U|C5|C6', 2],
    ['U|C2|O2', 2],
    ['U|C4|O4', 2],

    ['DG|C8|N7', 2],
    ['DG|C4|C5', 2],
    ['DG|C2|N3', 2],
    ['DG|C6|O6', 2],
    ['DC|C4|N3', 2],
    ['DC|C5|C6', 2],
    ['DC|C2|O2', 2],
    ['DA|C2|N3', 2],
    ['DA|C6|N1', 2],
    ['DA|C4|C5', 2],
    ['DA|C8|N7', 2],
    ['DT|C5|C6', 2],
    ['DT|C2|O2', 2],
    ['DT|C4|O4', 2]
]);

/**
 * Get order for bonds in aminoacids and nucleotides assuming standard IUPAC naming
 */
export function getIntraBondOrderFromTable (compId: string, atomId1: string, atomId2: string) {
    [ atomId1, atomId2 ] = atomId1 < atomId2 ? [ atomId1, atomId2 ] : [ atomId2, atomId1 ];
    if (AminoAcidNames.has(compId) && atomId1 === 'C' && atomId2 === 'O') return 2;
    if (BaseNames.has(compId) && atomId1 === 'OP1' && atomId2 === 'P') return 2;
    return IntraBondOrderTable.get(`${compId}|${atomId1}|${atomId2}`) || 1;
}

/**
 * Map of inter component bond orders assuming PDBx/mmCIF naming.
 * The key is constructed as `${compId1}|${compId2}|${atomId1}|${atomId2}` with `compId1 < compId2`.
 */
const InterBondOrderTable = new Map([
    ['LYS|NZ|RET|C15', 2] // Schiff base in Rhodopsin and Bacteriorhodopsin
]);

/**
 * Get order for bonds between component assuming PDBx/mmCIF naming.
 */
export function getInterBondOrderFromTable (compId1: string, atomId1: string, compId2: string, atomId2: string) {
    if (compId1 > compId2) {
        [ compId1, compId2 ] = [ compId2, compId1 ];
        [ atomId1, atomId2 ] = [ atomId2, atomId1 ];
    }
    return InterBondOrderTable.get(`${compId1}|${atomId1}|${compId2}|${atomId2}`) || 1;
}