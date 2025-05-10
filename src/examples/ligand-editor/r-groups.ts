/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { JSONCifLigandGraph, JSONCifLigandGraphAtom } from '../../extensions/json-cif/ligand-graph';
import { molfileToJSONCif } from '../../extensions/json-cif/utils';
import { Mat4, Vec3 } from '../../mol-math/linear-algebra';

export type RGroupName = keyof typeof RGroups;

export async function attachRGroup(pGraph: JSONCifLigandGraph, rgroupName: RGroupName, pAtomOrId: number | JSONCifLigandGraphAtom) {
    const pAtom = pGraph.getAtom(pAtomOrId);
    if (pAtom?.row?.type_symbol !== 'H') {
        throw new Error('R-group attachment point must be a hydrogen atom.');
    }

    const { molfile, jsoncif: rgroupData } = await molfileToJSONCif(RGroups[rgroupName]);
    const attachIdx = molfile.attachmentPoints?.[0].atomIdx;

    if (typeof attachIdx !== 'number') {
        throw new Error('R-group attachment point not specified.');
    }

    // Compute and apply rGroup transformation
    const pBonds = pGraph.getBonds(pAtom);
    if (pBonds.length !== 1) {
        throw new Error('R-group attachment point must have exactly 1 bond.');
    }
    const pDir = pGraph.getBondDirection(pBonds[0]);
    const pPivot = pBonds[0].atom_2;
    Vec3.negate(pDir, pDir);
    Vec3.normalize(pDir, pDir);

    const rGraph = new JSONCifLigandGraph(rgroupData.dataBlocks[0]);
    const rAtom = rGraph.getAtomAtIndex(attachIdx - 1);

    if (rAtom.row?.type_symbol !== 'R#') {
        throw new Error('R-group attachment point is not a R# atom.');
    }

    const rCoords = rGraph.getAtomCoords(rAtom);
    const rBonds = rGraph.getBonds(rAtom);
    if (rBonds.length !== 1) {
        throw new Error('R-group R# atom must have exactly 1 bond.');
    }
    const rPivot = rGraph.getAtom(rBonds[0].atom_2);
    const rDir = rGraph.getBondDirection(rBonds[0]);
    Vec3.normalize(rDir, rDir);

    const rotation = Vec3.makeRotation(Mat4(), rDir, pDir);
    const translation = Mat4.fromTranslation(Mat4(), Vec3.sub(Vec3(), pGraph.getAtomCoords(pPivot), rCoords));

    const C = Mat4.fromTranslation(Mat4(), Vec3.negate(Vec3(), rCoords));
    const CT = Mat4.fromTranslation(Mat4(), rCoords);
    const T0 = Mat4.mul3(Mat4(), CT, rotation, C);
    const T = Mat4.mul(Mat4(), translation, T0);

    rGraph.transformCoords(T);

    // Merge the two graphs
    pGraph.removeAtom(pAtom);
    rGraph.removeAtom(rAtom);

    const newAtomMap = new Map<string, JSONCifLigandGraphAtom>();

    // Add atoms
    for (const a of rGraph.atoms) {
        const newAtom = pGraph.addAtom(a.row);
        newAtomMap.set(a.key, newAtom);
        if (a === rPivot) {
            pGraph.addOrUpdateBond(pPivot, newAtom, rBonds[0].props);
        }
    }

    // Add bonds
    for (const a of rGraph.atoms) {
        if (a === rAtom) continue;
        const bonds = rGraph.getBonds(a);
        const atom1 = newAtomMap.get(a.key)!;
        for (const b of bonds) {
            if (b.atom_2 === rAtom) continue;
            const atom2 = newAtomMap.get(b.atom_2.key)!;
            pGraph.addOrUpdateBond(atom1, atom2, b.props);
        }
    }

    return pGraph;
}

// Assumes the "attachment point (M  APO)" points to a hydrogen atom that gets removed
// when the R-group is attached.
const RGroups = {
    CH3: `CH3
  -OEChem-05072507373D

  5  4  0     0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.5541    0.7996    0.4965 R#   0  0  0  0  0  0  0  0  0  0  0  0
    0.6833   -0.8134   -0.2536 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7782   -0.3735    0.6692 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4593    0.3874   -0.9121 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  1  3  1  0  0  0  0
  1  4  1  0  0  0  0
  1  5  1  0  0  0  0
M  APO  1   2   1
M  END`
};