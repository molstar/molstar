/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Sebastian Bittrich <sebastian.bittrich@rcsb.org>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ShrakeRupleyContext, VdWLookup } from './common';
import { getElementIdx, isHydrogen } from '../../../../mol-model/structure/structure/unit/bonds/common';
import { isPolymer, isNucleic, MoleculeType, ElementSymbol } from '../../../../mol-model/structure/model/types';
import { VdwRadius } from '../../../../mol-model/structure/model/properties/atomic';
import { StructureElement, StructureProperties } from '../../../../mol-model/structure/structure';
import { getElementMoleculeType } from '../../../../mol-model/structure/util';

const l = StructureElement.Location.create(void 0);

export function assignRadiusForHeavyAtoms(ctx: ShrakeRupleyContext) {
    const { key } = StructureProperties.residue;
    const { type_symbol, label_atom_id, label_comp_id } = StructureProperties.atom;
    const { structure, atomRadiusType, serialResidueIndex } = ctx;

    let prevResidueIdx = 0;
    let residueIdx = 0;
    let serialResidueIdx = -1;

    l.structure = structure;
    for (let i = 0, m = 0, il = structure.units.length; i < il; ++i) {
        const unit = structure.units[i];
        const { elements } = unit;
        l.unit = unit;

        prevResidueIdx = -1;

        for (let j = 0, jl = elements.length; j < jl; ++j) {
            const eI = elements[j];
            const mj = m + j;

            l.element = eI;
            residueIdx = key(l);

            if (prevResidueIdx !== residueIdx) ++serialResidueIdx;
            prevResidueIdx = residueIdx;

            const element = type_symbol(l);
            const elementIdx = getElementIdx(element);

            // skip hydrogen atoms
            if (isHydrogen(elementIdx)) {
                atomRadiusType[mj] = VdWLookup[0];
                serialResidueIndex[mj] = -1;
                continue;
            }

            const moleculeType = getElementMoleculeType(unit, eI);
            // skip water and optionally non-polymer groups
            if (moleculeType === MoleculeType.Water || (!ctx.nonPolymer && !isPolymer(moleculeType))) {
                atomRadiusType[mj] = VdWLookup[0];
                serialResidueIndex[mj] = -1;
                continue;
            }

            const atomId = label_atom_id(l);
            const compId = label_comp_id(l);

            if (isNucleic(moleculeType)) {
                atomRadiusType[mj] = determineRadiusNucl(atomId, element, compId);
            } else if (moleculeType === MoleculeType.Protein) {
                atomRadiusType[mj] = determineRadiusAmino(atomId, element, compId);
            } else {
                atomRadiusType[mj] = handleNonStandardCase(element);
            }
            serialResidueIndex[mj] = serialResidueIdx;
        }
        m += elements.length;
    }
}

/**
 * Gets the van der Waals radius of the given atom following the values defined by Chothia (1976)
 * J.Mol.Biol.105,1-14. NOTE: the vdw values defined by the paper assume no Hydrogens and thus "inflates" slightly
 * the heavy atoms to account for Hydrogens.
 */
function determineRadiusAmino(atomId: string, element: ElementSymbol, compId: string): number {
    switch (element) {
        case 'O':
            return 5;
        case 'S':
            return 6;
        case 'N':
            return atomId === 'NZ' ? 4 : 3;
        case 'C':
            switch (atomId) {
                case 'C': case 'CE1': case'CE2': case 'CE3': case 'CH2': case 'CZ': case 'CZ2': case 'CZ3':
                    return 1;
                case 'CA': case 'CB': case 'CE': case 'CG1': case 'CG2':
                    return 2;
                default:
                    switch (compId) {
                        case 'PHE': case 'TRP': case 'TYR': case 'HIS': case 'ASP': case 'ASN':
                            return 1;
                        case 'PRO': case 'LYS': case 'ARG': case 'MET': case 'ILE': case 'LEU':
                            return 2;
                        case 'GLU': case 'GLN':
                            return atomId === 'CD' ? 1 : 2;
                    }
            }
    }
    return handleNonStandardCase(element);
}

function determineRadiusNucl(atomId: string, element: ElementSymbol, compId: string): number {
    switch (element) {
        case 'C': return 7;
        case 'N': return 8;
        case 'P': return 9;
        case 'O': return 5;
    }
    return handleNonStandardCase(element);
}

function handleNonStandardCase(element: ElementSymbol): number {
    const radius = VdwRadius(element);
    let index = VdWLookup.indexOf(radius);
    if (index === -1) {
        // add novel value to lookup array
        index = VdWLookup.length;
        VdWLookup[index] = radius;
    }
    return index;
}