/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Sebastian Bittrich <sebastian.bittrich@rcsb.org>
 */

import { ShrakeRupleyContext, VdWLookup } from './common';
import { getElementIdx, isHydrogen } from '../../../../mol-model/structure/structure/unit/links/common';
import { isPolymer, isNucleic, MoleculeType, ElementSymbol } from '../../../../mol-model/structure/model/types';
import { VdwRadius } from '../../../../mol-model/structure/model/properties/atomic';

export async function assignRadiusForHeavyAtoms(ctx: ShrakeRupleyContext) {
    const { updateChunk, atomRadius } = ctx;
    for (let i = 0; i < atomRadius.length; i += updateChunk) {
        computeRange(ctx, i, Math.min(i + updateChunk, atomRadius.length));
    }
}

function computeRange(ctx: ShrakeRupleyContext, begin: number, end: number) {
    const { structure } = ctx;
    const { model } = structure;
    const { atoms: atomInfo, derived, residues, residueAtomSegments } = model.atomicHierarchy;
    const { label_comp_id } = residues;
    const { moleculeType } = derived.residue;
    const { type_symbol, label_atom_id } = atomInfo;

    for (let aI = begin; aI < end; ++aI) {
        const rI = residueAtomSegments.index[aI];
        const element = type_symbol.value(aI);
        const elementIdx = getElementIdx(element);
        // skip hydrogen atoms
        if (isHydrogen(elementIdx)) {
            ctx.atomRadius[aI] = VdWLookup[0];
            continue;
        }

        const residueType = moleculeType[rI];
        // skip non-polymer groups
        if (!ctx.nonPolymer) {
            if (!isPolymer(residueType)) {
                ctx.atomRadius[aI] = VdWLookup[0];
                continue;
            }
        }

        const atomId = label_atom_id.value(aI);
        let compId = label_comp_id.value(rI);

        // handle modified residues
        const parentId = model.properties.modifiedResidues.parentId.get(compId);
        if (parentId !== void 0) compId = parentId;

        if (isNucleic(residueType)) {
             ctx.atomRadius[aI] = determineRadiusNucl(atomId, element, compId);
        } else if (residueType === MoleculeType.Protein) {
            ctx.atomRadius[aI] = determineRadiusAmino(atomId, element, compId);
        } else {
            ctx.atomRadius[aI] = handleNonStandardCase(element);
        }
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
        case 'C':
        return 7;
        case 'N':
        return 8;
        case 'P':
        return 9;
        case 'O':
        return 5;
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