import { DnatcoTypes } from './types';
import { Segmentation } from '../../mol-data/int';
import { ElementIndex, ResidueIndex, StructureElement, StructureProperties } from '../../mol-model/structure';

const EmptyStepIndices = new Array<number>();

export namespace DnatcoUtil {
    export type Residue = Segmentation.Segment<ResidueIndex>;

    export function copyResidue(r?: Residue) {
        return r ? { index: r.index, start: r.start, end: r.end } : void 0;
    }

    export function getAtomIndex(loc: StructureElement.Location, residue: Residue, names: string[], altId: string): ElementIndex {
        for (let eI = residue.start; eI < residue.end; eI++) {
            loc.element = loc.unit.elements[eI];
            const elName = StructureProperties.atom.label_atom_id(loc);
            const elAltId = StructureProperties.atom.label_alt_id(loc);

            if (names.includes(elName) && (elAltId === altId || elAltId.length === 0))
                return loc.element;
        }

        return -1 as ElementIndex;
    }

    export function getStepIndices(data: DnatcoTypes.Steps, loc: StructureElement.Location, r: DnatcoUtil.Residue) {
        loc.element = loc.unit.elements[r.start];

        const modelIdx = StructureProperties.unit.model_num(loc) - 1;
        const chainId = StructureProperties.chain.auth_asym_id(loc);
        const seqId = StructureProperties.residue.auth_seq_id(loc);

        const chains = data.mapping[modelIdx];
        if (!chains) return EmptyStepIndices;
        const residues = chains.get(chainId);
        if (!residues) return EmptyStepIndices;
        return residues.get(seqId) ?? EmptyStepIndices;
    }
}
