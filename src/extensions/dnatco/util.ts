/**
 * Copyright (c) 2018-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Michal Malý <michal.maly@ibt.cas.cz>
 * @author Jiří Černý <jiri.cerny@ibt.cas.cz>
 */

import { DnatcoTypes } from './types';
import { OrderedSet, Segmentation } from '../../mol-data/int';
import { EmptyLoci } from '../../mol-model/loci';
import { ElementIndex, ResidueIndex, Structure, StructureElement, StructureProperties, Unit } from '../../mol-model/structure';

const EmptyStepIndices = new Array<number>();

export namespace DnatcoUtil {
    export type Residue = Segmentation.Segment<ResidueIndex>;

    export function copyResidue(r?: Residue) {
        return r ? { index: r.index, start: r.start, end: r.end } : void 0;
    }

    export function getAtomIndex(loc: StructureElement.Location, residue: Residue, names: string[], altId: string, insCode: string): ElementIndex {
        for (let eI = residue.start; eI < residue.end; eI++) {
            loc.element = loc.unit.elements[eI];
            const elName = StructureProperties.atom.label_atom_id(loc);
            const elAltId = StructureProperties.atom.label_alt_id(loc);
            const elInsCode = StructureProperties.residue.pdbx_PDB_ins_code(loc);

            if (names.includes(elName) && (elAltId === altId || elAltId.length === 0) && (elInsCode === insCode))
                return loc.element;
        }

        return -1 as ElementIndex;
    }

    export function getStepIndices(data: DnatcoTypes.Steps, loc: StructureElement.Location, r: DnatcoUtil.Residue) {
        loc.element = loc.unit.elements[r.start];

        const modelIdx = StructureProperties.unit.model_num(loc) - 1;
        const chainId = StructureProperties.chain.auth_asym_id(loc);
        const seqId = StructureProperties.residue.auth_seq_id(loc);
        const insCode = StructureProperties.residue.pdbx_PDB_ins_code(loc);

        const chains = data.mapping[modelIdx];
        if (!chains) return EmptyStepIndices;
        const residues = chains.get(chainId);
        if (!residues) return EmptyStepIndices;
        const indices = residues.get(seqId);
        if (!indices) return EmptyStepIndices;

        return insCode !== '' ? indices.filter(idx => data.steps[idx].PDB_ins_code_1 === insCode) : indices;
    }

    export function residueAltIds(structure: Structure, unit: Unit, residue: Residue) {
        const altIds = new Array<string>();
        const loc = StructureElement.Location.create(structure, unit);
        for (let eI = residue.start; eI < residue.end; eI++) {
            loc.element = OrderedSet.getAt(unit.elements, eI);
            const altId = StructureProperties.atom.label_alt_id(loc);
            if (altId !== '' && !altIds.includes(altId))
                altIds.push(altId);
        }

        return altIds;
    }

    const _loc = StructureElement.Location.create();
    export function residueToLoci(asymId: string, seqId: number, altId: string | undefined, insCode: string, loci: StructureElement.Loci, source: 'label' | 'auth') {
        _loc.structure = loci.structure;
        for (const e of loci.elements) {
            _loc.unit = e.unit;

            const getAsymId = source === 'label' ? StructureProperties.chain.label_asym_id : StructureProperties.chain.auth_asym_id;
            const getSeqId = source === 'label' ? StructureProperties.residue.label_seq_id : StructureProperties.residue.auth_seq_id;

            // Walk the entire unit and look for the requested residue
            const chainIt = Segmentation.transientSegments(e.unit.model.atomicHierarchy.chainAtomSegments, e.unit.elements);
            const residueIt = Segmentation.transientSegments(e.unit.model.atomicHierarchy.residueAtomSegments, e.unit.elements);

            const elemIndex = (idx: number) => OrderedSet.getAt(e.unit.elements, idx);
            while (chainIt.hasNext) {
                const chain = chainIt.move();
                _loc.element = elemIndex(chain.start);
                const _asymId = getAsymId(_loc);
                if (_asymId !== asymId)
                    continue; // Wrong chain, skip it

                residueIt.setSegment(chain);
                while (residueIt.hasNext) {
                    const residue = residueIt.move();
                    _loc.element = elemIndex(residue.start);

                    const _seqId = getSeqId(_loc);
                    if (_seqId === seqId) {
                        const _insCode = StructureProperties.residue.pdbx_PDB_ins_code(_loc);
                        if (_insCode !== insCode)
                            continue;
                        if (altId) {
                            const _altIds = residueAltIds(loci.structure, e.unit, residue);
                            if (!_altIds.includes(altId))
                                continue;
                        }

                        const start = residue.start as StructureElement.UnitIndex;
                        const end = residue.end as StructureElement.UnitIndex;
                        return StructureElement.Loci(
                            loci.structure,
                            [{ unit: e.unit, indices: OrderedSet.ofBounds(start, end) }]
                        );
                    }
                }
            }
        }

        return EmptyLoci;
    }
}
