/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Michal Malý <michal.maly@ibt.cas.cz>
 * @author Jiří Černý <jiri.cerny@ibt.cas.cz>
 */

import { ConfalPyramidsProvider } from './property';
import { ConfalPyramidsTypes as CPT } from './types';
import { DnatcoCommon as DC } from '../common';
import { OrderedSet, Segmentation } from '../../../mol-data/int';
import { Vec3 } from '../../../mol-math/linear-algebra';
import { Loci } from '../../../mol-model/loci';
import { ChainIndex, ElementIndex, ResidueIndex, Structure, StructureElement, StructureProperties, Unit } from '../../../mol-model/structure';

export namespace ConfalPyramidsUtil {
    const AllowedO3Names = [ 'O3\'', 'O3*' ];
    const AllowedOP1Names = [ 'OP1', 'O1P' ];

    type Residue = Segmentation.Segment<ResidueIndex>;

    export type AtomInfo = {
        pos: Vec3,
        index: ElementIndex,
        fakeAltId: string,
    };

    export type FirstResidueAtoms = {
        O3: AtomInfo,
    };

    export type SecondResidueAtoms = {
        OP1: AtomInfo,
        OP2: AtomInfo,
        O5: AtomInfo,
        P: AtomInfo,
    };

    type ResidueInfo = {
        PDB_model_num: number,
        asym_id: string,
        auth_asym_id: string,
        seq_id: number,
        auth_seq_id: number,
        comp_id: string,
        alt_id: string,
        ins_code: string,
    };

    export type Handler = (pyramid: CPT.Pyramid, first: FirstResidueAtoms, second: SecondResidueAtoms, firstLocIndex: number, secondLocIndex: number) => void;

    function residueInfoFromLocation(loc: StructureElement.Location): ResidueInfo {
        return {
            PDB_model_num: StructureProperties.unit.model_num(loc),
            asym_id: StructureProperties.chain.label_asym_id(loc),
            auth_asym_id: StructureProperties.chain.auth_asym_id(loc),
            seq_id: StructureProperties.residue.label_seq_id(loc),
            auth_seq_id: StructureProperties.residue.auth_seq_id(loc),
            comp_id: StructureProperties.atom.label_comp_id(loc),
            alt_id: StructureProperties.atom.label_alt_id(loc),
            ins_code: StructureProperties.residue.pdbx_PDB_ins_code(loc)
        };
    }

    export function hasMultipleModels(unit: Unit.Atomic): boolean {
        const prop = ConfalPyramidsProvider.get(unit.model).value;
        if (prop === undefined || prop.data === undefined) throw new Error('No custom properties data');
        return prop.data.hasMultipleModels;
    }

    function getPossibleAltIdsIndices(eIFirst: ElementIndex, eILast: ElementIndex, structure: Structure, unit: Unit.Atomic): string[] {
        const loc = StructureElement.Location.create(structure, unit, -1 as ElementIndex);

        const uIFirst = OrderedSet.indexOf(unit.elements, eIFirst);
        const uILast = OrderedSet.indexOf(unit.elements, eILast);

        const possibleAltIds: string[] = [];
        for (let uI = uIFirst; uI <= uILast; uI++) {
            loc.element = unit.elements[uI];
            const altId = StructureProperties.atom.label_alt_id(loc);
            if (altId !== '' && !possibleAltIds.includes(altId)) possibleAltIds.push(altId);
        }

        return possibleAltIds;
    }

    function getPossibleAltIdsResidue(residue: Residue, structure: Structure, unit: Unit.Atomic): string[] {
        return getPossibleAltIdsIndices(unit.elements[residue.start], unit.elements[residue.end], structure, unit);
    }

    class Utility {
        protected getAtomIndices(names: string[], residue: Residue, allowedAltIds: string[]): ElementIndex[] {
            let rI = residue.start;
            const rILast = residue.end;
            const indices: ElementIndex[] = [];

            for (; rI <= rILast; rI++) {
                const eI = this.unit.elements[rI];
                const loc = StructureElement.Location.create(this.structure, this.unit, eI);
                const thisName = StructureProperties.atom.label_atom_id(loc);
                if (allowedAltIds.length > 0) {
                    const thisAltId = StructureProperties.atom.label_alt_id(loc);
                    if (!allowedAltIds.includes(thisAltId)) continue;
                }
                if (names.includes(thisName)) indices.push(eI);
            }

            if (indices.length === 0) {
                let namesStr = '';
                for (const n of names)
                    namesStr += `${n} `;

                throw new Error(`Element [${namesStr}] not found on residue ${residue.index}`);
            }

            return indices;
        }

        private getAtomPositions(indices: ElementIndex[]): Vec3[] {
            const pos = this.unit.conformation.invariantPosition;
            const positions: Vec3[] = [];

            for (const eI of indices) {
                const v = Vec3.zero();
                pos(eI, v);
                positions.push(v);
            }

            return positions;
        }

        protected getPyramidByName(name: string): { pyramid: CPT.Pyramid | undefined, index: number } {
            const index = this.data.names.get(name);
            if (index === undefined) return { pyramid: undefined, index: -1 };

            return { pyramid: this.data.pyramids[index], index };
        }

        protected processFirstResidue(residue: Residue, possibleAltIds: string[]) {
            const indO3 = this.getAtomIndices(AllowedO3Names, residue, []);
            const posO3 = this.getAtomPositions(indO3);

            const altPos: FirstResidueAtoms[] = [
                { O3: { pos: posO3[0], index: indO3[0], fakeAltId: '' } }
            ];

            for (let i = 1; i < indO3.length; i++) {
                altPos.push({ O3: { pos: posO3[i], index: indO3[i], fakeAltId: '' } });
            }

            if (altPos.length === 1 && possibleAltIds.length > 1) {
                /* We have some alternate positions on the residue but O3 does not have any - fake them */
                altPos[0].O3.fakeAltId = possibleAltIds[0];

                for (let i = 1; i < possibleAltIds.length; i++)
                    altPos.push({ O3: { pos: posO3[0], index: indO3[0], fakeAltId: possibleAltIds[i] } });
            }

            return altPos;
        }

        protected processSecondResidue(residue: Residue, possibleAltIds: string[]) {
            const indOP1 = this.getAtomIndices(AllowedOP1Names, residue, []);
            const indOP2 = this.getAtomIndices(['OP2'], residue, []);
            const indO5 = this.getAtomIndices(['O5\'', 'O5*'], residue, []);
            const indP = this.getAtomIndices(['P'], residue, []);

            const posOP1 = this.getAtomPositions(indOP1);
            const posOP2 = this.getAtomPositions(indOP2);
            const posO5 = this.getAtomPositions(indO5);
            const posP = this.getAtomPositions(indP);

            const infoOP1: AtomInfo[] = [];
            /* We use OP1 as "pivotal" atom. There is no specific reason
             * to pick OP1, it is as good a choice as any other atom
             */
            if (indOP1.length === 1 && possibleAltIds.length > 1) {
                /* No altIds on OP1, fake them */
                for (const altId of possibleAltIds)
                    infoOP1.push({ pos: posOP1[0], index: indOP1[0], fakeAltId: altId });
            } else {
                for (let i = 0; i < indOP1.length; i++)
                    infoOP1.push({ pos: posOP1[i], index: indOP1[i], fakeAltId: '' });
            }

            const mkInfo = (i: number, indices: ElementIndex[], positions: Vec3[], altId: string) => {
                if (i >= indices.length) {
                    const last = indices.length - 1;
                    return { pos: positions[last], index: indices[last], fakeAltId: altId };
                }

                return { pos: positions[i], index: indices[i], fakeAltId: altId };
            };

            const altPos: SecondResidueAtoms[] = [];
            for (let i = 0; i < infoOP1.length; i++) {
                const altId = infoOP1[i].fakeAltId;

                const OP2 = mkInfo(i, indOP2, posOP2, altId);
                const O5 = mkInfo(i, indO5, posO5, altId);
                const P = mkInfo(i, indP, posP, altId);

                altPos.push({ OP1: infoOP1[i], OP2, O5, P });
            }

            return altPos;
        }

        protected stepToName(entry_id: string, modelNum: number, locFirst: StructureElement.Location, locSecond: StructureElement.Location, fakeAltId_1: string, fakeAltId_2: string) {
            const first = residueInfoFromLocation(locFirst);
            const second = residueInfoFromLocation(locSecond);
            const model_id = this.hasMultipleModels ? `-m${modelNum}` : '';
            const alt_id_1 =  fakeAltId_1 !== '' ? `.${fakeAltId_1}` : (first.alt_id.length ? `.${first.alt_id}` : '');
            const alt_id_2 = fakeAltId_2 !== '' ? `.${fakeAltId_2}` : (second.alt_id.length ? `.${second.alt_id}` : '');
            const ins_code_1 = first.ins_code.length ? `.${first.ins_code}` : '';
            const ins_code_2 = second.ins_code.length ? `.${second.ins_code}` : '';

            return `${entry_id}${model_id}_${first.auth_asym_id}_${first.comp_id}${alt_id_1}_${first.auth_seq_id}${ins_code_1}_${second.comp_id}${alt_id_2}_${second.auth_seq_id}${ins_code_2}`;
        }

        constructor(protected structure: Structure, protected unit: Unit.Atomic) {
            const prop = ConfalPyramidsProvider.get(unit.model).value;
            if (prop === undefined || prop.data === undefined) throw new Error('No custom properties data');

            this.data = prop.data;
            this.hasMultipleModels = hasMultipleModels(unit);

            this.entryId = unit.model.entryId.toLowerCase();
            this.modelNum = unit.model.modelNum;
        }

        protected readonly data: CPT.PyramidsData
        protected readonly hasMultipleModels: boolean;
        protected readonly entryId: string;
        protected readonly modelNum: number;
    }

    export class UnitWalker extends Utility {
        private handleStep(firstAtoms: FirstResidueAtoms[], secondAtoms: SecondResidueAtoms[]) {
            const modelNum = this.hasMultipleModels ? this.modelNum : -1;
            let ok = false;

            const firstLoc = StructureElement.Location.create(this.structure, this.unit, -1 as ElementIndex);
            const secondLoc = StructureElement.Location.create(this.structure, this.unit, -1 as ElementIndex);
            for (let i = 0; i < firstAtoms.length; i++) {
                const first = firstAtoms[i];
                for (let j = 0; j < secondAtoms.length; j++) {
                    const second = secondAtoms[j];
                    firstLoc.element = first.O3.index;
                    secondLoc.element = second.OP1.index;

                    const name = this.stepToName(this.entryId, modelNum, firstLoc, secondLoc, first.O3.fakeAltId, second.OP1.fakeAltId);
                    const { pyramid, index } = this.getPyramidByName(name);
                    if (pyramid !== undefined) {
                        const setLoc = (loc: CPT.Location, eI: ElementIndex) => {
                            loc.element.structure = this.structure;
                            loc.element.unit = this.unit;
                            loc.element.element = eI;
                        };

                        const locIndex = index * 2;
                        setLoc(this.data.locations[locIndex], firstLoc.element);
                        setLoc(this.data.locations[locIndex + 1], secondLoc.element);
                        this.handler(pyramid, first, second, locIndex, locIndex + 1);
                        ok = true;
                    }
                }
            }

            if (!ok) throw new Error('Bogus step');
        }

        private step(residue: Residue): { firstAtoms: FirstResidueAtoms[], secondAtoms: SecondResidueAtoms[] } {
            const firstPossibleAltIds = getPossibleAltIdsResidue(residue, this.structure, this.unit);
            const firstAtoms = this.processFirstResidue(residue, firstPossibleAltIds);

            residue = this.residueIt.move();

            const secondPossibleAltIds = getPossibleAltIdsResidue(residue, this.structure, this.unit);
            const secondAtoms = this.processSecondResidue(residue, secondPossibleAltIds);

            return { firstAtoms, secondAtoms };
        }

        walk() {
            while (this.chainIt.hasNext) {
                this.residueIt.setSegment(this.chainIt.move());

                let residue = this.residueIt.move();
                while (this.residueIt.hasNext) {
                    try {
                        const { firstAtoms, secondAtoms } = this.step(residue);

                        this.handleStep(firstAtoms, secondAtoms);
                    } catch (error) {
                        /* Skip and move along */
                        residue = this.residueIt.move();
                    }
                }
            }
        }

        constructor(structure: Structure, unit: Unit.Atomic, private handler: Handler) {
            super(structure, unit);

            this.chainIt = Segmentation.transientSegments(unit.model.atomicHierarchy.chainAtomSegments, unit.elements);
            this.residueIt = Segmentation.transientSegments(unit.model.atomicHierarchy.residueAtomSegments, unit.elements);
        }

        private chainIt: Segmentation.SegmentIterator<ChainIndex>;
        private residueIt: Segmentation.SegmentIterator<ResidueIndex>;
    }

    class PyramidGetter extends Utility {
        private residuesToAtoms(firstResidue: Residue, firstAltIds: string[], secondResidue: Residue, secondAltIds: string[]) {
            const firstAtoms = this.getAtomIndices(AllowedO3Names, firstResidue, firstAltIds);
            const secondAtoms = this.getAtomIndices(AllowedOP1Names, secondResidue, secondAltIds);

            return { firstAtoms, secondAtoms };
        }

        get() {
            const indices = this.loci.elements[0].indices;
            const NI = OrderedSet.size(indices);
            const unitElements = this.unit.elements;
            const { label_seq_id } = this.unit.model.atomicHierarchy.residues;
            const { label_alt_id } = this.unit.model.atomicHierarchy.atoms;
            const { index: residueIndex } = this.unit.model.atomicHierarchy.residueAtomSegments;
            const uIFirst = OrderedSet.getAt(indices, 0);
            const uILast = OrderedSet.getAt(indices, NI - 1);
            const resnoFirst = label_seq_id.value(residueIndex[unitElements[uIFirst]]);
            const resnoLast = label_seq_id.value(residueIndex[unitElements[uILast]]);
            const possibleFirstAltIds: string[] = [];
            const possibleSecondAltIds: string[] = [];

            if (resnoFirst + 1 !== resnoLast) return void 0; // Check that we have two consecutive residues

            let uISecond = uIFirst;
            let i = 1;
            for (; i < NI; i++) {
                const uI = OrderedSet.getAt(indices, i);
                const resno = label_seq_id.value(residueIndex[unitElements[uI]]);
                const altId = label_alt_id.value(unitElements[uI]);
                if (resno === resnoLast) {
                    uISecond = uI;
                    break;
                } else {
                    if (altId !== '' && !possibleFirstAltIds.includes(altId)) possibleFirstAltIds.push(altId);
                }
            }
            if (i === NI)
                return void 0; // No valid second residue

            // Fill possible second residue alt ids array
            for (i; i < NI; i++) {
                const uI = OrderedSet.getAt(indices, i);
                const altId = label_alt_id.value(unitElements[uI]);
                if (altId !== '' && !possibleSecondAltIds.includes(altId)) possibleSecondAltIds.push(altId);
            }

            if (possibleFirstAltIds.length > 1 || possibleSecondAltIds.length > 1) return void 0; // We got selection with more than one conformation

            const firstResidue = { index: residueIndex[unitElements[uIFirst]], start: uIFirst, end: uISecond - 1 };
            const secondResidue = { index: residueIndex[unitElements[uISecond]], start: uISecond, end: uILast };

            try {
                const { firstAtoms, secondAtoms } = this.residuesToAtoms(firstResidue, possibleFirstAltIds, secondResidue, possibleSecondAltIds);

                const modelNum = this.hasMultipleModels ? this.modelNum : -1;
                const firstLoc = StructureElement.Location.create(this.structure, this.unit, -1 as ElementIndex);
                const secondLoc = StructureElement.Location.create(this.structure, this.unit, -1 as ElementIndex);
                for (let i = 0; i < firstAtoms.length; i++) {
                    const first = firstAtoms[i];
                    for (let j = 0; j < secondAtoms.length; j++) {
                        const second = secondAtoms[j];
                        firstLoc.element = first;
                        secondLoc.element = second;

                        const name = this.stepToName(this.entryId, modelNum, firstLoc, secondLoc, '', '');
                        const { pyramid } = this.getPyramidByName(name);
                        if (pyramid !== undefined) return pyramid;
                    }
                }

                return void 0;
            } catch (e) {
                return void 0;
            }
        }

        constructor(unit: Unit.Atomic, private loci: StructureElement.Loci) {
            super(loci.structure, unit);
        }
    }

    export function lociToPyramid(loci: Loci) {
        if (loci.kind !== 'element-loci') return void 0;
        if (loci.elements.length < 1) return void 0;
        if (!DC.isApplicable(loci.structure.model)) return void 0;
        if (!Unit.isAtomic(loci.elements[0].unit)) return void 0;

        return (new PyramidGetter(loci.elements[0].unit, loci)).get();
    }
}
