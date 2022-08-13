/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Michal Malý <michal.maly@ibt.cas.cz>
 * @author Jiří Černý <jiri.cerny@ibt.cas.cz>
 */

import { ConfalPyramidsProvider } from './property';
import { ConfalPyramidsTypes as CPT } from './types';
import { Segmentation } from '../../../mol-data/int';
import { ChainIndex, ElementIndex, ResidueIndex, Structure, StructureElement, StructureProperties, Unit } from '../../../mol-model/structure';

type Residue = Segmentation.Segment<ResidueIndex>;

export type Pyramid = {
    O3: ElementIndex,
    P: ElementIndex,
    OP1: ElementIndex,
    OP2: ElementIndex,
    O5: ElementIndex,
    confalScore: number,
    stepIdx: number,
};

const EmptyStepIndices = new Array<number>();

function copyResidue(r?: Residue) {
    return r ? { index: r.index, start: r.start, end: r.end } : void 0;
}

function getAtomIndex(loc: StructureElement.Location, residue: Residue, names: string[], altId: string): ElementIndex {
    for (let eI = residue.start; eI < residue.end; eI++) {
        loc.element = loc.unit.elements[eI];
        const elName = StructureProperties.atom.label_atom_id(loc);
        const elAltId = StructureProperties.atom.label_alt_id(loc);

        if (names.includes(elName) && (elAltId === altId || elAltId.length === 0))
            return loc.element;
    }

    return -1 as ElementIndex;
}

function getPyramid(loc: StructureElement.Location, one: Residue, two: Residue, altIdOne: string, altIdTwo: string, confalScore: number, stepIdx: number): Pyramid {
    const O3 = getAtomIndex(loc, one, ['O3\'', 'O3*'], altIdOne);
    const P = getAtomIndex(loc, two, ['P'], altIdTwo);
    const OP1 = getAtomIndex(loc, two, ['OP1'], altIdTwo);
    const OP2 = getAtomIndex(loc, two, ['OP2'], altIdTwo);
    const O5 = getAtomIndex(loc, two, ['O5\'', 'O5*'], altIdTwo);

    return { O3, P, OP1, OP2, O5, confalScore, stepIdx };
}

export class ConfalPyramidsIterator {
    private chainIt: Segmentation.SegmentIterator<ChainIndex>;
    private residueIt: Segmentation.SegmentIterator<ResidueIndex>;
    private residueOne?: Residue;
    private residueTwo: Residue;
    private data?: CPT.Steps;
    private loc: StructureElement.Location;

    private getStepIndices(r: Residue) {
        this.loc.element = this.loc.unit.elements[r.start];

        const modelIdx = StructureProperties.unit.model_num(this.loc) - 1;
        const chainId = StructureProperties.chain.auth_asym_id(this.loc);
        const seqId = StructureProperties.residue.auth_seq_id(this.loc);

        const chains = this.data!.mapping[modelIdx];
        if (!chains) return EmptyStepIndices;
        const residues = chains.get(chainId);
        if (!residues) return EmptyStepIndices;
        return residues.get(seqId) ?? EmptyStepIndices;
    }

    private moveStep() {
        this.residueOne = copyResidue(this.residueTwo);
        this.residueTwo = copyResidue(this.residueIt.move())!;

        return this.toPyramids(this.residueOne!, this.residueTwo);
    }

    private toPyramids(one: Residue, two: Residue) {
        const indices = this.getStepIndices(one);

        const points = [];
        for (const idx of indices) {
            const step = this.data!.steps[idx];
            points.push(getPyramid(this.loc, one, two, step.label_alt_id_1, step.label_alt_id_2, step.confal_score, idx));
        }

        return points;
    }

    constructor(structure: Structure, unit: Unit) {
        this.chainIt = Segmentation.transientSegments(unit.model.atomicHierarchy.chainAtomSegments, unit.elements);
        this.residueIt = Segmentation.transientSegments(unit.model.atomicHierarchy.residueAtomSegments, unit.elements);

        const prop = ConfalPyramidsProvider.get(unit.model).value;
        this.data = prop?.data;

        if (this.chainIt.hasNext) {
            this.residueIt.setSegment(this.chainIt.move());
            if (this.residueIt.hasNext)
                this.residueTwo = this.residueIt.move();
        }

        this.loc = StructureElement.Location.create(structure, unit, -1 as ElementIndex);
    }

    get hasNext() {
        if (!this.data)
            return false;
        return this.residueIt.hasNext
            ? true
            : this.chainIt.hasNext;
    }

    move() {
        if (this.residueIt.hasNext) {
            return this.moveStep();
        } else {
            this.residueIt.setSegment(this.chainIt.move());
            return this.moveStep();
        }
    }
}
