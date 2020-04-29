/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { StructureSelection, StructureQuery, Structure, Queries, StructureProperties as SP, StructureElement, Unit } from '../../mol-model/structure';
import { SequenceWrapper, StructureUnit } from './wrapper';
import { OrderedSet, Interval, SortedArray } from '../../mol-data/int';
import { Loci } from '../../mol-model/loci';
import { Sequence } from '../../mol-model/sequence';
import { MissingResidues } from '../../mol-model/structure/model/properties/common';
import { ColorNames } from '../../mol-util/color/names';
import { MarkerAction, applyMarkerAction, applyMarkerActionAtPosition } from '../../mol-util/marker-action';

export class PolymerSequenceWrapper extends SequenceWrapper<StructureUnit> {
    private readonly unitMap: Map<number, Unit>
    private readonly sequence: Sequence
    private readonly missing: MissingResidues
    private readonly observed: OrderedSet // sequences indices

    private readonly modelNum: number
    private readonly asymId: string

    private seqId(seqIdx: number) {
        return this.sequence.seqId.value(seqIdx);
    }

    residueLabel(seqIdx: number) {
        return this.sequence.label.value(seqIdx);
    }

    residueColor(seqIdx: number) {
        return this.missing.has(this.modelNum, this.asymId, this.seqId(seqIdx))
            ? ColorNames.grey
            : ColorNames.black;
    }
    residueClass(seqIdx: number) {
        return this.missing.has(this.modelNum, this.asymId, this.seqId(seqIdx))
            ? 'msp-sequence-missing'
            : 'msp-sequence-present';
    }

    mark(loci: Loci, action: MarkerAction): boolean {
        let changed = false;
        const { structure } = this.data;
        const index = (seqId: number) => this.sequence.index(seqId);
        if (StructureElement.Loci.is(loci)) {
            if (!Structure.areRootsEquivalent(loci.structure, structure)) return false;
            loci = StructureElement.Loci.remap(loci, structure);

            for (const e of loci.elements) {
                if (!this.unitMap.has(e.unit.id)) continue;

                if (Unit.isAtomic(e.unit)) {
                    changed = applyMarkerAtomic(e, action, this.markerArray, index) || changed;
                } else {
                    changed = applyMarkerCoarse(e, action, this.markerArray, index) || changed;
                }
            }
        } else if (Structure.isLoci(loci)) {
            if (!Structure.areRootsEquivalent(loci.structure, structure)) return false;
            if (applyMarkerAction(this.markerArray, this.observed, action)) changed = true;
        }
        return changed;
    }

    getLoci(seqIdx: number) {
        const query = createResidueQuery(this.data.units[0].chainGroupId, this.data.units[0].conformation.operator.name, this.seqId(seqIdx));
        return StructureSelection.toLociWithSourceUnits(StructureQuery.run(query, this.data.structure));
    }

    constructor(data: StructureUnit) {
        const l = StructureElement.Location.create(data.structure, data.units[0], data.units[0].elements[0]);
        const entitySeq = data.units[0].model.sequence.byEntityKey[SP.entity.key(l)];

        const length = entitySeq.sequence.length;
        const markerArray = new Uint8Array(length);

        super(data, markerArray, length);

        this.unitMap = new Map();
        for (const unit of data.units) this.unitMap.set(unit.id, unit);

        this.sequence = entitySeq.sequence;
        this.missing = data.units[0].model.properties.missingResidues;

        this.modelNum = data.units[0].model.modelNum;
        this.asymId = Unit.isAtomic(data.units[0]) ? SP.chain.label_asym_id(l) : SP.coarse.asym_id(l);

        const missing: number[] = [];
        for (let i = 0; i < length; ++i) {
            if (this.missing.has(this.modelNum, this.asymId, this.seqId(i))) missing.push(i);
        }
        this.observed = OrderedSet.subtract(Interval.ofBounds(0, length),  SortedArray.ofSortedArray(missing));
    }
}

function createResidueQuery(chainGroupId: number, operatorName: string, label_seq_id: number) {
    return Queries.generators.atoms({
        unitTest: ctx => {
            return (
                SP.unit.chainGroupId(ctx.element) === chainGroupId &&
                SP.unit.operator_name(ctx.element) === operatorName
            );
        },
        residueTest: ctx => {
            if (ctx.element.unit.kind === Unit.Kind.Atomic) {
                return SP.residue.label_seq_id(ctx.element) === label_seq_id;
            } else {
                return (
                    SP.coarse.seq_id_begin(ctx.element) <= label_seq_id &&
                    SP.coarse.seq_id_end(ctx.element) >= label_seq_id
                );
            }
        }
    });
}

function applyMarkerAtomic(e: StructureElement.Loci.Element, action: MarkerAction, markerArray: Uint8Array, index: (seqId: number) => number) {
    const { model, elements } = e.unit;
    const { index: residueIndex } = model.atomicHierarchy.residueAtomSegments;
    const { label_seq_id } = model.atomicHierarchy.residues;

    let changed = false;
    OrderedSet.forEachSegment(e.indices, i => residueIndex[elements[i]], rI => {
        const seqId = label_seq_id.value(rI);
        changed = applyMarkerActionAtPosition(markerArray, index(seqId), action) || changed;
    });
    return changed;
}

function applyMarkerCoarse(e: StructureElement.Loci.Element, action: MarkerAction, markerArray: Uint8Array, index: (seqId: number) => number) {
    const { model, elements } = e.unit;
    const begin = Unit.isSpheres(e.unit) ? model.coarseHierarchy.spheres.seq_id_begin : model.coarseHierarchy.gaussians.seq_id_begin;
    const end = Unit.isSpheres(e.unit) ? model.coarseHierarchy.spheres.seq_id_end : model.coarseHierarchy.gaussians.seq_id_end;

    let changed = false;
    OrderedSet.forEach(e.indices, i => {
        const eI = elements[i];
        for (let s = index(begin.value(eI)), e = index(end.value(eI)); s <= e; s++) {
            changed = applyMarkerActionAtPosition(markerArray, s, action) || changed;
        }
    });
    return changed;
}