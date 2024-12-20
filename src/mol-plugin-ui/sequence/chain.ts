/**
 * Copyright (c) 2019-2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Adam Midlik <midlik@gmail.com>
 */

import { Interval, OrderedSet } from '../../mol-data/int';
import { Loci } from '../../mol-model/loci';
import { Structure, StructureElement, StructureProperties } from '../../mol-model/structure';
import { ColorNames } from '../../mol-util/color/names';
import { SequenceWrapper, StructureUnit } from './wrapper';

export class ChainSequenceWrapper extends SequenceWrapper<StructureUnit> {
    private label: string;
    private unitIndices: Map<number, Interval<StructureElement.UnitIndex>>;
    private loci: StructureElement.Loci;

    residueLabel(seqIdx: number) {
        return this.label;
    }
    residueColor(seqIdx: number) {
        return ColorNames.black;
    }
    residueClass(seqIdx: number) {
        return 'msp-sequence-present';
    }

    override getSeqIndices(loci: Loci): OrderedSet {
        const { structure } = this.data;
        if (StructureElement.Loci.is(loci)) {
            if (!Structure.areRootsEquivalent(loci.structure, structure)) return Interval.Empty;
            loci = StructureElement.Loci.remap(loci, structure);

            for (const e of loci.elements) {
                const indices = this.unitIndices.get(e.unit.id);
                if (indices) {
                    if (OrderedSet.isSubset(indices, e.indices)) {
                        return Interval.ofSingleton(0);
                    }
                }
            }
        } else if (Structure.isLoci(loci)) {
            if (!Structure.areRootsEquivalent(loci.structure, structure)) return Interval.Empty;
            return Interval.ofSingleton(0);
        }
        return Interval.Empty;
    }

    getLoci(seqIdx: number) {
        return this.loci;
    }

    constructor(data: StructureUnit) {
        let residueCount = 0;
        let elementCount = 0;
        const counts: string[] = [];
        const l = StructureElement.Location.create(data.structure);

        const unitIndices = new Map<number, Interval<StructureElement.UnitIndex>>();
        const lociElements: StructureElement.Loci['elements'][0][] = [];

        for (let i = 0, il = data.units.length; i < il; ++i) {
            const unit = data.units[i];
            StructureElement.Location.set(l, data.structure, unit, unit.elements[0]);
            const entitySeq = unit.model.sequence.byEntityKey[StructureProperties.entity.key(l)];
            if (entitySeq) residueCount += entitySeq.sequence.length;
            elementCount += unit.elements.length;

            const indices = Interval.ofBounds(0, unit.elements.length);
            unitIndices.set(unit.id, indices);
            lociElements.push({ unit, indices });
        }

        if (residueCount > 0) counts.push(`${residueCount} residues`);
        counts.push(`${elementCount} elements`);

        const length = 1;

        super(data, length);

        this.label = `Whole Chain (${counts.join(', ')})`;
        this.unitIndices = unitIndices;
        this.loci = StructureElement.Loci(this.data.structure, lociElements);
    }
}
