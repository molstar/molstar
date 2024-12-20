/**
 * Copyright (c) 2019-2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Adam Midlik <midlik@gmail.com>
 */

import { Interval, OrderedSet } from '../../mol-data/int';
import { Loci } from '../../mol-model/loci';
import { Structure, StructureElement } from '../../mol-model/structure';
import { ColorNames } from '../../mol-util/color/names';
import { SequenceWrapper, StructureUnit } from './wrapper';

export class ElementSequenceWrapper extends SequenceWrapper<StructureUnit> {
    private unitIndices: Map<number, Interval<StructureElement.UnitIndex>>;

    residueLabel(seqIdx: number) {
        return 'X';
    }
    residueColor(seqIdx: number) {
        return ColorNames.black;
    }
    residueClass(seqIdx: number) {
        return 'msp-sequence-present';
    }

    override getSeqIndices(loci: Loci): OrderedSet {
        const { structure, units } = this.data;
        if (StructureElement.Loci.is(loci)) {
            if (!Structure.areRootsEquivalent(loci.structure, structure)) return Interval.Empty;
            loci = StructureElement.Loci.remap(loci, structure);

            for (const e of loci.elements) {
                const indices = this.unitIndices.get(e.unit.id);
                if (indices) {
                    if (OrderedSet.isSubset(indices, e.indices)) {
                        // this assumes this SequnceWrapper has only a single unit, otherwise we'd need to collect indices from all units and apply offsets
                        return e.indices;
                    }
                }
            }
        } else if (Structure.isLoci(loci)) {
            if (!Structure.areRootsEquivalent(loci.structure, structure)) return Interval.Empty;

            for (let i = 0, il = units.length; i < il; ++i) {
                const indices = this.unitIndices.get(units[i].id)!;
                return indices;
            }
        }
        return Interval.Empty;
    }

    getLoci(seqIdx: number) {
        const { units } = this.data;
        const lociElements: StructureElement.Loci['elements'][0][] = [];
        let offset = 0;
        for (let i = 0, il = units.length; i < il; ++i) {
            const unit = units[i];
            if (seqIdx < offset + unit.elements.length) {
                lociElements.push({ unit, indices: Interval.ofSingleton(seqIdx - offset) });
                break;
            }
            offset += unit.elements.length;
        }
        return StructureElement.Loci(this.data.structure, lociElements);
    }

    constructor(data: StructureUnit) {
        let length = 0;

        const unitIndices = new Map<number, Interval<StructureElement.UnitIndex>>();
        const lociElements: StructureElement.Loci['elements'][0][] = [];

        for (let i = 0, il = data.units.length; i < il; ++i) {
            const unit = data.units[i];
            length += unit.elements.length;

            const indices = Interval.ofBounds(0, unit.elements.length);
            unitIndices.set(unit.id, indices);
            lociElements.push({ unit, indices });
        }

        super(data, length);

        this.unitIndices = unitIndices;
    }
}
