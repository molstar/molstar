/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Segmentation } from 'mol-data/int';
import StructureElement from 'mol-model/structure/structure/element';
import { StructureProperties as P, Unit } from '../../structure';
import Structure from '../../structure/structure';
import { StructureQuery } from '../query';
import { StructureSelection } from '../selection';

export function atomicSequence(): StructureQuery {
    return ctx => {
        const { inputStructure } = ctx;
        const l = StructureElement.create();

        const units: Unit[] = [];
        for (const unit of inputStructure.units) {
            if (unit.kind !== Unit.Kind.Atomic) continue;
            l.unit = unit;
            const elements = unit.elements;
            l.element = elements[0];
            if (P.entity.type(l) !== 'polymer') continue;

            const residuesIt = Segmentation.transientSegments(unit.model.atomicHierarchy.residueAtomSegments, elements);
            let residueCount = 0;
            while (residuesIt.hasNext) {
                residueCount++;
                residuesIt.move();
            }

            if (residueCount < 8) continue;

            units.push(unit);
        }
        return StructureSelection.Singletons(inputStructure, new Structure(units));
    };
}

export function water(): StructureQuery {
    return ctx => {
        const { inputStructure } = ctx;
        const l = StructureElement.create();

        const units: Unit[] = [];
        for (const unit of inputStructure.units) {
            if (unit.kind !== Unit.Kind.Atomic) continue;

            l.unit = unit;
            const elements = unit.elements;
            l.element = elements[0];
            if (P.entity.type(l) !== 'water') continue;
            units.push(unit);
        }
        return StructureSelection.Singletons(inputStructure, new Structure(units));
    };
}

export function atomicHet(): StructureQuery {
    return ctx => {
        const { inputStructure } = ctx;
        const l = StructureElement.create();

        const units: Unit[] = [];
        for (const unit of inputStructure.units) {
            if (unit.kind !== Unit.Kind.Atomic) continue;

            l.unit = unit;
            const elements = unit.elements;
            l.element = elements[0];
            if (P.entity.type(l) === 'water') continue;
            if (P.entity.type(l) === 'polymer') {
                const residuesIt = Segmentation.transientSegments(unit.model.atomicHierarchy.residueAtomSegments, elements);
                let residueCount = 0;
                while (residuesIt.hasNext) {
                    residueCount++;
                    residuesIt.move();
                }

                if (residueCount >= 8) continue;
            }

            units.push(unit);
        }
        return StructureSelection.Singletons(inputStructure, new Structure(units));
    };
}

export function spheres(): StructureQuery {
    return ctx => {
        const { inputStructure } = ctx;

        const units: Unit[] = [];
        for (const unit of inputStructure.units) {
            if (unit.kind !== Unit.Kind.Spheres) continue;
            units.push(unit);
        }
        return StructureSelection.Singletons(inputStructure, new Structure(units));
    };
}
