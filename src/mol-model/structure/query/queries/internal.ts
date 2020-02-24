/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Segmentation, SortedArray } from '../../../../mol-data/int';
import StructureElement from '../../../../mol-model/structure/structure/element';
import { StructureProperties as P, Unit } from '../../structure';
import Structure from '../../structure/structure';
import { StructureQuery } from '../query';
import { StructureSelection } from '../selection';
import { QueryContext } from '../context';
import { BondType } from '../../model/types';
import { BundleElement, Bundle } from '../../structure/element/bundle';
import { UnitIndex } from '../../structure/element/element';

export function defaultBondTest(ctx: QueryContext) {
    return BondType.isCovalent(ctx.atomicBond.type);
}

export function atomicSequence(): StructureQuery {
    return function query_atomicSequence(ctx) {
        const { inputStructure } = ctx;
        const l = StructureElement.Location.create(inputStructure);

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
        return StructureSelection.Singletons(inputStructure, new Structure(units, { parent: inputStructure }));
    };
}

export function water(): StructureQuery {
    return function query_water(ctx) {
        const { inputStructure } = ctx;
        const l = StructureElement.Location.create(inputStructure);

        const units: Unit[] = [];
        for (const unit of inputStructure.units) {
            if (unit.kind !== Unit.Kind.Atomic) continue;

            l.unit = unit;
            const elements = unit.elements;
            l.element = elements[0];
            if (P.entity.type(l) !== 'water') continue;
            units.push(unit);
        }
        return StructureSelection.Singletons(inputStructure, new Structure(units, { parent: inputStructure }));
    };
}

export function atomicHet(): StructureQuery {
    return function query_atomicHet(ctx) {
        const { inputStructure } = ctx;
        const l = StructureElement.Location.create(inputStructure);

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
        return StructureSelection.Singletons(inputStructure, new Structure(units, { parent: inputStructure }));
    };
}

export function spheres(): StructureQuery {
    return function query_spheres(ctx) {
        const { inputStructure } = ctx;

        const units: Unit[] = [];
        for (const unit of inputStructure.units) {
            if (unit.kind !== Unit.Kind.Spheres) continue;
            units.push(unit);
        }
        return StructureSelection.Singletons(inputStructure, new Structure(units, { parent: inputStructure }));
    };
}

export function bundleElementImpl(groupedUnits: number[][], ranges: number[], set: number[]): BundleElement {
    return {
        groupedUnits: groupedUnits as any as SortedArray<number>[],
        ranges: ranges as any as SortedArray<UnitIndex>,
        set: set as any as SortedArray<UnitIndex>
    };
}

export function bundleGenerator(elements: BundleElement[]): StructureQuery {
    return function query_bundleGenerator(ctx) {
        const bundle: Bundle = {
            hash: ctx.inputStructure.hashCode,
            elements
        };

        return StructureSelection.Sequence(ctx.inputStructure, [Bundle.toStructure(bundle, ctx.inputStructure)]);
    };
}