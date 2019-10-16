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
import { LinkType } from '../../model/types';
import { BundleElement, Bundle } from '../../structure/element/bundle';
import { UnitIndex } from '../../structure/element/element';

export function defaultLinkTest(ctx: QueryContext) {
    return LinkType.isCovalent(ctx.atomicLink.type);
}

export function atomicSequence(): StructureQuery {
    return ctx => {
        const { inputStructure } = ctx;
        const l = StructureElement.Location.create();

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
    return ctx => {
        const { inputStructure } = ctx;
        const l = StructureElement.Location.create();

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
    return ctx => {
        const { inputStructure } = ctx;
        const l = StructureElement.Location.create();

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
    return ctx => {
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
    return ctx => {
        const bundle: Bundle = {
            hash: ctx.inputStructure.hashCode,
            elements
        };

        return StructureSelection.Sequence(ctx.inputStructure, [Bundle.toStructure(bundle, ctx.inputStructure)]);
    };
}