/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { StructureQuery } from '../query';
import { StructureSelection } from '../selection';
import { none } from './generators';
import { HashSet } from '../../../../mol-data/generic';
import { Structure } from '../../structure';

export function merge(queries: ArrayLike<StructureQuery>): StructureQuery {
    if (queries.length === 0) {
        return none;
    } else if (queries.length === 1) {
        return queries[0];
    }
    return ctx => {
        const ret = StructureSelection.UniqueBuilder(ctx.inputStructure);
        for (let i = 0; i < queries.length; i++) {
            StructureSelection.forEach(queries[i](ctx), (s, j) => {
                ret.add(s);
                if (i % 100) ctx.throwIfTimedOut();
            });
        }
        return ret.getSelection();
    };
}

export function intersect(queries: ArrayLike<StructureQuery>): StructureQuery {
    if (queries.length === 0) {
        return none;
    } else if (queries.length === 1) {
        return queries[0];
    }
    return ctx => {
        const selections: StructureSelection[] = [];
        for (let i = 0; i < queries.length; i++) selections.push(queries[i](ctx));
        let pivotIndex = 0, pivotLength = StructureSelection.structureCount(selections[0]);
        for (let i = 1; i < selections.length; i++) {
            const len = StructureSelection.structureCount(selections[i]);
            if (len < pivotLength) {
                pivotIndex = i;
                pivotLength = len;
            }
        }

        ctx.throwIfTimedOut();
        const pivotSet = HashSet<Structure>(s => s.hashCode, Structure.areUnitAndIndicesEqual);
        StructureSelection.forEach(selections[pivotIndex], s => pivotSet.add(s));

        const ret = StructureSelection.UniqueBuilder(ctx.inputStructure);

        for (let pI = 0; pI < selections.length; pI++) {
            if (pI === pivotIndex) continue;
            StructureSelection.forEach(selections[pI], s => {
                if (pivotSet.has(s)) ret.add(s);
            });
            ctx.throwIfTimedOut();
        }

        return ret.getSelection();
    };
}

// TODO: distanceCluster