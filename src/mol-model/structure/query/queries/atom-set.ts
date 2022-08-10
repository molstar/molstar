/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Koya Sakuma
 * Adapted from MolQL implemtation of atom-set.ts
 *
 * Copyright (c) 2017 MolQL contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */
/*
import Environment from '../environment'
import Expression from '../expression'
import Context from '../context'
import AtomSet from '../../data/atom-set'
import AtomSelection from '../../data/atom-selection'
import ElementAddress from '../../data/element-address'
import { FastSet } from '../../../utils/collections'
import { getAtomSetProperties } from './filters'
*/

import { StructureQuery } from '../query';
import { StructureSelection } from '../selection';
import { getCurrentStructureProperties, UnitTypeProperties } from './filters';
import { QueryContext } from '../context';
// import { none } from './generators';
// import { HashSet } from '../../../../mol-data/generic';
// import { Structure } from '../../structure';


/*
import { SetUtils } from '../../../../mol-util/set';
import { Unit } from '../../structure';
import { QueryContext, QueryFn, QueryPredicate } from '../context';
import { StructureQuery } from '../query';
import { StructureSelection } from '../selection';
import { structureAreIntersecting } from '../utils/structure-set';
import { Vec3 } from '../../../../mol-math/linear-algebra';
import { checkStructureMaxRadiusDistance, checkStructureMinMaxDistance } from '../utils/structure-distance';
import { Structure } from '../../structure/structure';
import { StructureElement } from '../../structure/element';
import { SortedArray } from '../../../../mol-data/int';
*/
/*
export function pick(query: StructureQuery, pred: QueryPredicate): StructureQuery {
    return ctx => {
        const sel = query(ctx);
        const ret = StructureSelection.LinearBuilder(ctx.inputStructure);
        ctx.pushCurrentElement();
        StructureSelection.forEach(sel, (s, i) => {
            ctx.currentStructure = s;
            if (pred(ctx)) ret.add(s);
            if (i % 100) ctx.throwIfTimedOut();
        });
        ctx.popCurrentStructure();
        return ret.getSelection();
    };
}

export function pick(env: Environment, selection: Selection, pred: Expression<boolean>) {
    const sel = selection(env);
    const ret = AtomSelection.linearBuilder();

    Environment.lockSlot(env, 'atomSet');
    const { slots } = env;
    for (const atomSet of AtomSelection.atomSets(sel)) {
        slots.atomSet = atomSet;
        if (pred(env)) ret.add(atomSet);
    }
    Environment.unlockSlot(env, 'atomSet');
    return ret.getSelection();
}

*/




// xport function merge(queries: ArrayLike<StructureQuery>): S

// export function atomCount(env: Environment) {
//    return AtomSet.count(env.slots.atomSet);
// }

// export function atomCount(query : StructureSelection) : StructureQuery {
export function atomCount(ctx: QueryContext) {
    return (ctx: any) => {
        const x: number = StructureSelection.structureCount(ctx);
        return x;
    };
}

// export function countQuery(env: Environment, query: Expression<AtomSelection>) {
//    const sel = query(Environment(Context.ofAtomSet(env.context, env.slots.atomSet)))
//    return AtomSelection.atomSets(sel).length;
// }

export function countQuery(ctx: QueryContext, query: StructureQuery) {
    return (ctx: any) => {
        const sel = query(ctx);
        const x: number = StructureSelection.structureCount(sel);
        return x;
    };
}


export function propertySet(ctx: QueryContext, prop: UnitTypeProperties) {
    const set = new Set();
    const x = getCurrentStructureProperties(ctx, prop, set);
    return x;
}

/*
unction getSelectionProperties(ctx: QueryContext, query: StructureQuery, props: UnitTypeProperties) {



    const sel = query(ctx);
    ctx.pushCurrentElement();
    StructureSelection.forEach(sel, (s, i) => {
	ctx.currentStructure = s;
        getCurrentStructureProperties(ctx, props, set);

        if (i % 10) ctx.throwIfTimedOut();
    });
    ctx.popCurrentElement();
    return set;
*/
