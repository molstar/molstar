/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { SetUtils } from '../../../../mol-util/set';
import { Unit } from '../../structure';
import { QueryContext, QueryFn, QueryPredicate } from '../context';
import { StructureQuery } from '../query';
import { StructureSelection } from '../selection';
import { structureAreIntersecting } from '../utils/structure-set';
import { Vec3 } from '../../../../mol-math/linear-algebra';
import { checkStructureMaxRadiusDistance, checkStructureMinMaxDistance } from '../utils/structure-distance';
import Structure from '../../structure/structure';
import StructureElement from '../../structure/element';
import { SortedArray } from '../../../../mol-data/int';

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

export function first(query: StructureQuery): StructureQuery {
    return ctx => {
        const sel = query(ctx);
        const ret = StructureSelection.LinearBuilder(ctx.inputStructure);
        if (sel.kind === 'singletons') {
            if (sel.structure.elementCount > 0) {
                const u = sel.structure.units[0];
                const s = Structure.create([u.getChild(SortedArray.ofSingleton(u.elements[0]))], { parent: ctx.inputStructure });
                ret.add(s);
            }
        } else {
            if (sel.structures.length > 0) {
                ret.add(sel.structures[0]);
            }
        }
        return ret.getSelection();
    };
}

export interface UnitTypeProperties { atomic?: QueryFn, coarse?: QueryFn }

export function getCurrentStructureProperties(ctx: QueryContext, props: UnitTypeProperties, set: Set<any>) {
    const { units } = ctx.currentStructure;
    const l = ctx.pushCurrentElement();

    l.structure = ctx.currentStructure;
    for (const unit of units) {
        l.unit = unit;
        const elements = unit.elements;

        let fn;
        if (Unit.isAtomic(unit)) fn = props.atomic;
        else fn = props.coarse;
        if (!fn) continue;

        for (let j = 0, _j = elements.length; j < _j; j++) {
            l.element = elements[j];
            set.add(fn(ctx));
        }

        ctx.throwIfTimedOut();
    }
    ctx.popCurrentElement();
    return set;
}

function getSelectionProperties(ctx: QueryContext, query: StructureQuery, props: UnitTypeProperties) {
    const set = new Set();

    const sel = query(ctx);
    ctx.pushCurrentElement();
    StructureSelection.forEach(sel, (s, i) => {
        ctx.currentStructure = s;
        getCurrentStructureProperties(ctx, props, set);

        if (i % 10) ctx.throwIfTimedOut();
    });
    ctx.popCurrentElement();
    return set;
}

export function withSameAtomProperties(query: StructureQuery, propertySource: StructureQuery, props: UnitTypeProperties): StructureQuery {
    return ctx => {
        const sel = query(ctx);
        const propSet = getSelectionProperties(ctx, propertySource, props);

        const ret = StructureSelection.LinearBuilder(ctx.inputStructure);
        ctx.pushCurrentStructure();
        StructureSelection.forEach(sel, (s, i) => {
            ctx.currentStructure = s;
            const currentProps = getCurrentStructureProperties(ctx, props, new Set());
            if (SetUtils.isSuperset(currentProps, propSet)) {
                ret.add(s);
            }

            if (i % 10) ctx.throwIfTimedOut();
        });
        ctx.popCurrentStructure();
        return ret.getSelection();
    };
}

export function areIntersectedBy(query: StructureQuery, by: StructureQuery): StructureQuery {
    return ctx => {
        const mask = StructureSelection.unionStructure(by(ctx));
        const ret = StructureSelection.LinearBuilder(ctx.inputStructure);

        StructureSelection.forEach(query(ctx), (s, i) => {
            if (structureAreIntersecting(mask, s)) ret.add(s);
            if (i % 10) ctx.throwIfTimedOut();
        });
        return ret.getSelection();
    };
}

export interface WithinParams {
    query: StructureQuery,
    target: StructureQuery,
    minRadius?: number,
    maxRadius: number,
    elementRadius?: QueryFn<number>,
    invert?: boolean
}

export function within(params: WithinParams): StructureQuery {
    return queryCtx => {
        const ctx: WithinContext = {
            queryCtx,
            selection: params.query(queryCtx),
            target: params.target(queryCtx),
            maxRadius: params.maxRadius,
            minRadius: params.minRadius ? Math.max(0, params.minRadius) : 0,
            elementRadius: params.elementRadius!,
            invert: !!params.invert,
        };

        if (ctx.minRadius === 0 && typeof params.minRadius === 'undefined') {
            return withinMaxRadiusLookup(ctx);
        } else if (ctx.minRadius === 0) {
            return withinMaxRadius(ctx);
        } else {
            return withinMinMaxRadius(ctx);
        }
    };
}

interface WithinContext {
    queryCtx: QueryContext,
    selection: StructureSelection,
    target: StructureSelection,
    minRadius: number,
    maxRadius: number,
    invert: boolean,
    elementRadius: QueryFn<number>
}
function withinMaxRadiusLookup({ queryCtx, selection, target, maxRadius, invert }: WithinContext) {
    const targetLookup = StructureSelection.unionStructure(target).lookup3d;
    const ret = StructureSelection.LinearBuilder(queryCtx.inputStructure);

    const pos = Vec3.zero();
    StructureSelection.forEach(selection, (s, sI) => {
        const { units } = s;
        let withinRadius = false;
        for (let i = 0, _i = units.length; i < _i; i++) {
            const unit = units[i];
            const { elements, conformation: { position, r } } = unit;

            for (let i = 0, _i = elements.length; i < _i; i++) {
                const e = elements[i];
                position(e, pos);
                if (targetLookup.check(pos[0], pos[1], pos[2], maxRadius + r(e))) {
                    withinRadius = true;
                    break;
                }
            }
            if (withinRadius) break;
        }
        if (invert) withinRadius = !withinRadius;
        if (withinRadius) ret.add(s);
        if (sI % 10 === 0) queryCtx.throwIfTimedOut();
    });

    return ret.getSelection();
}

function withinMaxRadius({ queryCtx, selection, target, maxRadius, invert, elementRadius }: WithinContext) {
    const targetStructure = StructureSelection.unionStructure(target);
    const ret = StructureSelection.LinearBuilder(queryCtx.inputStructure);

    queryCtx.pushCurrentElement();
    StructureSelection.forEach(selection, (s, sI) => {
        let withinRadius = checkStructureMaxRadiusDistance(queryCtx, targetStructure, s, maxRadius, elementRadius);
        if (invert) withinRadius = !withinRadius;
        if (withinRadius) ret.add(s);
        if (sI % 10 === 0) queryCtx.throwIfTimedOut();
    });
    queryCtx.popCurrentElement();

    return ret.getSelection();
}

function withinMinMaxRadius({ queryCtx, selection, target, minRadius, maxRadius, invert, elementRadius }: WithinContext) {
    const targetStructure = StructureSelection.unionStructure(target);
    const ret = StructureSelection.LinearBuilder(queryCtx.inputStructure);

    queryCtx.pushCurrentElement();
    StructureSelection.forEach(selection, (s, sI) => {
        let withinRadius = checkStructureMinMaxDistance(queryCtx, targetStructure, s, minRadius, maxRadius, elementRadius);
        if (invert) withinRadius = !withinRadius;
        if (withinRadius) ret.add(s);
        if (sI % 10 === 0) queryCtx.throwIfTimedOut();
    });
    queryCtx.popCurrentElement();

    return ret.getSelection();
}

interface IsConnectedToCtx {
    queryCtx: QueryContext,
    disjunct: boolean,
    input: Structure,
    target: Structure
}

function checkConnected(ctx: IsConnectedToCtx, structure: Structure) {
    const { queryCtx, input, target, disjunct } = ctx;
    const atomicBond = queryCtx.atomicBond;

    const interBonds = input.interUnitBonds;

    atomicBond.setStructure(input);

    for (const unit of structure.units) {
        if (!Unit.isAtomic(unit)) continue;

        const inputUnit = input.unitMap.get(unit.id) as Unit.Atomic;

        const { offset, b, edgeProps: { flags, order } } = inputUnit.bonds;
        const bondedUnits = interBonds.getConnectedUnits(unit);
        const buCount = bondedUnits.length;

        const srcElements = unit.elements;
        const inputElements = inputUnit.elements;

        for (let i = 0 as StructureElement.UnitIndex, _i = srcElements.length; i < _i; i++) {
            const inputIndex = SortedArray.indexOf(inputElements, srcElements[i]) as StructureElement.UnitIndex;

            atomicBond.a.unit = inputUnit;
            atomicBond.b.unit = inputUnit;

            // tElement.unit = unit;
            for (let l = offset[inputIndex], _l = offset[inputIndex + 1]; l < _l; l++) {
                // tElement.element = inputElements[b[l]];
                atomicBond.b.element = inputUnit.elements[b[l]];
                if (disjunct && SortedArray.has(unit.elements, atomicBond.b.element)) continue;
                if (!target.hasElement(atomicBond.b)) continue;

                atomicBond.aIndex = inputIndex;
                atomicBond.a.element = srcElements[i];
                atomicBond.bIndex = b[l] as StructureElement.UnitIndex;
                atomicBond.type = flags[l];
                atomicBond.order = order[l];
                if (atomicBond.test(queryCtx, true)) return true;
            }

            for (let li = 0; li < buCount; li++) {
                const lu = bondedUnits[li];
                const bElements = lu.unitB.elements;
                const bonds = lu.getEdges(inputIndex);
                for (let bi = 0, _bi = bonds.length; bi < _bi; bi++) {
                    const bond = bonds[bi];
                    atomicBond.b.unit = lu.unitB;
                    atomicBond.b.element = bElements[bond.indexB];
                    if (!target.hasElement(atomicBond.b)) continue;
                    if (disjunct && structure.hasElement(atomicBond.b)) continue;

                    atomicBond.a.unit = inputUnit;
                    atomicBond.aIndex = inputIndex;
                    atomicBond.a.element = srcElements[i];

                    atomicBond.bIndex = bond.indexB;
                    atomicBond.type = bond.props.flag;
                    atomicBond.order = bond.props.order;
                    if (atomicBond.test(queryCtx, true)) return true;
                }
            }
        }
    }

    return false;
}

export interface IsConnectedToParams {
    query: StructureQuery,
    target: StructureQuery,
    bondTest?: QueryFn<boolean>,
    disjunct: boolean,
    invert: boolean
}

export function isConnectedTo({ query, target, disjunct, invert, bondTest }: IsConnectedToParams): StructureQuery {
    return ctx => {
        const targetSel = target(ctx);
        if (StructureSelection.isEmpty(targetSel)) return targetSel;
        const selection = query(ctx);
        if (StructureSelection.isEmpty(selection)) return selection;

        const connCtx: IsConnectedToCtx = {
            queryCtx: ctx,
            input: ctx.inputStructure,
            disjunct,
            target: StructureSelection.unionStructure(targetSel)
        };

        const ret = StructureSelection.LinearBuilder(ctx.inputStructure);
        ctx.pushCurrentBond();
        ctx.atomicBond.setTestFn(bondTest);

        StructureSelection.forEach(selection, (s, sI) => {
            if (checkConnected(connCtx, s)) {
                ret.add(s);
            } else if (invert) {
                ret.add(s);
            }
            if (sI % 5 === 0) ctx.throwIfTimedOut();
        });
        ctx.popCurrentBond();

        return ret.getSelection();
    };
}