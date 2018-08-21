/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { isSuperset } from 'mol-util/set';
import { Unit } from '../../structure';
import { QueryContext, QueryFn, QueryPredicate } from '../context';
import { StructureQuery } from '../query';
import { StructureSelection } from '../selection';
import { structureAreIntersecting } from '../utils/structure';
import { Vec3 } from 'mol-math/linear-algebra';

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

export interface UnitTypeProperties { atomic?: QueryFn, coarse?: QueryFn }

export function getCurrentStructureProperties(ctx: QueryContext, props: UnitTypeProperties, set: Set<any>) {
    const { units } = ctx.currentStructure;
    const l = ctx.pushCurrentElement();

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
            if (isSuperset(currentProps, propSet)) {
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
    atomRadius?: QueryFn<number>,
    invert?: boolean
}

function _zeroRadius(ctx: QueryContext) { return 0; }

export function within(params: WithinParams): StructureQuery {
    return queryCtx => {
        const ctx: WithinContext = {
            queryCtx,
            selection: params.query(queryCtx),
            target: params.target(queryCtx),
            maxRadius: params.maxRadius,
            minRadius: params.minRadius ? params.minRadius : 0,
            atomRadius: params.atomRadius || _zeroRadius,
            invert: !!params.invert,
        }

        if (ctx.minRadius === 0 && ctx.atomRadius === _zeroRadius) {
            return withinMaxRadius(ctx);
        } else {
            // TODO
            throw 'not implemented';
            // return withinMinMaxRadius(ctx);
        }
    }
}

interface WithinContext {
    queryCtx: QueryContext,
    selection: StructureSelection,
    target: StructureSelection,
    minRadius: number,
    maxRadius: number,
    invert: boolean,
    atomRadius: QueryFn<number>
}
function withinMaxRadius({ queryCtx, selection, target, maxRadius, invert }: WithinContext) {
    const targetLookup = StructureSelection.unionStructure(target).lookup3d;
    const ret = StructureSelection.LinearBuilder(queryCtx.inputStructure);

    const pos = Vec3.zero();
    StructureSelection.forEach(selection, (s, sI) => {
        const { units } = s;

        let withinRadius = false;
        for (let i = 0, _i = units.length; i < _i; i++) {
            const unit = units[i];
            const { elements, conformation } = unit;

            switch (unit.kind) {
                case Unit.Kind.Atomic:
                // TODO: assign radius to gaussian elements?
                case Unit.Kind.Gaussians:
                    for (let i = 0, _i = elements.length; i < _i; i++) {
                        conformation.position(elements[i], pos);
                        if (targetLookup.check(pos[0], pos[1], pos[2], maxRadius)) {
                            withinRadius = true;
                            break;
                        }
                    }
                    break;
                case Unit.Kind.Spheres:
                    const radius = unit.coarseConformation.radius;
                    for (let i = 0, _i = elements.length; i < _i; i++) {
                        conformation.position(elements[i], pos);
                        if (targetLookup.check(pos[0], pos[1], pos[2], maxRadius + radius[elements[i]])) {
                            withinRadius = true;
                            break;
                        }
                    }
                    break;
            }
            if (withinRadius) break;
        }
        if (invert) withinRadius = !withinRadius;
        if (withinRadius) ret.add(s);
        if (sI % 10 === 0) queryCtx.throwIfTimedOut();
    });

    return ret.getSelection();
}


// TODO: isConnectedTo