/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Structure, Unit } from '../../structure';
import { Vec3 } from '../../../../mol-math/linear-algebra';
import { QueryFn, QueryContext } from '../context';

export function checkStructureMinMaxDistance(ctx: QueryContext, a: Structure, b: Structure, minDist: number, maxDist: number, elementRadius: QueryFn<number>) {
    if (a.elementCount === 0 || b.elementCount === 0) return true;

    if (a.elementCount <= b.elementCount) return MinMaxDist.check(ctx, a, b, minDist, maxDist, elementRadius);
    return MinMaxDist.check(ctx, b, a, minDist, maxDist, elementRadius);
}

export function checkStructureMaxRadiusDistance(ctx: QueryContext, a: Structure, b: Structure, maxDist: number, elementRadius: QueryFn<number>) {
    if (a.elementCount === 0 || b.elementCount === 0) return true;

    if (a.elementCount <= b.elementCount) return MaxRadiusDist.check(ctx, a, b, maxDist, elementRadius);
    return MaxRadiusDist.check(ctx, b, a, maxDist, elementRadius);
}

namespace MinMaxDist {
    const enum Result {
        BelowMin,
        WithinMax,
        Miss
    }

    const distVec = Vec3.zero();
    function inUnit(ctx: QueryContext, unit: Unit, p: Vec3, eRadius: number, minDist: number, maxDist: number, elementRadius: QueryFn<number>) {
        const { elements, conformation: { position } } = unit, dV = distVec;
        ctx.element.unit = unit;
        let withinRange = false;
        for (let i = 0, _i = elements.length; i < _i; i++) {
            const e = elements[i];
            ctx.element.element = e;
            const d = Math.max(0, Vec3.distance(p, position(e, dV)) - eRadius - elementRadius(ctx));
            if (d < minDist) return Result.BelowMin;
            if (d < maxDist) withinRange = true;
        }
        return withinRange ? Result.WithinMax : Result.Miss;
    }

    export function toPoint(ctx: QueryContext, s: Structure, point: Vec3, radius: number, minDist: number, maxDist: number, elementRadius: QueryFn<number>) {
        const { units } = s;
        let withinRange = false;
        for (let i = 0, _i = units.length; i < _i; i++) {
            const iu = inUnit(ctx, units[i], point, radius, minDist, maxDist, elementRadius);
            if (iu === Result.BelowMin) return Result.BelowMin;
            if (iu === Result.WithinMax) withinRange = true;
        }
        return withinRange ? Result.WithinMax : Result.Miss;
    }

    const distPivot = Vec3.zero();
    export function check(ctx: QueryContext, a: Structure, b: Structure, minDist: number, maxDist: number, elementRadius: QueryFn<number>) {
        if (a.elementCount === 0 || b.elementCount === 0) return 0;

        const { units } = a;
        let withinRange = false;
        ctx.element.structure = a;
        for (let i = 0, _i = units.length; i < _i; i++) {
            const unit = units[i];
            const { elements, conformation: { position } } = unit;
            ctx.element.unit = unit;
            for (let i = 0, _i = elements.length; i < _i; i++) {
                const e = elements[i];
                ctx.element.element = e;
                const tp = toPoint(ctx, b, position(e, distPivot), elementRadius(ctx), minDist, maxDist, elementRadius);
                if (tp === Result.BelowMin) return Result.BelowMin;
                if (tp === Result.WithinMax) withinRange = true;
            }
        }
        return withinRange;
    }
}

namespace MaxRadiusDist {
    const distVec = Vec3.zero();
    function inUnit(ctx: QueryContext, unit: Unit, p: Vec3, eRadius: number, maxDist: number, elementRadius: QueryFn<number>) {
        const { elements, conformation: { position } } = unit, dV = distVec;
        ctx.element.unit = unit;
        for (let i = 0, _i = elements.length; i < _i; i++) {
            const e = elements[i];
            ctx.element.element = e;
            if (Math.max(0, Vec3.distance(p, position(e, dV)) - eRadius - elementRadius(ctx)) <= maxDist) return true;
        }
        return false;
    }

    export function toPoint(ctx: QueryContext, s: Structure, point: Vec3, radius: number, maxDist: number, elementRadius: QueryFn<number>) {
        const { units } = s;
        for (let i = 0, _i = units.length; i < _i; i++) {
            if (inUnit(ctx, units[i], point, radius, maxDist, elementRadius)) return true;
        }
        return false;
    }

    const distPivot = Vec3.zero();
    export function check(ctx: QueryContext, a: Structure, b: Structure, maxDist: number, elementRadius: QueryFn<number>) {
        if (a.elementCount === 0 || b.elementCount === 0) return 0;

        const { units } = a;
        ctx.element.structure = a;
        for (let i = 0, _i = units.length; i < _i; i++) {
            const unit = units[i];
            ctx.element.unit = unit;
            const { elements, conformation: { position } } = unit;
            for (let i = 0, _i = elements.length; i < _i; i++) {
                const e = elements[i];
                ctx.element.element = e;
                if (toPoint(ctx, b, position(e, distPivot), elementRadius(ctx), maxDist, elementRadius)) return true;
            }
        }
        return false;
    }
}