/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Segmentation } from 'mol-data/int';
import { Structure, Unit } from '../../structure';
import { StructureQuery } from '../query';
import { StructureSelection } from '../selection';
import { UniqueStructuresBuilder } from '../utils/builders';
import { StructureUniqueSubsetBuilder } from '../../structure/util/unique-subset-builder';
import { QueryContext } from '../context';
import { structureIntersect, structureSubtract } from '../utils/structure-set';

function getWholeResidues(ctx: QueryContext, source: Structure, structure: Structure) {
    const builder = source.subsetBuilder(true);
    for (const unit of structure.units) {
        if (unit.kind !== Unit.Kind.Atomic) {
            // just copy non-atomic units.
            builder.setUnit(unit.id, unit.elements);
            continue;
        }

        const { residueAtomSegments } = unit.model.atomicHierarchy;

        const elements = unit.elements;
        builder.beginUnit(unit.id);
        const residuesIt = Segmentation.transientSegments(residueAtomSegments, elements);
        while (residuesIt.hasNext) {
            const rI = residuesIt.move().index;
            for (let j = residueAtomSegments.offsets[rI], _j = residueAtomSegments.offsets[rI + 1]; j < _j; j++) {
                builder.addElement(j);
            }
        }
        builder.commitUnit();

        ctx.throwIfTimedOut();
    }
    return builder.getStructure();
}

export function wholeResidues(query: StructureQuery): StructureQuery {
    return ctx => {
        const inner = query(ctx);
        if (StructureSelection.isSingleton(inner)) {
            return StructureSelection.Singletons(ctx.inputStructure, getWholeResidues(ctx, ctx.inputStructure, inner.structure));
        } else {
            const builder = new UniqueStructuresBuilder(ctx.inputStructure);
            for (const s of inner.structures) {
                builder.add(getWholeResidues(ctx, ctx.inputStructure, s));
            }
            return builder.getSelection();
        }
    };
}

export interface IncludeSurroundingsParams {
    radius: number,
    // TODO
    // atomRadius?: Element.Property<number>,
    wholeResidues?: boolean
}

function getIncludeSurroundings(ctx: QueryContext, source: Structure, structure: Structure, params: IncludeSurroundingsParams) {
    const builder = new StructureUniqueSubsetBuilder(source);
    const lookup = source.lookup3d;
    const r = params.radius;

    for (const unit of structure.units) {
        const { x, y, z } = unit.conformation;
        const elements = unit.elements;
        for (let i = 0, _i = elements.length; i < _i; i++) {
            const e = elements[i];
            lookup.findIntoBuilder(x(e), y(e), z(e), r, builder);
        }

        ctx.throwIfTimedOut();
    }
    return !!params.wholeResidues ? getWholeResidues(ctx, source, builder.getStructure()) : builder.getStructure();
}

export function includeSurroundings(query: StructureQuery, params: IncludeSurroundingsParams): StructureQuery {
    return ctx => {
        const inner = query(ctx);
        if (StructureSelection.isSingleton(inner)) {
            const surr = getIncludeSurroundings(ctx, ctx.inputStructure, inner.structure, params);
            const ret = StructureSelection.Singletons(ctx.inputStructure, surr);
            return ret;
        } else {
            const builder = new UniqueStructuresBuilder(ctx.inputStructure);
            for (const s of inner.structures) {
                builder.add(getIncludeSurroundings(ctx, ctx.inputStructure, s, params));
            }
            return builder.getSelection();
        }
    };
}

export function querySelection(selection: StructureQuery, query: StructureQuery): StructureQuery {
    return ctx => {
        const targetSel = selection(ctx);
        if (StructureSelection.structureCount(targetSel) === 0) return targetSel;

        const ret = StructureSelection.UniqueBuilder(ctx.inputStructure);
        const add = (s: Structure) => ret.add(s);

        StructureSelection.forEach(targetSel, (s, sI) => {
            ctx.pushInputStructure(s);
            StructureSelection.forEach(query(ctx), add);
            ctx.popInputStructure();
            if (sI % 10 === 0) ctx.throwIfTimedOut();
        });
        return ret.getSelection();
    }
}

export function intersectBy(query: StructureQuery, by: StructureQuery): StructureQuery {
    return ctx => {
        const selection = query(ctx);
        if (StructureSelection.structureCount(selection) === 0) return selection;

        const bySel = by(ctx);
        if (StructureSelection.structureCount(bySel) === 0) return StructureSelection.Empty(ctx.inputStructure);
        const unionBy = StructureSelection.unionStructure(bySel);

        const ret = StructureSelection.UniqueBuilder(ctx.inputStructure);
        StructureSelection.forEach(selection, (s, sI) => {
            const ii = structureIntersect(unionBy, s);
            if (ii.elementCount !== 0) ret.add(ii);
        });

        return ret.getSelection();
    };
}

export function exceptBy(query: StructureQuery, by: StructureQuery): StructureQuery {
    return ctx => {
        const selection = query(ctx);
        if (StructureSelection.structureCount(selection) === 0) return selection;

        const bySel = by(ctx);
        if (StructureSelection.structureCount(bySel) === 0) return StructureSelection.Empty(ctx.inputStructure);
        const subtractBy = StructureSelection.unionStructure(bySel);

        const ret = StructureSelection.UniqueBuilder(ctx.inputStructure);
        StructureSelection.forEach(selection, (s, sI) => {
            const diff = structureSubtract(s, subtractBy);
            if (diff.elementCount !== 0) ret.add(diff);
        });

        return ret.getSelection();
    };
}

export function union(query: StructureQuery): StructureQuery {
    return ctx => {
        const ret = StructureSelection.LinearBuilder(ctx.inputStructure);
        ret.add(StructureSelection.unionStructure(query(ctx)));
        return ret.getSelection();
    };
}

// TODO: unionBy (skip this one?), cluster, includeConnected, includeSurroundings with "radii", expandProperty