/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Segmentation } from 'mol-data/int';
import { RuntimeContext } from 'mol-task';
import { Structure, Unit } from '../../structure';
import { StructureQuery } from '../query';
import { StructureSelection } from '../selection';
import { UniqueStructuresBuilder } from '../utils/builders';
import { StructureUniqueSubsetBuilder } from '../../structure/util/unique-subset-builder';

function getWholeResidues(ctx: RuntimeContext, source: Structure, structure: Structure) {
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
    }
    return builder.getStructure();
}

export function wholeResidues(query: StructureQuery, isFlat: boolean): StructureQuery {
    return async (ctx) => {
        const inner = await query(ctx);
        if (StructureSelection.isSingleton(inner)) {
            return StructureSelection.Singletons(ctx.inputStructure, getWholeResidues(ctx.taskCtx, ctx.inputStructure, inner.structure));
        } else {
            const builder = new UniqueStructuresBuilder(ctx.inputStructure);
            let progress = 0;
            for (const s of inner.structures) {
                builder.add(getWholeResidues(ctx.taskCtx, ctx.inputStructure, s));
                progress++;
                if (ctx.taskCtx.shouldUpdate) await ctx.taskCtx.update({ message: 'Whole Residues', current: progress, max: inner.structures.length });
            }
            return builder.getSelection();
        }
    };
}


// export function groupBy()  ...

export interface IncludeSurroundingsParams {
    radius: number,
    // atomRadius?: Element.Property<number>,
    wholeResidues?: boolean
}

async function getIncludeSurroundings(ctx: RuntimeContext, source: Structure, structure: Structure, params: IncludeSurroundingsParams) {
    const builder = new StructureUniqueSubsetBuilder(source);
    const lookup = source.lookup3d;
    const r = params.radius;

    let progress = 0;
    for (const unit of structure.units) {
        const { x, y, z } = unit.conformation;
        const elements = unit.elements;
        for (let i = 0, _i = elements.length; i < _i; i++) {
            const e = elements[i];
            lookup.findIntoBuilder(x(e), y(e), z(e), r, builder);
        }
        progress++;
        if (progress % 2500 === 0 && ctx.shouldUpdate) await ctx.update({ message: 'Include Surroudnings', isIndeterminate: true });
    }
    return !!params.wholeResidues ? getWholeResidues(ctx, source, builder.getStructure()) : builder.getStructure();
}

export function includeSurroundings(query: StructureQuery, params: IncludeSurroundingsParams): StructureQuery {
    return async (ctx) => {
        const inner = await query(ctx);
        if (StructureSelection.isSingleton(inner)) {
            const surr = await getIncludeSurroundings(ctx.taskCtx, ctx.inputStructure, inner.structure, params);
            const ret = StructureSelection.Singletons(ctx.inputStructure, surr);
            return ret;
        } else {
            const builder = new UniqueStructuresBuilder(ctx.inputStructure);
            for (const s of inner.structures) {
                builder.add(await getIncludeSurroundings(ctx.taskCtx, ctx.inputStructure, s, params));
            }
            return builder.getSelection();
        }
    };
}