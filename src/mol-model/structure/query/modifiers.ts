/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Segmentation } from 'mol-data/int';
import { RuntimeContext } from 'mol-task';
import { Structure, Unit, Element } from '../structure';
import Query from './query';
import Selection from './selection';
import { UniqueStructuresBuilder } from './utils/builders';
import { StructureUniqueSubsetBuilder } from '../structure/util/unique-subset-builder';

function getWholeResidues(ctx: RuntimeContext, source: Structure, structure: Structure) {
    const builder = source.subsetBuilder(true);
    for (const unit of structure.units) {
        if (unit.kind !== Unit.Kind.Atomic) {
            // just copy non-atomic units.
            builder.setUnit(unit.id, unit.elements);
            continue;
        }

        const { residueSegments } = unit.model.atomicHierarchy;

        const elements = unit.elements;
        builder.beginUnit(unit.id);
        const residuesIt = Segmentation.transientSegments(residueSegments, elements);
        while (residuesIt.hasNext) {
            const rI = residuesIt.move().index;
            for (let j = residueSegments.segments[rI], _j = residueSegments.segments[rI + 1]; j < _j; j++) {
                builder.addElement(j as Element);
            }
        }
        builder.commitUnit();
    }
    return builder.getStructure();
}

export function wholeResidues(query: Query.Provider, isFlat: boolean): Query.Provider {
    return async (structure, ctx) => {
        const inner = await query(structure, ctx);
        if (Selection.isSingleton(inner)) {
            return Selection.Singletons(structure, getWholeResidues(ctx, structure, inner.structure));
        } else {
            const builder = new UniqueStructuresBuilder(structure);
            let progress = 0;
            for (const s of inner.structures) {
                builder.add(getWholeResidues(ctx, structure, s));
                progress++;
                if (ctx.shouldUpdate) await ctx.update({ message: 'Whole Residues', current: progress, max: inner.structures.length });
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

export function includeSurroundings(query: Query.Provider, params: IncludeSurroundingsParams): Query.Provider {
    return async (structure, ctx) => {
        const inner = await query(structure, ctx);
        if (Selection.isSingleton(inner)) {
            const surr = await getIncludeSurroundings(ctx, structure, inner.structure, params);
            const ret = Selection.Singletons(structure, surr);
            return ret;
        } else {
            const builder = new UniqueStructuresBuilder(structure);
            for (const s of inner.structures) {
                builder.add(await getIncludeSurroundings(ctx, structure, s, params));
            }
            return builder.getSelection();
        }
    };
}