/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Segmentation, SortedArray } from '../../../../mol-data/int';
import { Structure, Unit } from '../../structure';
import { StructureQuery } from '../query';
import { StructureSelection } from '../selection';
import { UniqueStructuresBuilder } from '../utils/builders';
import { StructureUniqueSubsetBuilder } from '../../structure/util/unique-subset-builder';
import { QueryContext, QueryFn } from '../context';
import { structureIntersect, structureSubtract } from '../utils/structure-set';
import { UniqueArray } from '../../../../mol-data/generic';
import { StructureSubsetBuilder } from '../../structure/util/subset-builder';
import StructureElement from '../../structure/element';
import { defaultLinkTest } from './internal';

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
    elementRadius?: QueryFn<number>,
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

interface IncludeSurroundingsParamsWithRadius extends IncludeSurroundingsParams {
    elementRadius: QueryFn<number>,
    elementRadiusClosure: StructureElement.Property<number>,
    sourceMaxRadius: number
}

function getIncludeSurroundingsWithRadius(ctx: QueryContext, source: Structure, structure: Structure, params: IncludeSurroundingsParamsWithRadius) {
    const builder = new StructureUniqueSubsetBuilder(source);
    const lookup = source.lookup3d;
    const { elementRadius, elementRadiusClosure, sourceMaxRadius, radius } = params;

    ctx.pushCurrentElement();
    for (const unit of structure.units) {
        ctx.element.unit = unit;
        const { x, y, z } = unit.conformation;
        const elements = unit.elements;

        for (let i = 0, _i = elements.length; i < _i; i++) {
            const e = elements[i];
            ctx.element.element = e;
            const eRadius = elementRadius(ctx);
            lookup.findIntoBuilderWithRadius(x(e), y(e), z(e), eRadius, sourceMaxRadius, radius, elementRadiusClosure, builder);
        }

        ctx.throwIfTimedOut();
    }

    ctx.popCurrentElement();
    return !!params.wholeResidues ? getWholeResidues(ctx, source, builder.getStructure()) : builder.getStructure();
}

function createElementRadiusFn(ctx: QueryContext, eRadius: QueryFn<number>): StructureElement.Property<number> {
    return e => {
        ctx.element.unit = e.unit;
        ctx.element.element = e.element;
        return eRadius(ctx);
    }
}

function findStructureRadius(ctx: QueryContext, eRadius: QueryFn<number>) {
    let r = 0;
    for (const unit of ctx.inputStructure.units) {
        ctx.element.unit = unit;
        const elements = unit.elements;

        for (let i = 0, _i = elements.length; i < _i; i++) {
            const e = elements[i];
            ctx.element.element = e;
            const eR = eRadius(ctx);
            if (eR > r) r = eR;
        }

    }
    ctx.throwIfTimedOut();
    return r;
}

export function includeSurroundings(query: StructureQuery, params: IncludeSurroundingsParams): StructureQuery {
    return ctx => {
        const inner = query(ctx);

        if (params.elementRadius) {
            const prms: IncludeSurroundingsParamsWithRadius = {
                ...params,
                elementRadius: params.elementRadius,
                elementRadiusClosure: createElementRadiusFn(ctx, params.elementRadius),
                sourceMaxRadius: findStructureRadius(ctx, params.elementRadius)
            };

            if (StructureSelection.isSingleton(inner)) {
                const surr = getIncludeSurroundingsWithRadius(ctx, ctx.inputStructure, inner.structure, prms);
                const ret = StructureSelection.Singletons(ctx.inputStructure, surr);
                return ret;
            } else {
                const builder = new UniqueStructuresBuilder(ctx.inputStructure);
                for (const s of inner.structures) {
                    builder.add(getIncludeSurroundingsWithRadius(ctx, ctx.inputStructure, s, prms));
                }
                return builder.getSelection();
            }
        }

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
            if (sI % 50 === 0) ctx.throwIfTimedOut();
        });

        return ret.getSelection();
    };
}

export function exceptBy(query: StructureQuery, by: StructureQuery): StructureQuery {
    return ctx => {
        const selection = query(ctx);
        if (StructureSelection.structureCount(selection) === 0) return selection;

        const bySel = by(ctx);
        if (StructureSelection.structureCount(bySel) === 0) return selection;
        const subtractBy = StructureSelection.unionStructure(bySel);

        const ret = StructureSelection.UniqueBuilder(ctx.inputStructure);
        StructureSelection.forEach(selection, (s, sI) => {
            const diff = structureSubtract(s, subtractBy);
            if (diff.elementCount !== 0) ret.add(diff);
            if (sI % 50 === 0) ctx.throwIfTimedOut();
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

export function expandProperty(query: StructureQuery, property: QueryFn): StructureQuery {
    return ctx => {
        const src = query(ctx);
        const propertyToStructureIndexMap = new Map<any, UniqueArray<number>>();

        const builders: StructureSubsetBuilder[] = [];
        ctx.pushCurrentElement();
        StructureSelection.forEach(src, (s, sI) => {
            for (const unit of s.units) {
                ctx.element.unit = unit;
                const elements = unit.elements;
                for (let i = 0, _i = elements.length; i < _i; i++) {
                    ctx.element.element = elements[i];
                    const p = property(ctx);
                    let arr: UniqueArray<number>;
                    if (propertyToStructureIndexMap.has(p)) arr = propertyToStructureIndexMap.get(p)!;
                    else {
                        arr = UniqueArray.create<number>();
                        propertyToStructureIndexMap.set(p, arr);
                    }
                    UniqueArray.add(arr, sI, sI);
                }
            }
            builders[sI] = ctx.inputStructure.subsetBuilder(true);

            if (sI % 10 === 0) ctx.throwIfTimedOut();
        });

        for (const unit of ctx.inputStructure.units) {
            ctx.element.unit = unit;
            const elements = unit.elements;
            for (let i = 0, _i = elements.length; i < _i; i++) {
                ctx.element.element = elements[i];
                const p = property(ctx);
                if (!propertyToStructureIndexMap.has(p)) continue;
                const indices = propertyToStructureIndexMap.get(p)!.array;

                for (let _sI = 0, __sI = indices.length; _sI < __sI; _sI++) {
                    builders[indices[_sI]].addToUnit(unit.id, elements[i]);
                }
            }
        }

        ctx.popCurrentElement();

        const ret = StructureSelection.UniqueBuilder(ctx.inputStructure);
        for (const b of builders) ret.add(b.getStructure());

        return ret.getSelection();
    };
}

export interface IncludeConnectedParams {
    query: StructureQuery,
    linkTest?: QueryFn<boolean>,
    layerCount: number,
    wholeResidues: boolean
}

export function includeConnected({ query, layerCount, wholeResidues, linkTest }: IncludeConnectedParams): StructureQuery {
    const bt = linkTest || defaultLinkTest;
    const lc = Math.max(layerCount, 0);
    return ctx => {
        const builder = StructureSelection.UniqueBuilder(ctx.inputStructure);
        const src = query(ctx);
        ctx.pushCurrentLink();
        StructureSelection.forEach(src, (s, sI) => {
            let incl = s;
            for (let i = 0; i < lc; i++) {
                incl = includeConnectedStep(ctx, bt, wholeResidues, incl);
            }
            builder.add(incl);
            if (sI % 10 === 0) ctx.throwIfTimedOut();
        });
        ctx.popCurrentLink();
        return builder.getSelection();
    }
}

function includeConnectedStep(ctx: QueryContext, linkTest: QueryFn<boolean>, wholeResidues: boolean, structure: Structure) {
    const expanded = expandConnected(ctx, structure, linkTest);
    if (wholeResidues) return getWholeResidues(ctx, ctx.inputStructure, expanded);
    return expanded;
}

function expandConnected(ctx: QueryContext, structure: Structure, linkTest: QueryFn<boolean>) {
    const inputStructure = ctx.inputStructure;
    const interLinks = inputStructure.links;
    const builder = new StructureUniqueSubsetBuilder(inputStructure);

    const processedUnits = new Set<number>();

    const atomicLink = ctx.atomicLink;

    // Process intra unit links
    for (const unit of structure.units) {
        processedUnits.add(unit.id);

        if (unit.kind !== Unit.Kind.Atomic) {
            // add the whole unit
            builder.beginUnit(unit.id);
            for (let i = 0, _i = unit.elements.length; i < _i; i++) {
                builder.addElement(unit.elements[i]);
            }
            builder.commitUnit();
            continue;
        }

        const inputUnit = inputStructure.unitMap.get(unit.id) as Unit.Atomic;
        const { offset: intraLinkOffset, b: intraLinkB, edgeProps: { flags, order } } = inputUnit.links;

        // Process intra unit links
        atomicLink.link.aUnit = inputUnit;
        atomicLink.link.bUnit = inputUnit;
        for (let i = 0, _i = unit.elements.length; i < _i; i++) {
            // add the current element
            builder.addToUnit(unit.id, unit.elements[i]);

            const srcIndex = SortedArray.indexOf(inputUnit.elements, unit.elements[i]);
            atomicLink.link.aIndex = srcIndex as StructureElement.UnitIndex;

            // check intra unit links
            for (let lI = intraLinkOffset[srcIndex], _lI = intraLinkOffset[srcIndex + 1]; lI < _lI; lI++) {
                atomicLink.link.bIndex = intraLinkB[lI] as StructureElement.UnitIndex;
                atomicLink.type = flags[lI];
                atomicLink.order = order[lI];
                if (linkTest(ctx)) {
                    builder.addToUnit(unit.id, inputUnit.elements[intraLinkB[lI]]);
                }
            }
        }

        // Process inter unit links
        for (const linkedUnit of interLinks.getLinkedUnits(inputUnit)) {
            if (processedUnits.has(linkedUnit.unitB.id)) continue;

            atomicLink.link.bUnit = linkedUnit.unitB;
            for (const aI of linkedUnit.linkedElementIndices) {
                // check if the element is in the expanded structure
                if (!SortedArray.has(unit.elements, inputUnit.elements[aI])) continue;

                atomicLink.link.aIndex = aI;
                for (const bond of linkedUnit.getBonds(aI)) {
                    atomicLink.link.bIndex = bond.indexB;
                    atomicLink.type = bond.flag;
                    atomicLink.order = bond.order;
                    if (linkTest(ctx)) {
                        builder.addToUnit(linkedUnit.unitB.id, linkedUnit.unitB.elements[bond.indexB]);
                    }
                }
            }
        }
    }

    return builder.getStructure();
}

// TODO: unionBy (skip this one?), cluster