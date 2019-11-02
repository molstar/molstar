/**
 * Copyright (c) 2017-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { StructureQuery } from '../query'
import { StructureSelection } from '../selection'
import { Unit, StructureProperties as P, StructureElement } from '../../structure'
import { Segmentation, SortedArray } from '../../../../mol-data/int'
import { LinearGroupingBuilder } from '../utils/builders';
import { QueryPredicate, QueryFn, QueryContextView } from '../context';
import { UnitRing } from '../../structure/unit/rings';
import Structure from '../../structure/structure';
import { ElementIndex } from '../../model';
import { UniqueArray } from '../../../../mol-data/generic';
import { structureSubtract } from '../utils/structure-set';

export const none: StructureQuery = ctx => StructureSelection.Sequence(ctx.inputStructure, []);
export const all: StructureQuery = ctx => StructureSelection.Singletons(ctx.inputStructure, ctx.inputStructure);

export interface AtomsQueryParams {
    /** Query to be executed for each unit once */
    unitTest: QueryPredicate,
    /** Query to be executed for each entity once */
    entityTest: QueryPredicate,
    /** Query to be executed for each chain once */
    chainTest: QueryPredicate,
    /** Query to be executed for each residue (or coarse element) once */
    residueTest: QueryPredicate,
    /** Query to be executed for each atom */
    atomTest: QueryPredicate,
    groupBy: QueryFn
}

export function residues(params?: Partial<AtomsQueryParams>) { return atoms({ ...params, groupBy: ctx => P.residue.key(ctx.element) }); }
export function chains(params?: Partial<AtomsQueryParams>) { return atoms({ ...params, groupBy: ctx => P.chain.key(ctx.element) }); }

function _true(ctx: QueryContextView) { return true; }
function _zero(ctx: QueryContextView) { return 0; }

export function atoms(params?: Partial<AtomsQueryParams>): StructureQuery {
    if (!params || (!params.atomTest && !params.residueTest && !params.chainTest && !params.entityTest && !params.unitTest && !params.groupBy)) return all;
    if (!!params.atomTest && !params.residueTest && !params.chainTest && !params.entityTest && !params.unitTest && !params.groupBy) return atomGroupsLinear(params.atomTest);

    const normalized: AtomsQueryParams = {
        unitTest: params.unitTest || _true,
        entityTest: params.entityTest || _true,
        chainTest: params.chainTest || _true,
        residueTest: params.residueTest || _true,
        atomTest: params.atomTest || _true,
        groupBy: params.groupBy || _zero,
    };

    if (!params.groupBy) return atomGroupsSegmented(normalized)
    return atomGroupsGrouped(normalized);
}

function atomGroupsLinear(atomTest: QueryPredicate): StructureQuery {
    return function query_atomGroupsLinear(ctx) {
        const { inputStructure } = ctx;
        const { units } = inputStructure;
        const l = ctx.pushCurrentElement();
        const builder = inputStructure.subsetBuilder(true);

        for (const unit of units) {
            l.unit = unit;
            const elements = unit.elements;

            builder.beginUnit(unit.id);
            for (let j = 0, _j = elements.length; j < _j; j++) {
                l.element = elements[j];
                if (atomTest(ctx)) builder.addElement(l.element);
            }
            builder.commitUnit();

            ctx.throwIfTimedOut();
        }
        ctx.popCurrentElement();
        return StructureSelection.Singletons(inputStructure, builder.getStructure());
    };
}

function atomGroupsSegmented({ unitTest, entityTest, chainTest, residueTest, atomTest }: AtomsQueryParams): StructureQuery {
    return function query_atomGroupsSegmented(ctx) {
        const { inputStructure } = ctx;
        const { units } = inputStructure;
        const l = ctx.pushCurrentElement();
        const builder = inputStructure.subsetBuilder(true);

        for (const unit of units) {
            l.unit = unit;
            if (!unitTest(ctx)) continue;

            const { elements, model } = unit;
            builder.beginUnit(unit.id);

            if (unit.kind === Unit.Kind.Atomic) {
                const chainsIt = Segmentation.transientSegments(unit.model.atomicHierarchy.chainAtomSegments, elements);
                const residuesIt = Segmentation.transientSegments(unit.model.atomicHierarchy.residueAtomSegments, elements);

                while (chainsIt.hasNext) {
                    const chainSegment = chainsIt.move();
                    l.element = elements[chainSegment.start];
                    // test entity and chain
                    if (!entityTest(ctx) || !chainTest(ctx)) continue;

                    residuesIt.setSegment(chainSegment);
                    while (residuesIt.hasNext) {
                        const residueSegment = residuesIt.move();
                        l.element = elements[residueSegment.start];

                        // test residue
                        if (!residueTest(ctx)) continue;

                        for (let j = residueSegment.start, _j = residueSegment.end; j < _j; j++) {
                            l.element = elements[j];
                            // test atom
                            if (atomTest(ctx)) {
                                builder.addElement(l.element);
                            }
                        }
                    }
                }
            } else {
                const { chainElementSegments } = Unit.Kind.Spheres ? model.coarseHierarchy.spheres : model.coarseHierarchy.gaussians
                const chainsIt = Segmentation.transientSegments(chainElementSegments, elements);

                while (chainsIt.hasNext) {
                    const chainSegment = chainsIt.move();
                    l.element = elements[chainSegment.start];
                    // test entity and chain
                    if (!entityTest(ctx) || !chainTest(ctx)) continue;

                    for (let j = chainSegment.start, _j = chainSegment.end; j < _j; j++) {
                        l.element = elements[j];
                        // test residue/coarse element
                        if (residueTest(ctx)) {
                            builder.addElement(l.element);
                        }
                    }
                }
            }
            builder.commitUnit();

            ctx.throwIfTimedOut();
        }
        ctx.popCurrentElement();
        return StructureSelection.Singletons(inputStructure, builder.getStructure());
    };
}

function atomGroupsGrouped({ unitTest, entityTest, chainTest, residueTest, atomTest, groupBy }: AtomsQueryParams): StructureQuery {
    return function query_atomGroupsGrouped(ctx) {
        const { inputStructure } = ctx;
        const { units } = inputStructure;
        const l = ctx.pushCurrentElement();
        const builder = new LinearGroupingBuilder(inputStructure);

        for (const unit of units) {
            l.unit = unit;
            if (!unitTest(ctx)) continue;

            const { elements, model } = unit;

            if (unit.kind === Unit.Kind.Atomic) {
                const chainsIt = Segmentation.transientSegments(model.atomicHierarchy.chainAtomSegments, elements);
                const residuesIt = Segmentation.transientSegments(model.atomicHierarchy.residueAtomSegments, elements);

                while (chainsIt.hasNext) {
                    const chainSegment = chainsIt.move();
                    l.element = elements[chainSegment.start];
                    // test entity and chain
                    if (!entityTest(ctx) || !chainTest(ctx)) continue;

                    residuesIt.setSegment(chainSegment);
                    while (residuesIt.hasNext) {
                        const residueSegment = residuesIt.move();
                        l.element = elements[residueSegment.start];

                        // test residue
                        if (!residueTest(ctx)) continue;

                        for (let j = residueSegment.start, _j = residueSegment.end; j < _j; j++) {
                            l.element = elements[j];
                            // test atom
                            if (atomTest(ctx)) {
                                builder.add(groupBy(ctx), unit.id, l.element);
                            }
                        }
                    }
                }
            } else {
                const { chainElementSegments } = Unit.Kind.Spheres ? model.coarseHierarchy.spheres : model.coarseHierarchy.gaussians
                const chainsIt = Segmentation.transientSegments(chainElementSegments, elements);
                while (chainsIt.hasNext) {
                    const chainSegment = chainsIt.move();
                    l.element = elements[chainSegment.start];
                    // test entity and chain
                    if (!entityTest(ctx) || !chainTest(ctx)) continue;

                    for (let j = chainSegment.start, _j = chainSegment.end; j < _j; j++) {
                        l.element = elements[j];
                        // test residue/coarse element
                        if (residueTest(ctx)) {
                            builder.add(groupBy(ctx), unit.id, l.element);
                        }
                    }
                }
            }

            ctx.throwIfTimedOut();
        }
        ctx.popCurrentElement();
        return builder.getSelection();
    };
}

function getRingStructure(unit: Unit.Atomic, ring: UnitRing, inputStructure: Structure) {
    const elements = new Int32Array(ring.length) as any as ElementIndex[];
    for (let i = 0, _i = ring.length; i < _i; i++) elements[i] = unit.elements[ring[i]];
    return Structure.create([unit.getChild(SortedArray.ofSortedArray(elements))], { parent: inputStructure });
}

export function rings(fingerprints?: ArrayLike<UnitRing.Fingerprint>): StructureQuery {
    return function query_rings(ctx) {
        const { units } = ctx.inputStructure;
        const ret = StructureSelection.LinearBuilder(ctx.inputStructure);

        if (!fingerprints || fingerprints.length === 0) {
            for (const u of units) {
                if (!Unit.isAtomic(u)) continue;

                for (const r of u.rings.all) {
                    ret.add(getRingStructure(u, r, ctx.inputStructure));
                }
            }
        } else {
            const uniqueFps = UniqueArray.create<UnitRing.Fingerprint, UnitRing.Fingerprint>();
            for (let i = 0; i < fingerprints.length; i++) UniqueArray.add(uniqueFps, fingerprints[i], fingerprints[i]);

            for (const u of units) {
                if (!Unit.isAtomic(u)) continue;

                const rings = u.rings;
                for (const fp of uniqueFps.array) {
                    if (!rings.byFingerprint.has(fp)) continue;
                    for (const r of rings.byFingerprint.get(fp)!) {
                        ret.add(getRingStructure(u, rings.all[r], ctx.inputStructure));
                    }
                }
            }
        }

        return ret.getSelection();
    }
}

export function querySelection(selection: StructureQuery, query: StructureQuery, inComplement: boolean = false): StructureQuery {
    return function query_querySelection(ctx) {
        const targetSel = selection(ctx);
        if (StructureSelection.structureCount(targetSel) === 0) return targetSel;

        const target = inComplement
            ? structureSubtract(ctx.inputStructure, StructureSelection.unionStructure(targetSel))
            : StructureSelection.unionStructure(targetSel);

        if (target.elementCount === 0) return StructureSelection.Empty(ctx.inputStructure);
        ctx.throwIfTimedOut();

        ctx.pushInputStructure(target);
        const result = query(ctx);
        ctx.popInputStructure();
        return StructureSelection.withInputStructure(result, ctx.inputStructure);
    }
}

export function linkedAtomicPairs(linkTest?: QueryPredicate): StructureQuery {
    return function query_linkedAtomicPairs(ctx) {
        const structure = ctx.inputStructure;

        const interLinks = structure.links;
        // Note: each link is called twice, that's why we need the unique builder.
        const ret = StructureSelection.UniqueBuilder(ctx.inputStructure);

        ctx.pushCurrentLink();
        const atomicLink = ctx.atomicLink;
        atomicLink.setTestFn(linkTest);

        // Process intra unit links
        for (const unit of structure.units) {
            if (unit.kind !== Unit.Kind.Atomic) continue;

            const { offset: intraLinkOffset, b: intraLinkB, edgeProps: { flags, order } } = unit.links;
            atomicLink.a.unit = unit;
            atomicLink.b.unit = unit;
            for (let i = 0 as StructureElement.UnitIndex, _i = unit.elements.length; i < _i; i++) {
                atomicLink.aIndex = i as StructureElement.UnitIndex;
                atomicLink.a.element = unit.elements[i];

                // check intra unit links
                for (let lI = intraLinkOffset[i], _lI = intraLinkOffset[i + 1]; lI < _lI; lI++) {
                    atomicLink.bIndex = intraLinkB[lI] as StructureElement.UnitIndex;
                    atomicLink.b.element = unit.elements[intraLinkB[lI]];
                    atomicLink.type = flags[lI];
                    atomicLink.order = order[lI];
                    // No need to "swap test" because each bond direction will be visited eventually.
                    if (atomicLink.test(ctx, false)) {
                        const b = structure.subsetBuilder(false);
                        b.beginUnit(unit.id);
                        b.addElement(atomicLink.a.element);
                        b.addElement(atomicLink.b.element);
                        b.commitUnit();
                        ret.add(b.getStructure());
                    }
                }
            }
        }

        // Process inter unit links
        for (const bond of interLinks.bonds) {
            atomicLink.a.unit = bond.unitA;
            atomicLink.a.element = bond.unitA.elements[bond.indexA];
            atomicLink.aIndex = bond.indexA;
            atomicLink.b.unit = bond.unitB;
            atomicLink.b.element = bond.unitB.elements[bond.indexB];
            atomicLink.bIndex = bond.indexB;
            atomicLink.order = bond.order;
            atomicLink.type = bond.flag;

            // No need to "swap test" because each bond direction will be visited eventually.
            if (atomicLink.test(ctx, false)) {
                const b = structure.subsetBuilder(false);
                b.addToUnit(atomicLink.a.unit.id, atomicLink.a.element);
                b.addToUnit(atomicLink.b.unit.id, atomicLink.b.element);
                ret.add(b.getStructure());
            }
        }

        ctx.popCurrentLink();
        return ret.getSelection();
    };
}