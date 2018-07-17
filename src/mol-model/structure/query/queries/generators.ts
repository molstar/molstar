/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { StructureQuery } from '../query'
import { StructureSelection } from '../selection'
import { Unit, StructureProperties as P } from '../../structure'
import { Segmentation } from 'mol-data/int'
import { LinearGroupingBuilder } from '../utils/builders';
import { QueryPredicate, QueryFn, QueryContextView } from '../context';

export const all: StructureQuery = ctx => StructureSelection.Singletons(ctx.inputStructure, ctx.inputStructure);

export interface AtomsQueryParams {
    entityTest: QueryPredicate,
    chainTest: QueryPredicate,
    residueTest: QueryPredicate,
    atomTest: QueryPredicate,
    groupBy: QueryFn
}

export function residues(params?: Partial<AtomsQueryParams>) { return atoms({ ...params, groupBy: ctx => P.residue.key(ctx.element) }); }
export function chains(params?: Partial<AtomsQueryParams>) { return atoms({ ...params, groupBy: ctx => P.chain.key(ctx.element) }); }

function _true(ctx: QueryContextView) { return true; }
function _zero(ctx: QueryContextView) { return 0; }

export function atoms(params?: Partial<AtomsQueryParams>): StructureQuery {
    if (!params || (!params.atomTest && !params.residueTest && !params.chainTest && !params.entityTest && !params.groupBy)) return all;
    if (!!params.atomTest && !params.residueTest && !params.chainTest && !params.entityTest && !params.groupBy) return atomGroupsLinear(params.atomTest);

    const normalized: AtomsQueryParams = {
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
    return ctx => {
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
        }
        ctx.popCurrentElement();
        return StructureSelection.Singletons(inputStructure, builder.getStructure());
    };
}

function atomGroupsSegmented({ entityTest, chainTest, residueTest, atomTest }: AtomsQueryParams): StructureQuery {
    return ctx => {
        const { inputStructure } = ctx;
        const { units } = inputStructure;
        const l = ctx.pushCurrentElement();
        const builder = inputStructure.subsetBuilder(true);

        for (const unit of units) {
            if (unit.kind !== Unit.Kind.Atomic) continue;

            l.unit = unit;
            const elements = unit.elements;

            builder.beginUnit(unit.id);
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
                        if (atomTest(ctx)) {
                            builder.addElement(l.element);
                        }
                    }
                }
            }
            builder.commitUnit();
        }
        ctx.popCurrentElement();
        return StructureSelection.Singletons(inputStructure, builder.getStructure());
    };
}

function atomGroupsGrouped({ entityTest, chainTest, residueTest, atomTest, groupBy }: AtomsQueryParams): StructureQuery {
    return ctx => {
        const { inputStructure } = ctx;
        const { units } = inputStructure;
        const l = ctx.pushCurrentElement();
        const builder = new LinearGroupingBuilder(inputStructure);

        for (const unit of units) {
            if (unit.kind !== Unit.Kind.Atomic) continue;

            l.unit = unit;
            const elements = unit.elements;

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
                        if (atomTest(ctx)) builder.add(groupBy(ctx), unit.id, l.element);
                    }
                }
            }
        }
        ctx.popCurrentElement();
        return builder.getSelection();
    };
}