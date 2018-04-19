/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import Query from './query'
import Selection from './selection'
import P from './properties'
import { Structure, ElementSet, Element } from '../structure'
import { OrderedSet, Segmentation } from 'mol-data/int'

export const all: Query.Provider = async (s, ctx) => Selection.Singletons(s, s.elements);

export interface AtomQueryParams {
    entityTest: Element.Predicate,
    chainTest: Element.Predicate,
    residueTest: Element.Predicate,
    atomTest: Element.Predicate,
    groupBy: Element.Property<any>
}

export interface AtomGroupsQueryParams extends AtomQueryParams {
    groupBy: Element.Property<any>
}

export function residues(params?: Partial<AtomQueryParams>) { return atoms({ ...params, groupBy: P.residue.key }); }
export function chains(params?: Partial<AtomQueryParams>) { return atoms({ ...params, groupBy: P.chain.key }); }

export function atoms(params?: Partial<AtomGroupsQueryParams>): Query.Provider {
    if (!params || (!params.atomTest && !params.residueTest && !params.chainTest && !params.entityTest && !params.groupBy)) return all;
    if (!!params.atomTest && !params.residueTest && !params.chainTest && !params.entityTest && !params.groupBy) return atomGroupsLinear(params.atomTest);

    const normalized: AtomGroupsQueryParams = {
        entityTest: params.entityTest || P.constant.true,
        chainTest: params.chainTest || P.constant.true,
        residueTest: params.residueTest || P.constant.true,
        atomTest: params.atomTest || P.constant.true,
        groupBy: params.groupBy || P.constant.zero,
    };

    if (!params.groupBy) return atomGroupsSegmented(normalized)
    return atomGroupsGrouped(normalized);
}

function atomGroupsLinear(atomTest: Element.Predicate): Query.Provider {
    return async (structure, ctx) => {
        const { elements, units } = structure;
        const unitIds = ElementSet.unitIds(elements);
        const l = Element.Location();
        const builder = ElementSet.LinearBuilder(elements);

        for (let i = 0, _i = unitIds.length; i < _i; i++) {
            const unitId = unitIds[i];
            l.unit = units[unitId];
            const set = ElementSet.unitGetByIndex(elements, i).elements;

            builder.beginUnit();
            for (let j = 0, _j = OrderedSet.size(set); j < _j; j++) {
                l.element = OrderedSet.getAt(set, j);
                if (atomTest(l)) builder.addToUnit(l.element);
            }
            builder.commitUnit(unitId);

            if (ctx.shouldUpdate) await ctx.update({ message: 'Atom Groups', current: 0, max: unitIds.length });
        }

        return Selection.Singletons(structure, builder.getSet());
    };
}

function atomGroupsSegmented({ entityTest, chainTest, residueTest, atomTest }: AtomGroupsQueryParams): Query.Provider {
    return async (structure, ctx) => {
        const { elements, units } = structure;
        const unitIds = ElementSet.unitIds(elements);
        const l = Element.Location();
        const builder = ElementSet.LinearBuilder(elements);

        for (let i = 0, _i = unitIds.length; i < _i; i++) {
            const unitId = unitIds[i];
            const unit = units[unitId];
            l.unit = unit;
            const set = ElementSet.unitGetByIndex(elements, i).elements;

            builder.beginUnit();
            const chainsIt = Segmentation.transientSegments(unit.hierarchy.chainSegments, set);
            const residuesIt = Segmentation.transientSegments(unit.hierarchy.residueSegments, set);
            while (chainsIt.hasNext) {
                const chainSegment = chainsIt.move();
                l.element = OrderedSet.getAt(set, chainSegment.start);
                // test entity and chain
                if (!entityTest(l) || !chainTest(l)) continue;

                residuesIt.setSegment(chainSegment);
                while (residuesIt.hasNext) {
                    const residueSegment = residuesIt.move();
                    l.element = OrderedSet.getAt(set, residueSegment.start);

                    // test residue
                    if (!residueTest(l)) continue;

                    for (let j = residueSegment.start, _j = residueSegment.end; j < _j; j++) {
                        l.element = OrderedSet.getAt(set, j);
                        if (atomTest(l)) builder.addToUnit(l.element);
                    }
                }
            }
            builder.commitUnit(unitId);

            if (ctx.shouldUpdate) await ctx.update({ message: 'Atom Groups', current: 0, max: unitIds.length });
        }

        return Selection.Singletons(structure, builder.getSet());
    };
}

class LinearGroupingBuilder {
    private builders: ElementSet.Builder[] = [];
    private builderMap = new Map<string, ElementSet.Builder>();

    add(key: any, unit: number, atom: number) {
        let b = this.builderMap.get(key);
        if (!b) {
            b = ElementSet.LinearBuilder(this.structure.elements);
            this.builders[this.builders.length] = b;
            this.builderMap.set(key, b);
        }
        b.add(unit, atom);
    }

    private allSingletons() {
        for (let i = 0, _i = this.builders.length; i < _i; i++) {
            if (this.builders[i].elementCount > 1) return false;
        }
        return true;
    }

    private singletonSelection(): Selection {
        const atoms: Element[] = Element.createEmptyArray(this.builders.length);
        for (let i = 0, _i = this.builders.length; i < _i; i++) {
            atoms[i] = this.builders[i].singleton();
        }
        return Selection.Singletons(this.structure, ElementSet.ofAtoms(atoms, this.structure.elements));
    }

    private fullSelection() {
        const sets: ElementSet[] = new Array(this.builders.length);
        for (let i = 0, _i = this.builders.length; i < _i; i++) {
            sets[i] = this.builders[i].getSet();
        }
        return Selection.Sequence(this.structure, sets);
    }

    getSelection(): Selection {
        const len = this.builders.length;
        if (len === 0) return Selection.Empty(this.structure);
        if (this.allSingletons()) return this.singletonSelection();
        return this.fullSelection();
    }

    constructor(private structure: Structure) { }
}

function atomGroupsGrouped({ entityTest, chainTest, residueTest, atomTest, groupBy }: AtomGroupsQueryParams): Query.Provider {
    return async (structure, ctx) => {
        const { elements, units } = structure;
        const unitIds = ElementSet.unitIds(elements);
        const l = Element.Location();
        const builder = new LinearGroupingBuilder(structure);

        for (let i = 0, _i = unitIds.length; i < _i; i++) {
            const unitId = unitIds[i];
            const unit = units[unitId];
            l.unit = unit;
            const set = ElementSet.unitGetByIndex(elements, i).elements;

            const chainsIt = Segmentation.transientSegments(unit.hierarchy.chainSegments, set);
            const residuesIt = Segmentation.transientSegments(unit.hierarchy.residueSegments, set);
            while (chainsIt.hasNext) {
                const chainSegment = chainsIt.move();
                l.element = OrderedSet.getAt(set, chainSegment.start);
                // test entity and chain
                if (!entityTest(l) || !chainTest(l)) continue;

                residuesIt.setSegment(chainSegment);
                while (residuesIt.hasNext) {
                    const residueSegment = residuesIt.move();
                    l.element = OrderedSet.getAt(set, residueSegment.start);

                    // test residue
                    if (!residueTest(l)) continue;

                    for (let j = residueSegment.start, _j = residueSegment.end; j < _j; j++) {
                        l.element = OrderedSet.getAt(set, j);
                        if (atomTest(l)) builder.add(groupBy(l), unitId, l.element);
                    }
                }
            }

            if (ctx.shouldUpdate) await ctx.update({ message: 'Atom Groups', current: 0, max: unitIds.length });
        }

        return builder.getSelection();
    };
}