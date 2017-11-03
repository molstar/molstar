/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import Query from './query'
import Selection from './selection'
import * as P from './properties'
import { Structure, AtomSet, Atom } from '../structure'
import { OrderedSet, Segmentation } from 'mol-base/collections/integer'

export interface AtomGroupsSpec {
    entityTest: Atom.Predicate,
    chainTest: Atom.Predicate,
    residueTest: Atom.Predicate,
    atomTest: Atom.Predicate,
    groupBy: Atom.Property<any>
}

export const all: Query = s => s;

export function atomGroups(spec?: Partial<AtomGroupsSpec>): Query {
    if (!spec || (!spec.atomTest && !spec.residueTest && !spec.chainTest && !spec.entityTest && !spec.groupBy)) return all;
    if (!!spec.atomTest && !spec.residueTest && !spec.chainTest && !spec.entityTest && !spec.groupBy) return atomGroupsLinear(spec.atomTest);

    const normalized: AtomGroupsSpec = {
        entityTest: spec.entityTest || P.constant.true,
        chainTest: spec.entityTest || P.constant.true,
        residueTest: spec.residueTest || P.constant.true,
        atomTest: spec.atomTest || P.constant.true,
        groupBy: spec.entityTest || P.constant.zero,
    };

    if (!spec.groupBy) return atomGroupsSegmented(normalized)
    return atomGroupsGrouped(normalized);
}

function atomGroupsLinear(atomTest: Atom.Predicate): Query {
    return structure => {
        const { atoms, units } = structure;
        const unitIds = AtomSet.unitIds(atoms);
        const l = Atom.Location();
        const builder = AtomSet.LinearBuilder(atoms);

        for (let i = 0, _i = unitIds.length; i < _i; i++) {
            const unitId = unitIds[i];
            l.unit = units[unitId];
            const set = AtomSet.unitGetByIndex(atoms, i);

            builder.beginUnit();
            for (let j = 0, _j = OrderedSet.size(set); j < _j; j++) {
                l.atom = OrderedSet.getAt(set, j);
                if (atomTest(l)) builder.addToUnit(l.atom);
            }
            builder.commitUnit(unitId);
        }

        return Structure.create(units, builder.getSet());
    };
}

function atomGroupsSegmented({ entityTest, chainTest, residueTest, atomTest }: AtomGroupsSpec): Query {
    return structure => {
        const { atoms, units } = structure;
        const unitIds = AtomSet.unitIds(atoms);
        const l = Atom.Location();
        const builder = AtomSet.LinearBuilder(atoms);

        for (let i = 0, _i = unitIds.length; i < _i; i++) {
            const unitId = unitIds[i];
            const unit = units[unitId];
            l.unit = unit;
            const set = AtomSet.unitGetByIndex(atoms, i);

            builder.beginUnit();
            const chainsIt = Segmentation.transientSegments(unit.hierarchy.chainSegments, set);
            const residuesIt = Segmentation.transientSegments(unit.hierarchy.residueSegments, set);
            while (chainsIt.hasNext) {
                const chainSegment = chainsIt.move();
                l.atom = OrderedSet.getAt(set, chainSegment.start);
                // test entity and chain
                if (!entityTest(l) || !chainTest(l)) continue;

                residuesIt.setSegment(chainSegment);
                while (residuesIt.hasNext) {
                    const residueSegment = residuesIt.move();
                    l.atom = OrderedSet.getAt(set, residueSegment.start);

                    // test residue
                    if (!residueTest(l)) continue;

                    for (let j = residueSegment.start, _j = residueSegment.end; j < _j; j++) {
                        l.atom = OrderedSet.getAt(set, j);
                        if (atomTest(l)) builder.addToUnit(l.atom);
                    }
                }
            }
            builder.commitUnit(unitId);
        }

        return Structure.create(units, builder.getSet());
    };
}

class LinearGroupingBuilder {
    private builders: AtomSet.Builder[] = [];
    private builderMap: { [key: string]: AtomSet.Builder } = Object.create(null);

    add(key: any, unit: number, atom: number) {
        let b = this.builderMap[key];
        if (!b) {
            b = AtomSet.LinearBuilder(this.structure.atoms);
            this.builders[this.builders.length] = b;
            this.builderMap[key] = b;
        }
        b.add(unit, atom);
    }

    private allSingletons() {
        for (let i = 0, _i = this.builders.length; i < _i; i++) {
            if (this.builders[i].atomCount > 1) return false;
        }
        return true;
    }

    private singletonStructure(): Structure {
        const atoms: Atom[] = Atom.createEmptyArray(this.builders.length);
        for (let i = 0, _i = this.builders.length; i < _i; i++) {
            atoms[i] = this.builders[i].singleton();
        }
        return Structure.create(this.structure.units, AtomSet.create(atoms));
    }

    private fullSelection() {
        const ret: Structure[] = [];
        for (let i = 0, _i = this.builders.length; i < _i; i++) {
            ret[i] = Structure.create(this.structure.units, this.builders[i].getSet());
        }
        return ret;
    }

    getSelection(): Selection {
        const len = this.builders.length;
        if (len === 0) return Selection.Empty;
        if (len === 1) return Structure.create(this.structure.units, this.builders[0].getSet());
        if (this.allSingletons()) return this.singletonStructure();
        return this.fullSelection();
    }

    constructor(private structure: Structure) { }
}

function atomGroupsGrouped({ entityTest, chainTest, residueTest, atomTest, groupBy }: AtomGroupsSpec): Query {
    return structure => {
        throw 'nyi'
    };
}