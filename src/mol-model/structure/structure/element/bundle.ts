/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { OrderedSet, SortedArray, Interval } from '../../../../mol-data/int';
import { ElementIndex } from '../../model';
import Structure from '../structure';
import Unit from '../unit';
import { hashFnv32a, hash2 } from '../../../../mol-data/util';
import SortedRanges from '../../../../mol-data/int/sorted-ranges';
import { UnitIndex } from './element';
import { Loci } from './loci';
import Expression from '../../../../mol-script/language/expression';
import { MolScriptBuilder as MS } from '../../../../mol-script/language/builder';
import { StructureSelection } from '../../query';

export interface BundleElement {
    /**
     * Array (sorted by first element in sub-array) of
     * arrays of `Unit.id`s that share the same `Unit.invariantId`
     */
    groupedUnits: SortedArray<number>[],
    set: SortedArray<UnitIndex>
    ranges: SortedRanges<UnitIndex>
}

export interface Bundle {
    /** Hash of the structure with which the bundle is compatible */
    readonly hash: number
    /** Bundle elements */
    readonly elements: ReadonlyArray<Readonly<BundleElement>>
}

export namespace Bundle {
    export const Empty: Bundle = { hash: -1, elements: [] };

    export function fromSubStructure(parent: Structure, structure: Structure) {
        return fromLoci(StructureSelection.toLociWithSourceUnits(StructureSelection.Singletons(parent, structure)));
    }

    export function fromSelection(selection: StructureSelection) {
        return fromLoci(StructureSelection.toLociWithSourceUnits(selection));
    }

    export function fromLoci(loci: Loci): Bundle {
        const _elements: {
            unit: Unit
            set: SortedArray<UnitIndex>
            ranges: SortedRanges<UnitIndex>
        }[] = [];
        for (const e of loci.elements) {
            const { unit, indices } = e;
            if (OrderedSet.size(indices) === 0) continue;

            const ranges: UnitIndex[] = [];
            const set: UnitIndex[] = [];

            if (OrderedSet.isInterval(indices)) {
                if (OrderedSet.size(indices) === 1) {
                    set.push(Interval.min(indices));
                } else {
                    ranges.push(Interval.min(indices), Interval.max(indices));
                }
            } else {
                let i = 0, len = indices.length;
                while (i < len) {
                    const start = i;
                    i++;
                    while (i < len && indices[i - 1] + 1 === indices[i]) i++;
                    const end = i;
                    if (end - start > 2) {
                        ranges.push(indices[start], indices[end - 1]);
                    } else {
                        for (let j = start; j < end; j++) {
                            set[set.length] = indices[j];
                        }
                    }
                }
            }

            _elements.push({
                unit,
                set: SortedArray.ofSortedArray(set),
                ranges: SortedRanges.ofSortedRanges(ranges)
            });
        }

        const elementGroups = new Map<number, {
            groupedUnits: Map<number, number[]>
            set: SortedArray<UnitIndex>
            ranges: SortedRanges<UnitIndex>
        }>();
        for (let i = 0, il = _elements.length; i < il; ++i) {
            const e = _elements[i];
            const key = hash2(hashFnv32a(e.ranges), hashFnv32a(e.set));
            if (elementGroups.has(key)) {
                const { groupedUnits } = elementGroups.get(key)!;
                if (groupedUnits.has(e.unit.invariantId)) {
                    groupedUnits.get(e.unit.invariantId)!.push(e.unit.id);
                } else {
                    groupedUnits.set(e.unit.invariantId, [e.unit.id]);
                }
            } else {
                const groupedUnits = new Map<number, number[]>();
                groupedUnits.set(e.unit.invariantId, [e.unit.id]);
                elementGroups.set(key, { groupedUnits, set: e.set, ranges: e.ranges });
            }
        }

        const elements: BundleElement[] = [];
        elementGroups.forEach(e => {
            const groupedUnits: SortedArray<number>[] = [];
            e.groupedUnits.forEach(g => groupedUnits.push(SortedArray.ofUnsortedArray(g)));
            groupedUnits.sort((a, b) => a[0] - b[0]); // sort by first unit id of each group
            elements.push({ groupedUnits, set: e.set, ranges: e.ranges });
        });

        return { hash: loci.structure.hashCode, elements };
    }

    function getUnitsFromIds(unitIds: ArrayLike<number>, structure: Structure) {
        const units: Unit[] = [];
        for (let i = 0, il = unitIds.length; i < il; ++i) {
            const unitId = unitIds[i];
            if (structure.unitMap.has(unitId)) units.push(structure.unitMap.get(unitId));
        }
        return units;
    }

    export function toLoci(bundle: Bundle, structure: Structure): Loci {
        if (bundle.hash !== -1 && bundle.hash !== structure.hashCode) {
            new Error('Bundle not compatible with given structure');
        }
        const elements: Loci['elements'][0][] = [];
        for (const e of bundle.elements) {
            for (const g of e.groupedUnits) {
                const units = getUnitsFromIds(g, structure);
                if (units.length === 0) continue;

                let indices: OrderedSet<UnitIndex>;
                if (e.ranges.length === 0) {
                    indices = e.set;
                } else if (e.set.length === 0) {
                    if (e.ranges.length === 2) {
                        indices = Interval.ofRange(e.ranges[0], e.ranges[1]);
                    } else {
                        const _indices = new Int32Array(SortedRanges.size(e.ranges));
                        SortedRanges.forEach(e.ranges, (v, i) => _indices[i] = v);
                        indices = SortedArray.ofSortedArray(_indices);
                    }
                } else {
                    const rangesSize = SortedRanges.size(e.ranges);
                    const _indices = new Int32Array(e.set.length + rangesSize);
                    SortedRanges.forEach(e.ranges, (v, i) => _indices[i] = v);
                    _indices.set(e.set, rangesSize);
                    indices = SortedArray.ofUnsortedArray(_indices); // requires sort
                }

                for (const unit of units) {
                    elements.push({ unit, indices });
                }
            }
        }
        return Loci(structure, elements);
    }

    export function toStructure(bundle: Bundle, parent: Structure): Structure {
        if (bundle.hash !== -1 && bundle.hash !== parent.hashCode) {
            new Error('Bundle not compatible with given structure');
        }
        const units: Unit[] = [];
        for (const e of bundle.elements) {
            for (const g of e.groupedUnits) {
                const _units = getUnitsFromIds(g, parent);
                if (_units.length === 0) continue;

                const ue = _units[0].elements; // the elements are grouped by unit.invariantId
                const rangesSize = SortedRanges.size(e.ranges);
                const setSize = e.set.length;
                const _indices = new Int32Array(setSize + rangesSize);
                let indices: SortedArray<ElementIndex>;
                if (rangesSize === 0) {
                    for (let i = 0, il = setSize; i < il; ++i) {
                        _indices[i] = ue[e.set[i]];
                    }
                    indices = SortedArray.ofSortedArray(_indices);
                } else if (setSize === 0) {
                    SortedRanges.forEach(e.ranges, (v, i) => _indices[i] = ue[v]);
                    indices = SortedArray.ofSortedArray(_indices);
                } else {
                    if (SortedArray.min(e.set) > SortedRanges.max(e.ranges)) {
                        SortedRanges.forEach(e.ranges, (v, i) => _indices[i] = ue[v]);
                        for (let i = 0, il = setSize; i < il; ++i) {
                            _indices[i + rangesSize] = ue[e.set[i]];
                        }
                        indices = SortedArray.ofSortedArray(_indices);
                    } else if (SortedRanges.min(e.ranges) > SortedArray.max(e.set)) {
                        for (let i = 0, il = setSize; i < il; ++i) {
                            _indices[i] = ue[e.set[i]];
                        }
                        SortedRanges.forEach(e.ranges, (v, i) => _indices[i + setSize] = ue[v]);
                        indices = SortedArray.ofSortedArray(_indices);
                    } else {
                        SortedRanges.forEach(e.ranges, (v, i) => _indices[i] = ue[v]);
                        for (let i = 0, il = setSize; i < il; ++i) {
                            _indices[i + rangesSize] = ue[e.set[i]];
                        }
                        indices = SortedArray.ofUnsortedArray(_indices); // requires sort
                    }
                }

                for (const unit of _units) {
                    units.push(unit.getChild(indices));
                }
            }
        }
        return Structure.create(units, { parent });
    }


    function elementToExpression(e: BundleElement): Expression {
        return MS.internal.generator.bundleElement({
            groupedUnits: MS.core.type.list(e.groupedUnits.map(u => MS.core.type.list(u))),
            ranges: MS.core.type.list(e.ranges),
            set: MS.core.type.list(e.set),
        });
    }

    export function toExpression(bundle: Bundle): Expression {
        return MS.internal.generator.bundle({
            elements: MS.core.type.list(bundle.elements.map(elementToExpression))
        });
    }

    export function areEqual(a: Bundle, b: Bundle) {
        if (a.elements.length !== b.elements.length) return false;
        for (let i = 0, il = a.elements.length; i < il; ++i) {
            const elementA = a.elements[i], elementB = b.elements[i];
            if (elementA.groupedUnits.length !== elementB.groupedUnits.length) return false;
            for (let j = 0, jl = elementB.groupedUnits.length; j < jl; ++j) {
                if (!SortedArray.areEqual(elementA.groupedUnits[j], elementB.groupedUnits[j])) return false;
            }
            if (!SortedArray.areEqual(elementA.set, elementB.set)) return false;
            if (!SortedRanges.areEqual(elementA.ranges, elementB.ranges)) return false;
        }
        return true;
    }
}