/**
 * Copyright (c) 2017-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { UniqueArray } from '../../../../mol-data/generic';
import { OrderedSet, SortedArray, Interval } from '../../../../mol-data/int';
import { BoundaryHelper } from '../../../../mol-math/geometry/boundary-helper';
import { Vec3 } from '../../../../mol-math/linear-algebra';
import { MolScriptBuilder as MS } from '../../../../mol-script/language/builder';
import Structure from '../structure';
import Unit from '../unit';
import { Boundary } from '../util/boundary';
import { sortArray, hashFnv32a, hash2 } from '../../../../mol-data/util';
import Expression from '../../../../mol-script/language/expression';
import { ElementIndex } from '../../model';
import { UnitIndex } from './element';
import { Location } from './location';
import { ChainIndex } from '../../model/indexing';

/** Represents multiple element index locations */
export interface Loci {
    readonly kind: 'element-loci',
    readonly structure: Structure,
    /** Access i-th element as unit.elements[indices[i]] */
    readonly elements: ReadonlyArray<{
        unit: Unit,
        /**
         * Indices into the unit.elements array.
         * Can use OrderedSet.forEach to iterate (or OrderedSet.size + OrderedSet.getAt)
         */
        indices: OrderedSet<UnitIndex>
    }>
}

export function Loci(structure: Structure, elements: ArrayLike<{ unit: Unit, indices: OrderedSet<UnitIndex> }>): Loci {
    return { kind: 'element-loci', structure, elements: elements as Loci['elements'] };
}

export namespace Loci {
    export function is(x: any): x is Loci {
        return !!x && x.kind === 'element-loci';
    }

    export function areEqual(a: Loci, b: Loci) {
        if (a.elements.length !== b.elements.length) return false
        for (let i = 0, il = a.elements.length; i < il; ++i) {
            const elementA = a.elements[i]
            const elementB = b.elements[i]
            if (elementA.unit.id !== elementB.unit.id) return false
            if (!OrderedSet.areEqual(elementA.indices, elementB.indices)) return false
        }
        return true
    }

    export function isEmpty(loci: Loci) {
        return size(loci) === 0
    }

    export function isWholeStructure(loci: Loci) {
        return size(loci) === loci.structure.elementCount
    }

    export function size(loci: Loci) {
        let s = 0;
        for (const u of loci.elements) s += OrderedSet.size(u.indices);
        return s;
    }

    export function all(structure: Structure): Loci {
        return Loci(structure, structure.units.map(unit => ({
            unit,
            indices: OrderedSet.ofBounds<UnitIndex>(0 as UnitIndex, unit.elements.length as UnitIndex)
        })));
    }

    export function none(structure: Structure): Loci {
        return Loci(structure, []);
    }

    export function getFirstLocation(loci: Loci, e?: Location): Location | undefined {
        if (isEmpty(loci)) return void 0;
        const unit = loci.elements[0].unit;
        const element = unit.elements[OrderedSet.getAt(loci.elements[0].indices, 0)];
        if (e) {
            e.unit = loci.elements[0].unit;
            e.element = element;
            return e;
        }
        return Location.create(unit, element);
    }

    export function toStructure(loci: Loci): Structure {
        const units: Unit[] = []
        for (const e of loci.elements) {
            const { unit, indices } = e
            const elements = new Int32Array(OrderedSet.size(indices))
            OrderedSet.forEach(indices, (v, i) => elements[i] = unit.elements[v])
            units.push(unit.getChild(SortedArray.ofSortedArray(elements)))
        }
        return Structure.create(units, { parent: loci.structure.parent })
    }

    export function remap(loci: Loci, structure: Structure): Loci {
        if (structure === loci.structure) return loci

        const elements: Loci['elements'][0][] = [];
        loci.elements.forEach(e => {
            if (!structure.unitMap.has(e.unit.id)) return
            const unit = structure.unitMap.get(e.unit.id)

            if (SortedArray.areEqual(e.unit.elements, unit.elements)) {
                elements.push({ unit, indices: e.indices })
            } else {
                const _indices: UnitIndex[] = []
                const end = unit.elements.length
                let start = 0
                for (let i = 0; i < OrderedSet.size(e.indices); ++i) {
                    const v = OrderedSet.getAt(e.indices, i)
                    const eI = e.unit.elements[v]
                    const uI = SortedArray.indexOfInRange(unit.elements, eI, start, end) as UnitIndex | -1
                    if (uI !== -1) {
                        _indices.push(uI)
                        start = uI
                    }
                }

                let indices: OrderedSet<UnitIndex>
                if (_indices.length > 12 && _indices[_indices.length - 1] - _indices[0] === _indices.length - 1) {
                    indices = Interval.ofRange(_indices[0], _indices[_indices.length - 1])
                } else {
                    indices = SortedArray.ofSortedArray(_indices)
                }

                if (OrderedSet.size(indices) > 0) elements.push({ unit, indices })
            }
        });

        return Loci(structure, elements);
    }

    /** Create union of `xs` and `ys` */
    export function union(xs: Loci, ys: Loci): Loci {
        if (xs.elements.length > ys.elements.length) return union(ys, xs);
        if (Loci.isEmpty(xs)) return ys;

        const map = new Map<number, OrderedSet<UnitIndex>>();

        for (const e of xs.elements) map.set(e.unit.id, e.indices);

        const elements: Loci['elements'][0][] = [];
        for (const e of ys.elements) {
            if (map.has(e.unit.id)) {
                elements[elements.length] = { unit: e.unit, indices: OrderedSet.union(map.get(e.unit.id)!, e.indices) };
                map.delete(e.unit.id)
            } else {
                elements[elements.length] = e;
            }
        }

        map.forEach((indices, id) => {
            elements[elements.length] = { unit: xs.structure.unitMap.get(id)!, indices };
        });

        return Loci(xs.structure, elements);
    }

    /** Subtract `ys` from `xs` */
    export function subtract(xs: Loci, ys: Loci): Loci {
        const map = new Map<number, OrderedSet<UnitIndex>>();
        for (const e of ys.elements) map.set(e.unit.id, e.indices);

        const elements: Loci['elements'][0][] = [];
        for (const e of xs.elements) {
            if (map.has(e.unit.id)) {
                const indices = OrderedSet.subtract(e.indices, map.get(e.unit.id)!);
                if (OrderedSet.size(indices) === 0) continue;
                elements[elements.length] = { unit: e.unit, indices };
            } else {
                elements[elements.length] = e;
            }
        }

        return Loci(xs.structure, elements);
    }

    export function areIntersecting(xs: Loci, ys: Loci): boolean {
        if (xs.elements.length > ys.elements.length) return areIntersecting(ys, xs);
        if (Loci.isEmpty(xs)) return Loci.isEmpty(ys);

        const map = new Map<number, OrderedSet<UnitIndex>>();

        for (const e of xs.elements) map.set(e.unit.id, e.indices);
        for (const e of ys.elements) {
            if (!map.has(e.unit.id)) continue;
            if (OrderedSet.areIntersecting(map.get(e.unit.id)!, e.indices)) return true;
        }

        return false;
    }

    /** Check if second loci is a subset of the first */
    export function isSubset(xs: Loci, ys: Loci): boolean {
        if (Loci.isEmpty(xs)) return Loci.isEmpty(ys);

        const map = new Map<number, OrderedSet<UnitIndex>>();

        for (const e of xs.elements) map.set(e.unit.id, e.indices);
        for (const e of ys.elements) {
            if (!map.has(e.unit.id)) continue;
            if (!OrderedSet.isSubset(map.get(e.unit.id)!, e.indices)) return false;
        }

        return true;
    }

    export function extendToWholeResidues(loci: Loci, restrictToConformation?: boolean): Loci {
        const elements: Loci['elements'][0][] = [];
        const residueAltIds = new Set<string>()

        for (const lociElement of loci.elements) {
            if (lociElement.unit.kind === Unit.Kind.Atomic) {
                const unitElements = lociElement.unit.elements;
                const h = lociElement.unit.model.atomicHierarchy;
                const { label_alt_id } = lociElement.unit.model.atomicHierarchy.atoms;

                const { index: residueIndex, offsets: residueOffsets } = h.residueAtomSegments;

                const newIndices: UnitIndex[] = [];
                const indices = lociElement.indices, len = OrderedSet.size(indices);
                let i = 0;
                while (i < len) {
                    residueAltIds.clear()
                    const eI = unitElements[OrderedSet.getAt(indices, i)]
                    const rI = residueIndex[eI];
                    residueAltIds.add(label_alt_id.value(eI))
                    i++;
                    while (i < len) {
                        const eI = unitElements[OrderedSet.getAt(indices, i)]
                        if (residueIndex[eI] !== rI) break;
                        residueAltIds.add(label_alt_id.value(eI))
                        i++;
                    }
                    const hasSharedAltId = residueAltIds.has('')
                    for (let j = residueOffsets[rI], _j = residueOffsets[rI + 1]; j < _j; j++) {
                        const idx = OrderedSet.indexOf(unitElements, j);
                        if (idx >= 0) {
                            const altId = label_alt_id.value(j)
                            if (!restrictToConformation || hasSharedAltId || !altId || residueAltIds.has(altId)) {
                                newIndices[newIndices.length] = idx as UnitIndex;
                            }
                        }
                    }
                }

                elements[elements.length] = { unit: lociElement.unit, indices: SortedArray.ofSortedArray(newIndices) };
            } else {
                // coarse elements are already by-residue
                elements[elements.length] = lociElement;
            }
        }

        return Loci(loci.structure, elements);
    }

    function getChainSegments(unit: Unit) {
        switch (unit.kind) {
            case Unit.Kind.Atomic: return unit.model.atomicHierarchy.chainAtomSegments
            case Unit.Kind.Spheres: return unit.model.coarseHierarchy.spheres.chainElementSegments
            case Unit.Kind.Gaussians: return unit.model.coarseHierarchy.gaussians.chainElementSegments
        }
    }

    function makeIndexSet(newIndices: number[]): OrderedSet<UnitIndex> {
        if (newIndices.length > 12 && newIndices[newIndices.length - 1] - newIndices[0] === newIndices.length - 1) {
            return Interval.ofRange(newIndices[0], newIndices[newIndices.length - 1])
        } else {
            return  SortedArray.ofSortedArray(newIndices)
        }
    }

    function collectChains(unit: Unit, chainIndices: Set<ChainIndex>, elements: Loci['elements'][0][]) {
        const { index } = getChainSegments(unit);
        const xs = unit.elements;
        const newIndices: UnitIndex[] = [];
        for (let i = 0 as UnitIndex, _i = xs.length; i < _i; i++) {
            const eI = xs[i];
            const cI = index[eI];
            if (!chainIndices.has(cI)) continue;
            newIndices[newIndices.length] = i;
        }

        if (newIndices.length > 0) {
            elements[elements.length] = { unit, indices: makeIndexSet(newIndices) };
        }
    }

    function extendGroupToWholeChains(loci: Loci, start: number, end: number, isPartitioned: boolean, elements: Loci['elements'][0][]) {
        const { index: chainIndex } = getChainSegments(loci.elements[0].unit);

        const chainIndices = new Set<ChainIndex>();

        for (let lI = start; lI < end; lI++) {
            const lociElement = loci.elements[lI];
            const indices = lociElement.indices;
            const unitElements = lociElement.unit.elements;
            for (let i = 0, _i = OrderedSet.size(indices); i < _i; i++) {
                chainIndices.add(chainIndex[unitElements[OrderedSet.getAt(indices, i)]]);
            }
        }

        if (isPartitioned) {
            const groupId = loci.elements[0].unit.chainGroupId, operator = loci.elements[0].unit.conformation.operator;
            // TODO: check for accidental quadratic for really large structures (but should be ok).
            for (const unit of loci.structure.units) {
                if (unit.chainGroupId !== groupId || unit.conformation.operator !== operator) continue;
                collectChains(unit, chainIndices, elements);
            }
        } else {
            for (let lI = start; lI < end; lI++) {
                collectChains(loci.elements[lI].unit, chainIndices, elements);
            }
        }
    }

    export function extendToWholeChains(loci: Loci): Loci {
        const elements: Loci['elements'][0][] = [];

        for (let i = 0, len = loci.elements.length; i < len; i++) {
            const e = loci.elements[i];
            if (Unit.Traits.is(e.unit.traits, Unit.Trait.Patitioned)) {
                const start = i;
                while (i < len
                    && loci.elements[i].unit.chainGroupId === e.unit.chainGroupId
                    && loci.elements[i].unit.conformation.operator === e.unit.conformation.operator) {
                    i++;
                }
                const end = i;
                i--;
                extendGroupToWholeChains(loci, start, end, true, elements);
            } else {
                extendGroupToWholeChains(loci, i, i + 1, false, elements);
            }
        }

        return Loci(loci.structure, elements);
    }

    const boundaryHelper = new BoundaryHelper(), tempPos = Vec3.zero();
    export function getBoundary(loci: Loci): Boundary {
        boundaryHelper.reset(0);

        for (const e of loci.elements) {
            const { indices } = e;
            const pos = e.unit.conformation.position, r = e.unit.conformation.r;
            const { elements } = e.unit;
            for (let i = 0, _i = OrderedSet.size(indices); i < _i; i++) {
                const eI = elements[OrderedSet.getAt(indices, i)];
                pos(eI, tempPos);
                boundaryHelper.boundaryStep(tempPos, r(eI));
            }
        }
        boundaryHelper.finishBoundaryStep();
        for (const e of loci.elements) {
            const { indices } = e;
            const pos = e.unit.conformation.position, r = e.unit.conformation.r;
            const { elements } = e.unit;
            for (let i = 0, _i = OrderedSet.size(indices); i < _i; i++) {
                const eI = elements[OrderedSet.getAt(indices, i)];
                pos(eI, tempPos);
                boundaryHelper.extendStep(tempPos, r(eI));
            }
        }

        return { box: boundaryHelper.getBox(), sphere: boundaryHelper.getSphere() };
    }

    function sourceIndex(unit: Unit, element: ElementIndex) {
        return Unit.isAtomic(unit)
            ? unit.model.atomicHierarchy.atoms.sourceIndex.value(element)
            // TODO: when implemented, this should map to the source index.
            : element
    }

    export function toExpression(loci: Loci) {
        if (Loci.isEmpty(loci)) return MS.struct.generator.empty();

        const models = loci.structure.models;
        const sourceIndexMap = new Map<string, { modelLabel: string, modelIndex: number, xs: UniqueArray<number, number> }>();
        for (const e of loci.elements) {
            const { indices } = e;
            const { elements } = e.unit;

            const key = models.length === 1
                ? e.unit.conformation.operator.name
                : `${e.unit.conformation.operator.name} ${e.unit.model.label} ${e.unit.model.modelNum}`;

            let sourceIndices: UniqueArray<number, number>;
            if (sourceIndexMap.has(key)) sourceIndices = sourceIndexMap.get(key)!.xs;
            else {
                sourceIndices = UniqueArray.create<number, number>();
                sourceIndexMap.set(key, { modelLabel: e.unit.model.label, modelIndex: e.unit.model.modelNum, xs: sourceIndices });
            }

            for (let i = 0, _i = OrderedSet.size(indices); i < _i; i++) {
                const idx = sourceIndex(e.unit, elements[OrderedSet.getAt(indices, i)]);
                UniqueArray.add(sourceIndices, idx, idx);
            }
        }

        const opData: OpData[] = [];
        const keys = sourceIndexMap.keys();
        while (true) {
            const k = keys.next();
            if (k.done) break;
            const e = sourceIndexMap.get(k.value)!;
            opData.push(getOpData(k.value, e.xs.array, models.length > 1, e.modelLabel, e.modelIndex));
        }

        const opGroups = new Map<string, OpData>();
        for (let i = 0, il = opData.length; i < il; ++i) {
            const d = opData[i]
            const hash = hash2(hashFnv32a(d.atom.ranges), hashFnv32a(d.atom.set))
            const key = `${hash}|${d.entity ? (d.entity.modelLabel + d.entity.modelIndex) : ''}`
            if (opGroups.has(key)) {
                opGroups.get(key)!.chain.opName.push(...d.chain.opName)
            } else {
                opGroups.set(key, d)
            }
        }

        const opQueries: Expression[] = [];
        opGroups.forEach(d => {
            const { ranges, set } = d.atom
            const { opName } = d.chain

            const opProp = MS.struct.atomProperty.core.operatorName()
            const siProp = MS.struct.atomProperty.core.sourceIndex();
            const tests: Expression[] = [];

            // TODO: add set.ofRanges constructor to MolQL???
            if (set.length > 0) {
                tests[tests.length] = MS.core.set.has([MS.core.type.set(set), siProp]);
            }
            for (let rI = 0, _rI = ranges.length / 2; rI < _rI; rI++) {
                tests[tests.length] = MS.core.rel.inRange([siProp, ranges[2 * rI], ranges[2 * rI + 1]]);
            }

            if (d.entity) {
                const { modelLabel, modelIndex } = d.entity
                opQueries.push(MS.struct.generator.atomGroups({
                    'atom-test': tests.length > 1 ? MS.core.logic.or(tests) : tests[0],
                    'chain-test': opName.length > 1
                        ? MS.core.set.has([MS.core.type.set(opName), opProp])
                        : MS.core.rel.eq([opProp, opName[0]]),
                    'entity-test': MS.core.logic.and([
                        MS.core.rel.eq([MS.struct.atomProperty.core.modelLabel(), modelLabel]),
                        MS.core.rel.eq([MS.struct.atomProperty.core.modelIndex(), modelIndex]),
                    ])
                }))
            } else {
                opQueries.push(MS.struct.generator.atomGroups({
                    'atom-test': tests.length > 1 ? MS.core.logic.or(tests) : tests[0],
                    'chain-test': opName.length > 1
                        ? MS.core.set.has([MS.core.type.set(opName), opProp])
                        : MS.core.rel.eq([opProp, opName[0]])
                }))
            }
        })

        return MS.struct.modifier.union([
            opQueries.length === 1
                ? opQueries[0]
                // Need to union before merge for fast performance
                : MS.struct.combinator.merge(opQueries.map(q => MS.struct.modifier.union([ q ])))
        ]);
    }

    type OpData = {
        atom: { set: number[], ranges: number[] },
        chain: { opName: string[] },
        entity?: { modelLabel: string, modelIndex: number }
    }

    function getOpData(opName: string, xs: number[], multimodel: boolean, modelLabel: string, modelIndex: number): OpData {
        sortArray(xs);

        const ranges: number[] = [];
        const set: number[] = [];

        let i = 0, len = xs.length;
        while (i < len) {
            const start = i;
            i++;
            while (i < len && xs[i - 1] + 1 === xs[i]) i++;
            const end = i;
            // TODO: is this a good value?
            if (end - start > 12) {
                ranges[ranges.length] = xs[start];
                ranges[ranges.length] = xs[end - 1];
            } else {
                for (let j = start; j < end; j++) {
                    set[set.length] = xs[j];
                }
            }
        }

        return multimodel
            ? {
                atom: { set, ranges },
                chain: { opName: [ opName ] },
                entity: { modelLabel, modelIndex }
            }
            : {
                atom: { set, ranges },
                chain: { opName: [ opName ] },
            }
    }
}