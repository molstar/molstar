/**
 * Copyright (c) 2017-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { UniqueArray } from '../../../mol-data/generic';
import { OrderedSet, SortedArray } from '../../../mol-data/int';
import { BoundaryHelper } from '../../../mol-math/geometry/boundary-helper';
import { Vec3 } from '../../../mol-math/linear-algebra';
import { MolScriptBuilder as MS } from '../../../mol-script/language/builder';
import { ElementIndex } from '../model';
import { ChainIndex, ResidueIndex } from '../model/indexing';
import Structure from './structure';
import Unit from './unit';
import { Boundary } from './util/boundary';
import { StructureProperties } from '../structure';
import { sortArray, hashFnv32a, hash2 } from '../../../mol-data/util';
import Expression from '../../../mol-script/language/expression';

interface StructureElement<U = Unit> {
    readonly kind: 'element-location',
    unit: U,
    /** Index into element (atomic/coarse) properties of unit.model */
    element: ElementIndex
}

namespace StructureElement {
    export function create(unit?: Unit, element?: ElementIndex): StructureElement {
        return { kind: 'element-location', unit: unit!, element: element || (0 as ElementIndex) };
    }

    export function set(a: StructureElement, unit?: Unit, element?: ElementIndex): StructureElement {
        if (unit) a.unit = unit
        if (element !== undefined) a.element = element
        return a;
    }

    // TODO: when nominal types are available, make this indexed by UnitIndex
    export type Set = SortedArray<ElementIndex>

    /** Index into Unit.elements */
    export type UnitIndex = { readonly '@type': 'structure-element-index' } & number

    export interface Property<T> { (location: StructureElement): T }
    export interface Predicate extends Property<boolean> { }

    export function property<T>(p: Property<T>) { return p; }

    function _wrongUnitKind(kind: string) { throw new Error(`Property only available for ${kind} models.`); }
    export function atomicProperty<T>(p: (location: StructureElement<Unit.Atomic>) => T) {
        return property(l => Unit.isAtomic(l.unit) ? p(l as StructureElement<Unit.Atomic>) : _wrongUnitKind('atomic'));
    }

    export function coarseProperty<T>(p: (location: StructureElement<Unit.Spheres | Unit.Gaussians>) => T) {
        return property(l => Unit.isCoarse(l.unit) ? p(l as StructureElement<Unit.Spheres | Unit.Gaussians>) : _wrongUnitKind('coarse'));
    }

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

    export function isLoci(x: any): x is Loci {
        return !!x && x.kind === 'element-loci';
    }

    export function areLociEqual(a: Loci, b: Loci) {
        if (a.elements.length !== b.elements.length) return false
        for (let i = 0, il = a.elements.length; i < il; ++i) {
            const elementA = a.elements[i]
            const elementB = b.elements[i]
            if (elementA.unit.id !== elementB.unit.id) return false
            if (!OrderedSet.areEqual(elementA.indices, elementB.indices)) return false
        }
        return true
    }

    export function isLocation(x: any): x is StructureElement {
        return !!x && x.kind === 'element-location';
    }

    export function residueIndex(e: StructureElement) {
        if (Unit.isAtomic(e.unit)) {
            return e.unit.residueIndex[e.element];
        } else {
            // TODO: throw error instead?
            return -1 as ResidueIndex;
        }
    }

    export function chainIndex(e: StructureElement) {
        if (Unit.isAtomic(e.unit)) {
            return e.unit.chainIndex[e.element];
        } else {
            // TODO: throw error instead?
            return -1 as ChainIndex;
        }
    }

    export function entityIndex(l: StructureElement) {
        return StructureProperties.entity.key(l)
        }

    export namespace Loci {
        export function size(loci: Loci) {
            let s = 0;
            for (const u of loci.elements) s += OrderedSet.size(u.indices);
            return s;
        }

        export function all(structure: Structure): Loci {
            return Loci(structure, structure.units.map(unit => ({
                unit,
                indices: OrderedSet.ofRange<UnitIndex>(0 as UnitIndex, unit.elements.length as UnitIndex)
            })));
        }

        export function remap(loci: Loci, structure: Structure): Loci {
            if (structure === loci.structure) return loci

            const elements: Loci['elements'][0][] = [];
            loci.elements.forEach(e => {
                const unit = structure.unitMap.get(e.unit.id)
                if (!unit) return

                if (SortedArray.areEqual(e.unit.elements, unit.elements)) {
                    elements.push({ unit, indices: e.indices })
                } else {
                    // TODO optimize
                    const indices: UnitIndex[] = []
                    OrderedSet.forEach(e.indices, (v) => {
                        const eI = e.unit.elements[v]
                        const uI = SortedArray.indexOf(unit.elements, eI) as UnitIndex | -1
                        if (uI !== -1) indices.push(uI)
                    })
                    elements.push({
                        unit,
                        indices: SortedArray.ofSortedArray(indices)
                    })
                }
            });

            return Loci(structure, elements);
        }

        /** Create union of `xs` and `ys` */
        export function union(xs: Loci, ys: Loci): Loci {
            if (xs.elements.length > ys.elements.length) return union(ys, xs);
            if (xs.elements.length === 0) return ys;

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
            if (xs.elements.length === 0) return ys.elements.length === 0;

            const map = new Map<number, OrderedSet<UnitIndex>>();

            for (const e of xs.elements) map.set(e.unit.id, e.indices);
            for (const e of ys.elements) {
                if (!map.has(e.unit.id)) continue;
                if (OrderedSet.areIntersecting(map.get(e.unit.id)!, e.indices)) return true;
            }

            return false;
        }

        export function extendToWholeResidues(loci: Loci): Loci {
            const elements: Loci['elements'][0][] = [];

            for (const lociElement of loci.elements) {
                if (lociElement.unit.kind === Unit.Kind.Atomic) {
                    const unitElements = lociElement.unit.elements;
                    const h = lociElement.unit.model.atomicHierarchy;

                    const { index: residueIndex, offsets: residueOffsets } = h.residueAtomSegments;

                    const newIndices: UnitIndex[] = [];
                    const indices = lociElement.indices, len = OrderedSet.size(indices);
                    let i = 0;
                    while (i < len) {
                        const rI = residueIndex[unitElements[OrderedSet.getAt(indices, i)]];
                        i++;
                        while (i < len && residueIndex[unitElements[OrderedSet.getAt(indices, i)]] === rI) {
                            i++;
                        }

                        for (let j = residueOffsets[rI], _j = residueOffsets[rI + 1]; j < _j; j++) {
                            const idx = OrderedSet.indexOf(unitElements, j);
                            if (idx >= 0) newIndices[newIndices.length] = idx as UnitIndex;
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

        export function extendToWholeChains(loci: Loci): Loci {
            const elements: Loci['elements'][0][] = [];

            for (const lociElement of loci.elements) {
                const newIndices: UnitIndex[] = [];
                const unitElements = lociElement.unit.elements;

                const { index: chainIndex, offsets: chainOffsets } = getChainSegments(lociElement.unit)

                const indices = lociElement.indices, len = OrderedSet.size(indices);
                let i = 0;
                while (i < len) {
                    const cI = chainIndex[unitElements[OrderedSet.getAt(indices, i)]];
                    i++;
                    while (i < len && chainIndex[unitElements[OrderedSet.getAt(indices, i)]] === cI) {
                        i++;
                    }

                    for (let j = chainOffsets[cI], _j = chainOffsets[cI + 1]; j < _j; j++) {
                        const idx = OrderedSet.indexOf(unitElements, j);
                        if (idx >= 0) newIndices[newIndices.length] = idx as UnitIndex;
                    }
                }

                elements[elements.length] = { unit: lociElement.unit, indices: SortedArray.ofSortedArray(newIndices) };
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

        export function toScriptExpression(loci: Loci) {
            if (loci.elements.length === 0) return MS.struct.generator.empty();

            const models = loci.structure.models;
            const sourceIndexMap = new Map<string, { modelLabel: string, modelIndex: number, xs: UniqueArray<number, number> }>();
            const el = StructureElement.create(), p = StructureProperties.atom.sourceIndex;
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

                el.unit = e.unit;
                for (let i = 0, _i = OrderedSet.size(indices); i < _i; i++) {
                    el.element = elements[OrderedSet.getAt(indices, i)];
                    const idx = p(el);
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
}

export default StructureElement