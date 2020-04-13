/**
 * Copyright (c) 2017-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Unit, StructureElement } from '../../structure';
import Structure from '../structure';
import { BondType } from '../../model/types';
import { SortedArray, Iterator } from '../../../../mol-data/int';
import { CentroidHelper } from '../../../../mol-math/geometry/centroid-helper';
import { Sphere3D } from '../../../../mol-math/geometry';

export * from './bonds/data';
export * from './bonds/intra-compute';
export * from './bonds/inter-compute';

namespace Bond {
    export interface Location<U extends Unit = Unit> {
        readonly kind: 'bond-location',

        aStructure: Structure,
        aUnit: U,
        /** Index into aUnit.elements */
        aIndex: StructureElement.UnitIndex,

        bStructure: Structure,
        bUnit: U,
        /** Index into bUnit.elements */
        bIndex: StructureElement.UnitIndex,
    }

    export function Location(aStructure?: Structure, aUnit?: Unit, aIndex?: StructureElement.UnitIndex, bStructure?: Structure, bUnit?: Unit, bIndex?: StructureElement.UnitIndex): Location {
        return {
            kind: 'bond-location',
            aStructure: aStructure as any,
            aUnit: aUnit as any,
            aIndex: aIndex as any,
            bStructure: bStructure as any,
            bUnit: bUnit as any,
            bIndex: bIndex as any
        };
    }

    export function isLocation(x: any): x is Location {
        return !!x && x.kind === 'bond-location';
    }

    export function areLocationsEqual(locA: Location, locB: Location) {
        return (
            locA.aStructure.label === locB.aStructure.label && locA.bStructure.label === locB.bStructure.label &&
            locA.aIndex === locB.aIndex && locA.bIndex === locB.bIndex &&
            locA.aUnit.id === locB.aUnit.id && locA.bUnit.id === locB.bUnit.id
        );
    }

    export interface Loci {
        readonly kind: 'bond-loci',
        readonly structure: Structure
        readonly bonds: ReadonlyArray<Location>
    }

    export function Loci(structure: Structure, bonds: ArrayLike<Location>): Loci {
        return { kind: 'bond-loci', structure, bonds: bonds as Loci['bonds'] };
    }

    export function isLoci(x: any): x is Loci {
        return !!x && x.kind === 'bond-loci';
    }

    export function areLociEqual(a: Loci, b: Loci) {
        if (a.structure !== b.structure) return false;
        if (a.bonds.length !== b.bonds.length) return false;
        for (let i = 0, il = a.bonds.length; i < il; ++i) {
            if (!areLocationsEqual(a.bonds[i], b.bonds[i])) return false;
        }
        return true;
    }

    export function isLociEmpty(loci: Loci) {
        return loci.bonds.length === 0 ? true : false;
    }

    export function remapLoci(loci: Loci, structure: Structure): Loci {
        if (structure === loci.structure) return loci;

        const bonds: Loci['bonds'][0][] = [];
        loci.bonds.forEach(l => {
            const unitA = structure.unitMap.get(l.aUnit.id);
            if (!unitA) return;
            const unitB = structure.unitMap.get(l.bUnit.id);
            if (!unitB) return;

            const elementA = l.aUnit.elements[l.aIndex];
            const indexA = SortedArray.indexOf(unitA.elements, elementA) as StructureElement.UnitIndex | -1;
            if (indexA === -1) return;
            const elementB = l.bUnit.elements[l.bIndex];
            const indexB = SortedArray.indexOf(unitB.elements, elementB) as StructureElement.UnitIndex | -1;
            if (indexB === -1) return;

            bonds.push(Location(loci.structure, unitA, indexA, loci.structure, unitB, indexB));
        });

        return Loci(structure, bonds);
    }

    export function toStructureElementLoci(loci: Loci): StructureElement.Loci {
        const elements: StructureElement.Loci['elements'][0][] = [];
        const map = new Map<number, number[]>();

        for (const lociBond of loci.bonds) {
            const { aIndex, aUnit, bIndex, bUnit } = lociBond;
            if (aUnit === bUnit) {
                if (map.has(aUnit.id)) map.get(aUnit.id)!.push(aIndex, bIndex);
                else map.set(aUnit.id, [aIndex, bIndex]);
            } else {
                if (map.has(aUnit.id)) map.get(aUnit.id)!.push(aIndex);
                else map.set(aUnit.id, [aIndex]);
                if (map.has(bUnit.id)) map.get(bUnit.id)!.push(bIndex);
                else map.set(bUnit.id, [bIndex]);
            }
        }

        map.forEach((indices: number[], id: number) => {
            elements.push({
                unit: loci.structure.unitMap.get(id)!,
                indices: SortedArray.deduplicate(SortedArray.ofUnsortedArray(indices))
            });
        });

        return StructureElement.Loci(loci.structure, elements);
    }

    export function getType(structure: Structure, location: Location<Unit.Atomic>): BondType {
        if (location.aUnit === location.bUnit) {
            const bonds = location.aUnit.bonds;
            const idx = bonds.getEdgeIndex(location.aIndex, location.bIndex);
            if (idx < 0) return BondType.create(BondType.Flag.None);
            return BondType.create(bonds.edgeProps.flags[idx]);
        } else {
            const bond = structure.interUnitBonds.getBondFromLocation(location);
            if (bond) return BondType.create(bond.props.flag);
            return BondType.create(BondType.Flag.None);
        }
    }

    export function getOrder(structure: Structure, location: Location<Unit.Atomic>): number {
        if (location.aUnit === location.bUnit) {
            const bonds = location.aUnit.bonds;
            const idx = bonds.getEdgeIndex(location.aIndex, location.bIndex);
            if (idx < 0) return 0;
            return bonds.edgeProps.order[idx];
        } else {
            const bond = structure.interUnitBonds.getBondFromLocation(location);
            if (bond) return bond.props.order;
            return 0;
        }
    }

    export function getIntraUnitBondCount(structure: Structure) {
        let count = 0;
        for (let i = 0, il = structure.units.length; i < il; ++i) {
            const u = structure.units[i];
            if (Unit.isAtomic(u)) count += u.bonds.edgeCount / 2; // only count one direction
        }
        return count;
    }

    export interface ElementBondData {
        otherUnit: Unit.Atomic
        otherIndex: StructureElement.UnitIndex
        type: BondType
        order: number
    }

    export class ElementBondIterator implements Iterator<ElementBondData> {
        private current: ElementBondData = {} as any

        private structure: Structure
        private unit: Unit.Atomic
        private index: StructureElement.UnitIndex

        private interBondIndices: ReadonlyArray<number>
        private interBondCount: number
        private interBondIndex: number

        private intraBondEnd: number
        private intraBondIndex: number

        hasNext: boolean;
        move(): ElementBondData {
            this.advance();
            return this.current;
        }

        setElement(structure: Structure, unit: Unit.Atomic, index: StructureElement.UnitIndex) {
            this.structure = structure;
            this.unit = unit;
            this.index = index;

            this.interBondIndices = structure.interUnitBonds.getEdgeIndices(index, unit);
            this.interBondCount = this.interBondIndices.length;
            this.interBondIndex = 0;

            this.intraBondEnd = unit.bonds.offset[index + 1];
            this.intraBondIndex = unit.bonds.offset[index];

            this.hasNext = this.interBondIndex < this.interBondCount || this.intraBondIndex < this.intraBondEnd;
        }

        private advance() {
            if (this.intraBondIndex < this.intraBondEnd) {
                this.current.otherUnit = this.unit;
                this.current.otherIndex = this.unit.bonds.b[this.intraBondIndex] as StructureElement.UnitIndex;
                this.current.type = this.unit.bonds.edgeProps.flags[this.intraBondIndex];
                this.current.order = this.unit.bonds.edgeProps.order[this.intraBondIndex];
                this.intraBondIndex += 1;
            } else if (this.interBondIndex < this.interBondCount) {
                const b = this.structure.interUnitBonds.edges[this.interBondIndex];
                this.current.otherUnit = b.unitA !== this.unit ? b.unitA : b.unitB;
                this.current.otherIndex = b.indexA !== this.index ? b.indexA : b.indexB;
                this.current.type = b.props.flag;
                this.current.order = b.props.order;
                this.interBondIndex += 1;
            } else {
                this.hasNext = false;
                return;
            }
            this.hasNext = this.interBondIndex < this.interBondCount || this.intraBondIndex < this.intraBondEnd;
        }

        constructor() {
            this.hasNext = false;
        }
    }

    export function getBoundingSphere(loci: Loci, boundingSphere: Sphere3D) {
        return CentroidHelper.fromPairProvider(loci.bonds.length, (i, pA, pB) => {
            const { aUnit, aIndex, bUnit, bIndex } = loci.bonds[i];
            aUnit.conformation.position(aUnit.elements[aIndex], pA);
            bUnit.conformation.position(bUnit.elements[bIndex], pB);
        }, boundingSphere);
    }
}

export { Bond };