/**
 * Copyright (c) 2017-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Unit, StructureElement } from '../../structure'
import Structure from '../structure';
import { LinkType } from '../../model/types';
import { SortedArray } from '../../../../mol-data/int';

export * from './links/data'
export * from './links/intra-compute'
export * from './links/inter-compute'

namespace Link {
    export interface Location<U extends Unit = Unit> {
        readonly kind: 'link-location',
        aUnit: U,
        /** Index into aUnit.elements */
        aIndex: StructureElement.UnitIndex,
        bUnit: U,
        /** Index into bUnit.elements */
        bIndex: StructureElement.UnitIndex,
    }

    export function Location(aUnit?: Unit, aIndex?: StructureElement.UnitIndex, bUnit?: Unit, bIndex?: StructureElement.UnitIndex): Location {
        return { kind: 'link-location', aUnit: aUnit as any, aIndex: aIndex as any, bUnit: bUnit as any, bIndex: bIndex as any };
    }

    export function isLocation(x: any): x is Location {
        return !!x && x.kind === 'link-location';
    }

    export function areLocationsEqual(locA: Location, locB: Location) {
        return (
            locA.aIndex === locB.aIndex && locA.bIndex === locB.bIndex &&
            locA.aUnit.id === locB.aUnit.id && locA.bUnit.id === locB.bUnit.id
        )
    }

    export interface Loci {
        readonly kind: 'link-loci',
        readonly structure: Structure
        readonly links: ReadonlyArray<Location>
    }

    export function Loci(structure: Structure, links: ArrayLike<Location>): Loci {
        return { kind: 'link-loci', structure, links: links as Loci['links'] };
    }

    export function isLoci(x: any): x is Loci {
        return !!x && x.kind === 'link-loci';
    }

    export function areLociEqual(a: Loci, b: Loci) {
        if (a.links.length !== b.links.length) return false
        for (let i = 0, il = a.links.length; i < il; ++i) {
            if (!areLocationsEqual(a.links[i], b.links[i])) return false
        }
        return true
    }

    export function isLociEmpty(loci: Loci) {
        return loci.links.length === 0 ? true : false
    }

    export function remapLoci(loci: Loci, structure: Structure): Loci {
        if (structure === loci.structure) return loci

        const links: Loci['links'][0][] = [];
        loci.links.forEach(l => {
            const unitA = structure.unitMap.get(l.aUnit.id)
            if (!unitA) return
            const unitB = structure.unitMap.get(l.bUnit.id)
            if (!unitB) return

            const elementA = l.aUnit.elements[l.aIndex]
            const indexA = SortedArray.indexOf(unitA.elements, elementA) as StructureElement.UnitIndex | -1
            if (indexA === -1) return
            const elementB = l.bUnit.elements[l.bIndex]
            const indexB = SortedArray.indexOf(unitB.elements, elementB) as StructureElement.UnitIndex | -1
            if (indexB === -1) return

            links.push(Location(unitA, indexA, unitB, indexB))
        });

        return Loci(structure, links);
    }

    export function toStructureElementLoci(loci: Loci): StructureElement.Loci {
        const elements: StructureElement.Loci['elements'][0][] = []
        const map = new Map<number, number[]>()

        for (const lociLink of loci.links) {
            const { aIndex, aUnit, bIndex, bUnit } = lociLink
            if (aUnit === bUnit) {
                if (map.has(aUnit.id)) map.get(aUnit.id)!.push(aIndex, bIndex)
                else map.set(aUnit.id, [aIndex, bIndex])
            } else {
                if (map.has(aUnit.id)) map.get(aUnit.id)!.push(aIndex)
                else map.set(aUnit.id, [aIndex])
                if (map.has(bUnit.id)) map.get(bUnit.id)!.push(bIndex)
                else map.set(bUnit.id, [bIndex])
            }
        }

        map.forEach((indices: number[], id: number) => {
            elements.push({
                unit: loci.structure.unitMap.get(id)!,
                indices: SortedArray.deduplicate(SortedArray.ofUnsortedArray(indices))
            })
        })

        return StructureElement.Loci(loci.structure, elements);
    }

    export function getType(structure: Structure, link: Location<Unit.Atomic>): LinkType {
        if (link.aUnit === link.bUnit) {
            const links = link.aUnit.links;
            const idx = links.getEdgeIndex(link.aIndex, link.bIndex);
            if (idx < 0) return LinkType.create(LinkType.Flag.None);
            return LinkType.create(links.edgeProps.flags[idx]);
        } else {
            const bond = structure.links.getBondFromLocation(link);
            if (bond) return LinkType.create(bond.flag);
            return LinkType.create(LinkType.Flag.None);
        }
    }

    export function getOrder(structure: Structure, link: Location<Unit.Atomic>): number {
        if (link.aUnit === link.bUnit) {
            const links = link.aUnit.links;
            const idx = links.getEdgeIndex(link.aIndex, link.bIndex);
            if (idx < 0) return 0;
            return links.edgeProps.order[idx];
        } else {
            const bond = structure.links.getBondFromLocation(link);
            if (bond) return bond.order;
            return 0;
        }
    }
}

export { Link }