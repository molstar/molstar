/**
 * Copyright (c) 2017-2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Unit, StructureElement } from '../../structure'
import Structure from '../structure';
import { LinkType } from '../../model/types';

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
        readonly links: ReadonlyArray<Location>
    }

    export function Loci(links: ArrayLike<Location>): Loci {
        return { kind: 'link-loci', links: links as Loci['links'] };
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