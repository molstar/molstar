/**
 * Copyright (c) 2017-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ElementIndex } from '../../model';
import Unit from '../unit';
import { Vec3 } from '../../../../mol-math/linear-algebra';
import Structure from '../structure';

export { Location };

interface Location<U = Unit> {
    readonly kind: 'element-location',
    structure: Structure,
    unit: U,
    /** Index into element (atomic/coarse) properties of unit.model */
    element: ElementIndex
}

namespace Location {
    export function create<U extends Unit>(structure?: Structure, unit?: U, element?: ElementIndex): Location<U> {
        return {
            kind: 'element-location',
            structure: structure as any,
            unit: unit as any,
            element: element || (0 as ElementIndex)
        };
    }

    export function clone<U extends Unit>(l: Location<U>): Location<U> {
        return create(l.structure, l.unit, l.element);
    }

    export function set(a: Location, structure?: Structure, unit?: Unit, element?: ElementIndex): Location {
        if (structure) a.structure = structure;
        if (unit) a.unit = unit;
        if (element !== undefined) a.element = element;
        return a;
    }

    export function copy(out: Location, a: Location): Location {
        out.unit = a.unit;
        out.element = a.element;
        return out;
    }

    export function is(x: any): x is Location {
        return !!x && x.kind === 'element-location';
    }

    export function areEqual(a: Location, b: Location) {
        return a.unit === b.unit && a.element === b.element;
    }

    const pA = Vec3.zero(), pB = Vec3.zero();
    export function distance(a: Location, b: Location) {
        a.unit.conformation.position(a.element, pA);
        b.unit.conformation.position(b.element, pB);
        return Vec3.distance(pA, pB);
    }

    export function residueIndex(l: Location) {
        return l.unit.model.atomicHierarchy.residueAtomSegments.index[l.element];
    }

    export function chainIndex(l: Location) {
        return l.unit.model.atomicHierarchy.chainAtomSegments.index[l.element];
    }
}