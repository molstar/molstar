/**
 * Copyright (c) 2017-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ElementIndex } from '../../model';
import Unit from '../unit';

export { Location }

interface Location<U = Unit> {
    readonly kind: 'element-location',
    unit: U,
    /** Index into element (atomic/coarse) properties of unit.model */
    element: ElementIndex
}

namespace Location {
    export function create(unit?: Unit, element?: ElementIndex): Location {
        return { kind: 'element-location', unit: unit!, element: element || (0 as ElementIndex) };
    }

    export function set(a: Location, unit?: Unit, element?: ElementIndex): Location {
        if (unit) a.unit = unit
        if (element !== undefined) a.element = element
        return a;
    }

    export function copy(out: Location, a: Location): Location {
        out.unit = a.unit
        out.element = a.element
        return out
    }

    export function is(x: any): x is Location {
        return !!x && x.kind === 'element-location';
    }
}