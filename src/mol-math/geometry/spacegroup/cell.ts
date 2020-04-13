/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Vec3 } from '../../linear-algebra';

export { Cell };

interface Cell {
    readonly size: Vec3
    readonly anglesInRadians: Vec3
}

function Cell() {
    return Cell.empty();
}

namespace Cell {
    export function create(size: Vec3, anglesInRadians: Vec3): Cell { return { size, anglesInRadians }; }
    export function empty(): Cell { return { size: Vec3(), anglesInRadians: Vec3() }; }
}