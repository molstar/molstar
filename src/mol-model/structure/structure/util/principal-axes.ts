/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { PrincipalAxes } from '../../../../mol-math/linear-algebra/matrix/principal-axes';
import Unit from '../unit';
import { Vec3 } from '../../../../mol-math/linear-algebra';

const tempPos = Vec3.zero();
export function toPositionsArray(unit: Unit) {
    const pos = unit.conformation.invariantPosition;
    const { elements } = unit;
    const positions = new Float32Array(elements.length * 3);
    for (let i = 0, il = elements.length; i < il; i++) {
        pos(elements[i], tempPos);
        Vec3.toArray(tempPos, positions, i * 3);
    }
    return positions;
}

export function getPrincipalAxes(unit: Unit): PrincipalAxes {
    const positions = toPositionsArray(unit);
    return PrincipalAxes.ofPositions(positions);
}