/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Spacegroup, SpacegroupCell } from '../spacegroup/construction';
import { Vec3 } from '../../linear-algebra';

function getSpacegroup(name: string) {
    const size = Vec3.create(1, 1, 1);
    const anglesInRadians = Vec3.create(Math.PI / 2, Math.PI / 2, Math.PI / 2);
    const cell = SpacegroupCell.create(name, size, anglesInRadians);
    return Spacegroup.create(cell);
}

function checkOperatorsXyz(name: string, expected: string[]) {
    const spacegroup = getSpacegroup(name);
    for (let i = 0, il = spacegroup.operators.length; i < il; ++i) {
        const op = spacegroup.operators[i];
        const actual = Spacegroup.getOperatorXyz(op);
        expect(actual).toBe(expected[i]);
    }
}

describe('Spacegroup', () => {
    it('operators xyz', () => {
        checkOperatorsXyz('P 1', ['X,Y,Z']);
        checkOperatorsXyz('P -1', ['X,Y,Z', '-X,-Y,-Z']);
        checkOperatorsXyz('P 1 21 1', ['X,Y,Z', '-X,1/2+Y,-Z']);
        checkOperatorsXyz('P 1 21/m 1', ['X,Y,Z', '-X,1/2+Y,-Z', '-X,-Y,-Z', 'X,1/2-Y,Z']);
        checkOperatorsXyz('P 41', ['X,Y,Z', '-X,-Y,1/2+Z', '-Y,X,1/4+Z', 'Y,-X,3/4+Z']);
        checkOperatorsXyz('P 41 21 2', ['X,Y,Z', '-X,-Y,1/2+Z', '1/2-Y,1/2+X,1/4+Z', '1/2+Y,1/2-X,3/4+Z', '1/2-X,1/2+Y,1/4-Z', '1/2+X,1/2-Y,3/4-Z', 'Y,X,-Z', '-Y,-X,1/2-Z']);
        checkOperatorsXyz('P 3', ['X,Y,Z', '-Y,X-Y,Z', 'Y-X,-X,Z']);
    });
});
