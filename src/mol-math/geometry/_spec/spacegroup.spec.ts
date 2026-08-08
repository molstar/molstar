/**
 * Copyright (c) 2020-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Spacegroup, SpacegroupCell } from '../spacegroup/construction';
import { Cell } from '../spacegroup/cell';
import { operatorsForEntry, orderForEntry, SpacegroupData } from '../spacegroup/tables';
import { Vec3 } from '../../linear-algebra';

function getSpacegroup(name: string) {
    const size = Vec3.create(1, 1, 1);
    const anglesInRadians = Vec3.create(Math.PI / 2, Math.PI / 2, Math.PI / 2);
    const cell = SpacegroupCell.create(name, size, anglesInRadians);
    return Spacegroup.create(cell);
}

function checkOperatorsXyz(name: string, expected: string[]) {
    const spacegroup = getSpacegroup(name);
    const actual: string[] = [];
    for (let i = 0, il = spacegroup.operators.length; i < il; ++i) {
        actual.push(Spacegroup.getOperatorXyz(spacegroup.operators[i]));
    }
    // Compare as sets: Spacegroup.create now derives operators from Hall symbols,
    // which yields the same operators but in a different order than the legacy tables.
    expect(actual.slice().sort()).toEqual(expected.slice().sort());
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

    it('orderForEntry matches the generated operator count', () => {
        for (const entry of SpacegroupData) {
            expect([entry.names[0], orderForEntry(entry)]).toEqual([entry.names[0], operatorsForEntry(entry).length]);
        }
    });
});

describe('Cell', () => {
    const angles90 = Vec3.create(Math.PI / 2, Math.PI / 2, Math.PI / 2);

    it('volume', () => {
        expect(Cell.create(Vec3.create(2, 3, 4), angles90).volume).toBeCloseTo(24, 6);
        // triclinic, all angles 60 degrees
        const angles60 = Vec3.create(Math.PI / 3, Math.PI / 3, Math.PI / 3);
        expect(Cell.create(Vec3.create(10, 10, 10), angles60).volume).toBeCloseTo(707.10678, 4);
        expect(Cell.empty().volume).toBe(0);
    });

    it('order', () => {
        expect(Cell.create(Vec3.create(1, 1, 1), angles90).order).toBe(1);
        expect(SpacegroupCell.create('P 1', Vec3.create(10, 10, 10), angles90).order).toBe(1);
        expect(SpacegroupCell.create('P 21 21 21', Vec3.create(10, 10, 10), angles90).order).toBe(4);
        expect(SpacegroupCell.create('P 41 21 2', Vec3.create(10, 10, 10), angles90).order).toBe(8);
    });
});
