/**
 * Copyright (c) 2019-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Panagiotis Tourlas <panagiot_tourlov@hotmail.com>
 */

import { parseMol, formalChargeMapper } from '../mol/parser';

const MolString = `2244
  -OEChem-04072009073D

 21 21  0     0  0  0  0  0  0999 V2000
    1.2333    0.5540    0.7792 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6952   -2.7148   -0.7502 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.7958   -2.1843    0.8685 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.7813    0.8105   -1.4821 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0857    0.6088    0.4403 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7927   -0.5515    0.1244 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7288    1.8464    0.4133 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1426   -0.4741   -0.2184 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.0787    1.9238    0.0706 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.7855    0.7636   -0.2453 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.1409   -1.8536    0.1477 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.1094    0.6715   -0.3113 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.5305    0.5996    0.1635 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.1851    2.7545    0.6593 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.7247   -1.3605   -0.4564 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.5797    2.8872    0.0506 H   0  0  0  0  0  0  0  0  0  0  0  0
   -3.8374    0.8238   -0.5090 H   0  0  0  0  0  0  0  0  0  0  0  0
    3.7290    1.4184    0.8593 H   0  0  0  0  0  0  0  0  0  0  0  0
    4.2045    0.6969   -0.6924 H   0  0  0  0  0  0  0  0  0  0  0  0
    3.7105   -0.3659    0.6426 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2555   -3.5916   -0.7337 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  5  1  0  0  0  0
  1 12  1  0  0  0  0
  2 11  1  0  0  0  0
  2 21  1  0  0  0  0
  3 11  2  0  0  0  0
  4 12  2  0  0  0  0
  5  6  1  0  0  0  0
  5  7  2  0  0  0  0
  6  8  2  0  0  0  0
  6 11  1  0  0  0  0
  7  9  1  0  0  0  0
  7 14  1  0  0  0  0
  8 10  1  0  0  0  0
  8 15  1  0  0  0  0
  9 10  2  0  0  0  0
  9 16  1  0  0  0  0
 10 17  1  0  0  0  0
 12 13  1  0  0  0  0
 13 18  1  0  0  0  0
 13 19  1  0  0  0  0
 13 20  1  0  0  0  0
M  END`;

const MolStringWithAtomBlockCharge = `
  Ketcher  1 72215442D 1   1.00000     0.00000     0

  4  3  0  0  0  0            999 V2000
    0.0000    0.0000    0.0000 C   0  1  0  0  0  0  0  0  0  0  0  0
    0.8660    0.5000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8660    0.5000    0.0000 S   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000   -1.0000    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
  1  4  2  0  0  0  0
  3  1  1  0  0  0  0
  2  1  1  0  0  0  0
M  END`;

const MolStringWithPropertyBlockCharge = `
  Ketcher  1 72215442D 1   1.00000     0.00000     0

  4  3  0  0  0  0            999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.8660    0.5000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8660    0.5000    0.0000 S   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000   -1.0000    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
  1  4  2  0  0  0  0
  3  1  1  0  0  0  0
  2  1  1  0  0  0  0
M  CHG  3   2  -1   3   1   4   1
M  END`;

const MolStringWithMultipleChargeLines = `
  Ketcher  1 72215442D 1   1.00000     0.00000     0

  4  3  0  0  0  0            999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.8660    0.5000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8660    0.5000    0.0000 S   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000   -1.0000    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
  1  4  2  0  0  0  0
  3  1  1  0  0  0  0
  2  1  1  0  0  0  0
M  CHG  1   2  -1
M  CHG  2   3   1   4   1
M  END`;

describe('mol reader', () => {
    it('basic', async () => {
        const parsed = await parseMol(MolString).run();
        if (parsed.isError) {
            throw new Error(parsed.message);
        }
        const { atoms, bonds } = parsed.result;

        // number of structures
        expect(atoms.count).toBe(21);
        expect(bonds.count).toBe(21);

        expect(atoms.x.value(0)).toBeCloseTo(1.2333, 0.001);
        expect(atoms.y.value(0)).toBeCloseTo(0.5540, 0.0001);
        expect(atoms.z.value(0)).toBeCloseTo(0.7792, 0.0001);
        expect(atoms.type_symbol.value(0)).toBe('O');

        expect(bonds.atomIdxA.value(20)).toBe(13);
        expect(bonds.atomIdxB.value(20)).toBe(20);
        expect(bonds.order.value(20)).toBe(1);
    });
    it('property block charges', async () => {
        const parsed = await parseMol(MolStringWithPropertyBlockCharge).run();
        if (parsed.isError) {
            throw new Error(parsed.message);
        }
        const { formalCharges } = parsed.result;

        expect(formalCharges.atomIdx.rowCount).toBe(3);
        expect(formalCharges.charge.rowCount).toBe(3);

        expect(formalCharges.atomIdx.value(0)).toBe(2);
        expect(formalCharges.atomIdx.value(1)).toBe(3);

        expect(formalCharges.charge.value(0)).toBe(-1);
        expect(formalCharges.charge.value(1)).toBe(1);
    });
    it('multiple charge lines', async () => {
        const parsed = await parseMol(MolStringWithMultipleChargeLines).run();
        if (parsed.isError) {
            throw new Error(parsed.message);
        }
        const { formalCharges } = parsed.result;

        expect(formalCharges.atomIdx.rowCount).toBe(3);
        expect(formalCharges.charge.rowCount).toBe(3);

        expect(formalCharges.atomIdx.value(0)).toBe(2);
        expect(formalCharges.atomIdx.value(1)).toBe(3);

        expect(formalCharges.charge.value(0)).toBe(-1);
        expect(formalCharges.charge.value(1)).toBe(1);
    });

    it('atom block charge mapping', async () => {
        expect(formalChargeMapper(7)).toBe(-3);
        expect(formalChargeMapper(6)).toBe(-2);
        expect(formalChargeMapper(5)).toBe(-1);
        expect(formalChargeMapper(0)).toBe(0);
        expect(formalChargeMapper(3)).toBe(1);
        expect(formalChargeMapper(2)).toBe(2);
        expect(formalChargeMapper(1)).toBe(3);
        expect(formalChargeMapper(4)).toBe(0);
    });
    it('atom block charges', async () => {
        const parsed = await parseMol(MolStringWithAtomBlockCharge).run();
        if (parsed.isError) {
            throw new Error(parsed.message);
        }
        const { atoms, formalCharges } = parsed.result;

        /* No property block charges */
        expect(formalCharges.atomIdx.rowCount).toBe(0);
        expect(formalCharges.charge.rowCount).toBe(0);

        expect(atoms.formal_charge.value(0)).toBe(1);
        expect(atoms.formal_charge.value(1)).toBe(0);
        expect(atoms.formal_charge.value(2)).toBe(0);
        expect(atoms.formal_charge.value(3)).toBe(0);
    });
});
