/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { CIF } from '../../../mol-io/reader/cif';

const cifCoreString = `data_n1379
_audit_block_doi                 10.5517/ccy42jn
_database_code_depnum_ccdc_archive 'CCDC 867861'
loop_
_citation_id
_citation_doi
_citation_year
1 10.1002/chem.201202070 2012
_audit_update_record
;
2012-02-20 deposited with the CCDC.
2016-10-08 downloaded from the CCDC.
;

loop_
_atom_type_symbol
_atom_type_description
_atom_type_scat_dispersion_real
_atom_type_scat_dispersion_imag
_atom_type_scat_source
C C 0.0181 0.0091 'International Tables Vol C Tables 4.2.6.8 and 6.1.1.4'
H H 0.0000 0.0000 'International Tables Vol C Tables 4.2.6.8 and 6.1.1.4'
N N 0.0311 0.0180 'International Tables Vol C Tables 4.2.6.8 and 6.1.1.4'
O O 0.0492 0.0322 'International Tables Vol C Tables 4.2.6.8 and 6.1.1.4'
F F 0.0727 0.0534 'International Tables Vol C Tables 4.2.6.8 and 6.1.1.4'

_cell_length_a                   11.0829(8)
_cell_length_b                   14.6829(10)
_cell_length_c                   16.8532(17)
_cell_angle_alpha                105.728(6)
_cell_angle_beta                 100.310(6)
_cell_angle_gamma                110.620(4)
_cell_volume                     2353.3(3)
_cell_formula_units_Z            1
_cell_measurement_temperature    100(2)
_cell_measurement_reflns_used    5934
_cell_measurement_theta_min      2.86
_cell_measurement_theta_max      64.30

loop_
_atom_site_aniso_label
_atom_site_aniso_U_11
_atom_site_aniso_U_22
_atom_site_aniso_U_33
_atom_site_aniso_U_23
_atom_site_aniso_U_13
_atom_site_aniso_U_12
Pt1 0.0425(2) 0.0423(2) 0.0375(2) 0.00066(13) 0.01515(13) 0.00089(12)
K1 0.0605(15) 0.0687(17) 0.0559(17) 0.000 0.0203(13) 0.000
Cl2 0.0511(11) 0.0554(11) 0.0533(13) 0.0078(10) 0.0225(9) 0.0027(9)
Cl3 0.0708(13) 0.0484(11) 0.0605(13) -0.0053(10) 0.0276(10) 0.0026(10)
Cl1 0.0950(16) 0.0442(11) 0.0942(18) -0.0051(12) 0.0526(14) 0.0035(12)
N9 0.045(3) 0.047(4) 0.035(4) 0.004(3) 0.014(3) -0.003(3)
N7 0.040(3) 0.048(4) 0.036(3) 0.008(3) 0.004(3) -0.004(3)
O2 0.052(3) 0.098(4) 0.046(4) -0.012(4) 0.006(3) -0.016(3)
N3 0.041(3) 0.044(3) 0.044(4) 0.001(3) 0.008(3) -0.002(3)
O6 0.053(3) 0.093(4) 0.052(3) 0.008(3) 0.021(3) -0.019(3)
C4 0.044(4) 0.032(4) 0.050(5) 0.004(4) 0.011(4) 0.003(3)
N1 0.049(4) 0.049(4) 0.040(4) 0.004(3) 0.014(3) -0.005(3)
C8 0.050(4) 0.045(4) 0.033(4) -0.007(4) 0.000(3) -0.004(4)
C5 0.036(4) 0.039(4) 0.045(5) 0.003(4) 0.013(3) -0.001(3)
C2 0.047(4) 0.045(4) 0.039(5) -0.007(4) 0.011(4) -0.004(4)
C7 0.041(4) 0.072(5) 0.055(5) 0.013(5) 0.006(4) -0.015(4)
C1 0.061(5) 0.067(5) 0.043(5) -0.002(4) 0.017(4) -0.005(4)
C3 0.038(4) 0.090(6) 0.054(5) 0.003(5) 0.013(4) -0.018(4)
C6 0.045(4) 0.043(4) 0.038(4) 0.004(4) 0.008(3) -0.002(4)
`;

describe('cif-core read', () => {
    it('frame', async () => {
        const parsed = await CIF.parseText(cifCoreString).run();
        if (parsed.isError) return;
        const cifFile = parsed.result;
        const block = cifFile.blocks[0];

        expect(block.getField('cell_length_a')!.float(0)).toBe(11.0829);
        expect.assertions(1);
    });

    it('schema', async () => {
        const parsed = await CIF.parseText(cifCoreString).run();
        if (parsed.isError) return;
        const cifFile = parsed.result;
        const block = cifFile.blocks[0];
        const cifCore = CIF.schema.cifCore(block);

        expect(cifCore.cell.length_a.value(0)).toBe(11.0829);
        expect(cifCore.atom_site_aniso.U.value(0)).toEqual(new Float64Array([ 0.0425, 0, 0, 0.00089, 0.0423, 0, 0.01515, 0.00066, 0.0375 ]));
        expect.assertions(2);
    });
});