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
`

describe('cif-core read', () => {
    it('frame', async () => {
        const parsed = await CIF.parseText(cifCoreString).run()
        if (parsed.isError) return
        const cifFile = parsed.result
        const block = cifFile.blocks[0]

        expect(block.getField('cell_length_a')!.float(0)).toBe(11.0829)
        expect.assertions(1)
    });

    it('schema', async () => {
        const parsed = await CIF.parseText(cifCoreString).run()
        if (parsed.isError) return
        const cifFile = parsed.result
        const block = cifFile.blocks[0]
        const cifCore = CIF.schema.cifCore(block)

        expect(cifCore.cell.length_a.value(0)).toBe(11.0829)
        expect.assertions(1)
    });
});