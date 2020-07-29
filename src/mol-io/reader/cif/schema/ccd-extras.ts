/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { CCD_Schema } from './ccd';

// a reduced chem_comp_atom schema that provides charge and stereo_config information
export const ccd_chemCompAtom_schema = {
    comp_id: CCD_Schema.chem_comp_atom.comp_id,
    atom_id: CCD_Schema.chem_comp_atom.atom_id,
    charge: CCD_Schema.chem_comp_atom.charge,
    pdbx_stereo_config: CCD_Schema.chem_comp_atom.pdbx_stereo_config
};