/*
 * Copyright (c) 2017 MolQL contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { transpiler as jmol } from './jmol/parser';
import { transpiler as pymol } from './pymol/parser';
import { transpiler as vmd } from './vmd/parser';
import { transpiler as rasmol } from './rasmol/parser';

export const _transpiler = {
    pymol,
    vmd,
    jmol,
    rasmol
};
