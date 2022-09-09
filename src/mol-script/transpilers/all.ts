/**
 * Copyright (c) 2017-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 *
 * Adapted from MolQL project
 */

import { transpiler as jmol } from './jmol/parser';
import { transpiler as pymol } from './pymol/parser';
import { transpiler as vmd } from './vmd/parser';

export const _transpiler = {
    pymol,
    vmd,
    jmol,
};
