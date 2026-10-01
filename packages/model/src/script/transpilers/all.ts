/**
 * Copyright (c) 2017-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 *
 * Adapted from MolQL project
 */

import { transpiler as jmol } from './jmol/parser.js';
import { transpiler as pymol } from './pymol/parser.js';
import { transpiler as vmd } from './vmd/parser.js';

export const _transpiler = {
    pymol,
    vmd,
    jmol,
};
