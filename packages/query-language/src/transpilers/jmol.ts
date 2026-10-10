/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 */

import { registerTranspiler } from '../transpile.js';
import { transpiler } from './jmol/parser.js';

registerTranspiler('jmol', transpiler);
