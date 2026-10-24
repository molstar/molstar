/**
 * Copyright (c) 2017-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Panagiotis Tourlas <panagiot_tourlov@hotmail.com>
 *
 * Adapted from MolQL project
 */

import { Expression } from '../language/expression.js';

/** What a transpiler is told about the names its text may use besides the language's own words. */
export interface TranspileOptions {
  /**
   * A name that is not a keyword, property or operator (PyMOL: an object or a named selection, `lig`, `/1abc//A`) as
   * MolQL; undefined when there is no such name. Without it, such names fail as before.
   */
  resolveName?(name: string): Expression | undefined;
}

export type Transpiler = (source: string, options?: TranspileOptions) => Expression;

export const Transpiler = (source: string) => Expression;
