/**
 * Copyright (c) 2017-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Panagiotis Tourlas <panagiot_tourlov@hotmail.com>
 *
 * Adapted from MolQL project
 */

import { Expression } from '../language/expression';

export type Transpiler = (source: string) => Expression

export const Transpiler = (source: string) => Expression;
