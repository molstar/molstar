/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

/** Strings by code unit, as `<` orders them: the one string order of the queries' values. */
export const compareStrings = (a: string, b: string): number => (a < b ? -1 : a > b ? 1 : 0);
