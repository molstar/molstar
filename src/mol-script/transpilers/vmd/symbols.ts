/**
 * Copyright (c) 2017-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Panagiotis Tourlas <panagiot_tourlov@hotmail.com>
 *
 * Adapted from MolQL project
 */

import { properties } from './properties';
import { operators } from './operators';
import { keywords } from './keywords';
import { functions } from './functions';

export const Properties: string[] = [];
for (const name in properties) {
    if (properties[name].isUnsupported) continue;
    Properties.push(name);
    if (properties[name].abbr) Properties.push(...properties[name].abbr!);
}

export const Operators: string[] = [];
operators.forEach(o => {
    if (o.isUnsupported) return;
    Operators.push(o.name);
    if (o.abbr) Operators.push(...o.abbr);
});

export const Keywords: string[] = [];
for (const name in keywords) {
    if (!keywords[name].map) continue;
    Keywords.push(name);
    if (keywords[name].abbr) Keywords.push(...keywords[name].abbr!);
}

export const Functions: string[] = [];
for (const name in functions) {
    if (!functions[name].map) continue;
    Functions.push(name);
}

export const all = { Properties, Operators: [...Operators, ...Functions], Keywords };
