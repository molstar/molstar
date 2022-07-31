/*                                                                                                                                           
 * Copyright (c) 2017-2021 mol* contributors, licensed under MIT, See LICENSE file for more info.                                           
 * @author Alexander Rose <alexander.rose@weirdbyte.de>                                                                                     
 * @author Panagiotis Tourlas <panagiot_tourlov@hotmail.com>                                                                                 
 *                                                                                                                                           
 * @author Koya Sakuma                                                                                                                       
 * This module was taken from MolQL and modified in similar manner as pymol and vmd tranpilers.                                             \
*/

import { properties } from './properties';
import { operators } from './operators';
import { keywords } from './keywords';

export const Properties: string[] = []
for (const name in properties) {
    if (properties[name].isUnsupported) continue
    Properties.push(name)
    if (properties[name].abbr) Properties.push(...properties[name].abbr!)
}

export const Operators: string[] = []
operators.forEach(o => {
    if (o.isUnsupported) return
    Operators.push(o.name)
    if (o.abbr) Operators.push(...o.abbr)
})

export const Keywords: string[] = []
for (const name in keywords) {
    if (!keywords[name].map) continue
    Keywords.push(name)
    if (keywords[name].abbr) Keywords.push(...keywords[name].abbr!)
}

const _all = { Properties, Operators, Keywords }
export default _all
