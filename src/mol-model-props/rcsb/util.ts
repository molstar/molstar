/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Structure } from '../../mol-model/structure';
import { ParamDefinition as PD } from '../../mol-util/param-definition'
import {  AssemblySymmetryProvider } from './assembly-symmetry';

export function getSymmetrySelectParam(structure?: Structure) {
    const param = PD.Select<number>(0, [[0, 'No Symmetries']])
    if (structure) {
        const assemblySymmetry = AssemblySymmetryProvider.getValue(structure).value
        if (assemblySymmetry) {
            const options: [number, string][] = []
            for (let i = 0, il = assemblySymmetry.length; i < il; ++i) {
                const { symbol, kind } = assemblySymmetry[i]
                options.push([ i, `${i + 1}: ${symbol} ${kind}` ])
            }
            if (options.length) {
                param.options = options
                param.defaultValue = options[0][0]
            }
        }
    }
    return param
}