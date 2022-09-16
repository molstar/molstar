/**
 * Copyright (c) 2017-2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Panagiotis Tourlas <panagiot_tourlov@hotmail.com>
 *
 * Adapted from MolQL project
 */

import { MolScriptBuilder } from '../../../mol-script/language/builder';
const B = MolScriptBuilder;
import { FunctionDict } from '../types';

export const functions: FunctionDict = {
    'sqr': {
        '@desc': 'square of x',
        '@examples': ['sqr(2)'],
        map: x => B.core.math.pow([x, 2]),
    },
    'sqrt': {
        '@desc': 'square root of x',
        '@examples': ['sqrt(2)'],
        map: x => B.core.math.sqrt([x]),
    },
    'abs': {
        '@desc': 'absolute value of x',
        '@examples': ['abs(2)'],
        map: x => B.core.math.abs([x]),
    },
    'floor': {
        '@desc': 'largest integer not greater than x',
        '@examples': ['floor(2)'],
        map: x => B.core.math.floor([x]),
    },
    'ceil': {
        '@desc': 'smallest integer not less than x',
        '@examples': ['ceil(2)'],
        map: x => B.core.math.ceil([x]),
    },
    'sin': {
        '@desc': 'sine of x',
        '@examples': ['sin(2)'],
        map: x => B.core.math.sin([x]),
    },
    'cos': {
        '@desc': 'cosine of x',
        '@examples': ['cos(2)'],
        map: x => B.core.math.cos([x]),
    },
    'tan': {
        '@desc': 'tangent of x',
        '@examples': ['tan(2)'],
        map: x => B.core.math.tan([x]),
    },
    'atan': {
        '@desc': 'arctangent of x',
        '@examples': ['atan(2)'],
        map: x => B.core.math.atan([x]),
    },
    'asin': {
        '@desc': 'arcsin of x',
        '@examples': ['asin(2)'],
        map: x => B.core.math.asin([x]),
    },
    'acos': {
        '@desc': 'arccos of x',
        '@examples': ['acos(2)'],
        map: x => B.core.math.acos([x]),
    },
    'sinh': {
        '@desc': 'hyperbolic sine of x',
        '@examples': ['sinh(2)'],
        map: x => B.core.math.sinh([x]),
    },
    'cosh': {
        '@desc': 'hyperbolic cosine of x',
        '@examples': ['cosh(2)'],
        map: x => B.core.math.cosh([x]),
    },
    'tanh': {
        '@desc': 'hyperbolic tangent of x',
        '@examples': ['tanh(2)'],
        map: x => B.core.math.tanh([x]),
    },
    'exp': {
        '@desc': 'e to the power x',
        '@examples': ['exp(2)'],
        map: x => B.core.math.exp([x]),
    },
    'log': {
        '@desc': 'natural log of x',
        '@examples': ['log(2)'],
        map: x => B.core.math.log([x]),
    },
    'log10': {
        '@desc': 'log base 10 of x',
        '@examples': ['log10(2)'],
        map: x => B.core.math.log10([x]),
    }
};
