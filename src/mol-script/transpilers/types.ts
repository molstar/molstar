/**
 * Copyright (c) 2017-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Panagiotis Tourlas <panagiot_tourlov@hotmail.com>
 *
 * Adapted from MolQL project
 */

import * as P from '../../mol-util/monadic-parser';
import { Expression } from '../language/expression';

export interface AtomGroupArgs {
    [index: string]: any
    'entity-test'?: Expression
    'chain-test'?: Expression
    'residue-test'?: Expression
    'atom-test'?: Expression
    'groupBy'?: Expression
}

export interface Keyword {
    '@desc': string
    abbr?: string[]
    map?: () => Expression /* not given means the keyword is unsupported */
}

export type KeywordDict = { [name: string]: Keyword }

export interface Property {
    '@desc': string
    '@examples': string[]
    isUnsupported?: boolean
    isNumeric?: boolean
    abbr?: string[]
    regex: RegExp
    map: (s: string) => any
    level: 'atom-test' | 'residue-test' | 'chain-test' | 'entity-test'
    property?: Expression
}

export type PropertyDict = { [name: string]: Property }

export interface Operator {
    '@desc': string
    '@examples': string[]
    name: string
    abbr?: string[]
    isUnsupported?: boolean
    type: (p1: P.MonadicParser<any>, p2: P.MonadicParser<any>, fn: any) => P.MonadicParser<any>
    rule: P.MonadicParser<any>
    map: (x: any, y: any, z?: any) => Expression | Expression[]
}

export type OperatorList = Operator[]

export interface Function {
    '@desc': string
    '@examples': string[]
    map?: (x: any) => Expression /* not given means the keyword is unsupported */
}

export type FunctionDict = { [name: string]: Function }
