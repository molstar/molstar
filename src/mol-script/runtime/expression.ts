/*
 * Copyright (c) 2018 Mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import Environment from './environment'

type RuntimeExpression<C = any, T = any> = (env: Environment<C>) => T

export interface ExpressionInfo {
    isConst?: boolean
}

namespace RuntimeExpression {
    export function constant<C, T>(c: T): RuntimeExpression<C, T> {
        return env => c;
    }

    export function func<C, T>(f: (env: Environment<C>) => T): RuntimeExpression<C, T> {
        return env => f(env);
    }
}

export default RuntimeExpression