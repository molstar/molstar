/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import REGL = require('regl');
import { InitializationOptions } from 'regl'

export function create(params: InitializationOptions) {
    return REGL(params)
}
