/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

/**
 * on node `process.env.NODE_ENV` is available, in webpack build it is automatically set
 * by the DefinePlugin to the webpack `mode` value
 */
const isProductionMode = process.env.NODE_ENV === 'production'

/**
 * set to true to enable more comprehensive checks and assertions,
 * mostly used in `mol-gl`
 */
const isDebugMode = process.env.DEBUG === '*' || process.env.DEBUG === 'molstar'

export { isProductionMode, isDebugMode }