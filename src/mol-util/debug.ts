/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

/**
 * on node `process.env.NODE_ENV` is available, in webpack build it is automatically set
 * by the DefinePlugin to the webpack `mode` value
 */
let isProductionMode = process.env.NODE_ENV === 'production'

/**
 * set to true to enable more comprehensive checks and assertions,
 * mostly used in `mol-gl` and in valence-model calculation
 */
let isDebugMode = process.env.DEBUG === '*' || process.env.DEBUG === 'molstar'

if (typeof window !== 'undefined' && !(window as any).setMolStarDebugMode) {
    try {
        (window as any).setMolStarDebugMode = function setMolStarDebugMode(isDebug?: boolean, isProduction?: boolean) {
            if (typeof isDebug !== 'undefined') isDebugMode = isDebug;
            if (typeof isProduction !== 'undefined') isProductionMode = isProduction;
        };
    } catch { }
}

export { isProductionMode, isDebugMode }