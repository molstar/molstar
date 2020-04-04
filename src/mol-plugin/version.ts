/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

/** version from package.json, to be filled in at bundle build time */
declare const __VERSION__: string
export const PLUGIN_VERSION = __VERSION__;

/** unix time stamp, to be filled in at bundle build time */
declare const __VERSION_TIMESTAMP__: number
export const PLUGIN_VERSION_TIMESTAMP = __VERSION_TIMESTAMP__;
export const PLUGIN_VERSION_DATE = new Date(PLUGIN_VERSION_TIMESTAMP);