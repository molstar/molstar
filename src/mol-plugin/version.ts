/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

export const PLUGIN_VERSION = '0.5.0';
/** unix time stamp, to be filled in at bundle build time */
declare const __PLUGIN_VERSION_TIMESTAMP__: number
export const PLUGIN_VERSION_TIMESTAMP = __PLUGIN_VERSION_TIMESTAMP__;
export const PLUGIN_VERSION_DATE = new Date(PLUGIN_VERSION_TIMESTAMP);