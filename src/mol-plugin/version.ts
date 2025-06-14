/**
 * Copyright (c) 2018-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

declare const __MOLSTAR_BUILD_TIMESTAMP__: string | number | undefined;
declare const __MOLSTAR_PLUGIN_VERSION__: string | undefined;

/** version from package.json, to be filled in at build time */
export const PLUGIN_VERSION = typeof __MOLSTAR_PLUGIN_VERSION__ !== 'undefined' ? __MOLSTAR_PLUGIN_VERSION__ : '(development)';

/** to be filled in at build time */
export const PLUGIN_VERSION_DATE = new Date(typeof __MOLSTAR_BUILD_TIMESTAMP__ !== 'undefined' ? +__MOLSTAR_BUILD_TIMESTAMP__ : Date.now());