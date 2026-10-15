/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// `@molstar/prism-viewer/api` (api.md §6): the API's types and helpers, and the module an extension augments with its
// areas: `declare module '@molstar/prism-viewer/api' { interface ExtensionApis { 'geo-export': GeoExportApi } }`.

export * from './api/index.js';

/** Extensions' areas, by id: an extension augments this interface here, so `viewer.api['geo-export']` is typed. */
export interface ExtensionApis {}
