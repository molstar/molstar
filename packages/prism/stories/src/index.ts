/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Story mode, opt-in: an app gets it only by depending on this package (nothing in @molstar/prism-ui imports it). A
// story is a document (schemas and MolViewSpec) the canvas renders; story mode is a mode of the app (design/modes.md)
// with its own tabs and commands in place of Mol*'s.

export * from './model.js';
export { StoryTab } from './views/story-tab.js';
