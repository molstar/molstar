/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The stories, opt-in: an app gets them only by depending on this package (nothing in @molstar/prism-ui imports it). A
// story is a document (schemas and MolViewSpec) the canvas renders. No modes: the stories are an app of their own on
// prism-ui (the Prism viewer's `_spec/README.md`, decision 3), with their own tabs and commands in place of Mol*'s.

export * from './model.js';
export { StoryTab } from './views/story-tab.js';
