/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// In-memory fakes of the view models, for stories and tests here and in packages built on prism-ui; the whole mock
// shell (`ShellStory`) takes a module's contributions as an add-on.

export { story } from './story-root.js';
export { CanvasPlaceholder } from '../views/shell/canvas-placeholder.js';
export { mockSnapshots, mockTrajectory, mockTransition } from './playback.js';
export {
  type MockFrameContext,
  type MockMode,
  type MockModeParts,
  mockBlankMode,
  ShellStory,
  type ShellStoryProps,
} from './shell.js';
