/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type MockMode, ShellStory, type ShellStoryProps } from '@molstar/prism-ui/testing/index';
import type { MockScenePreset } from './scene.js';
import { mockViewerMode } from './viewer-mode.js';

/** Options of the viewer mode's mocks: the scene preset, remote states, what the data has to play. */
export interface ViewerMockOptions {
  preset?: MockScenePreset;
  /** The spec turns remote states (a state server) on; a number fills the server with that many. */
  remoteStates?: boolean | number;
  /** What the loaded data has to play: a session's snapshots, a trajectory's models (plain 1TQN has neither). */
  playback?: { snapshots?: boolean; trajectory?: boolean };
}

/** The mock frame with the viewer mode (and, given, other modes after it): the app as Mol*'s viewer. */
export function ViewerShellStory(props: ShellStoryProps & ViewerMockOptions & { more?: readonly MockMode[] }) {
  return <ShellStory {...props} modes={[mockViewerMode(props), ...(props.more ?? [])]} />;
}
