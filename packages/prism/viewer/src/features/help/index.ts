/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { ViewerExtension } from '../../spec.js';
import { helpShell } from './shell.js';

// The `help` feature (viewer-spec.md §3.1): Help, at the foot of the rail: the basics in prose, then every command
// with its keys and the canvas's mouse and key bindings; `help.show` (?) opens it.

export const HelpFeature: ViewerExtension = {
  id: 'help',
  title: 'Help',
  builtIn: true,
  shell: (ctx) => helpShell(ctx),
};
