/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { CommandDefinition } from '@molstar/prism-ui/core';
import type { ViewerApi } from '../api/area.js';
import type { ViewerCommandContext } from '../api/managers.js';
import { operationCommand } from './operations.js';
import { selectionCommands } from './selection.js';

// The core's commands (operations.md §3), over the API bound to the person; the features bring theirs (the
// structure, volume, trajectory and annotations features theirs: viewer-spec.md §3.1).

export function coreCommands(api: ViewerApi): CommandDefinition<ViewerCommandContext>[] {
  return [
    operationCommand(api, api.load.files, {
      category: 'Data',
      icon: 'folder-open',
      title: 'Open files…',
      alias: ['open', 'load file', 'import'],
      action: 'Open',
    }),
    operationCommand(api, api.session.open, {
      category: 'Data',
      icon: 'folder-open',
      title: 'Open a session…',
      alias: ['open session', 'molx', 'molj'],
      action: 'Open',
    }),
    operationCommand(api, api.session.clear, {
      category: 'Data',
      icon: 'trash-2',
      alias: ['clear', 'reset', 'remove everything', 'new session'],
    }),
    operationCommand(api, api.canvas.background, {
      category: 'Canvas',
      icon: 'palette',
      alias: ['background', 'bg'],
      action: 'Set',
    }),
    operationCommand(api, api.canvas.palette, {
      title: 'Use colorblind-safe colors',
      category: 'Canvas',
      icon: 'eye',
      alias: ['colorblind', 'color blind', 'okabe ito', 'palette', 'accessible colors'],
    }),
    operationCommand(api, api.camera.reset, {
      category: 'View',
      icon: 'focus',
      alias: ['reset view', 'fit', 'recenter'],
      values: { animate: true },
    }),
    operationCommand(api, api.camera.zoom, {
      category: 'View',
      icon: 'zoom-in',
      alias: ['zoom in', 'zoom out', 'closer', 'farther'],
      action: 'Zoom',
      values: { animate: true },
    }),
    operationCommand(api, api.camera.orient, {
      category: 'View',
      icon: 'rotate-3d',
      alias: ['orient', 'principal axes', 'align view'],
      action: 'Orient',
      values: { animate: true },
    }),
    operationCommand(api, api.camera.resetAxes, {
      category: 'View',
      icon: 'axis-3d',
      alias: ['reset axes', 'align axes', 'cartesian axes', 'xyz'],
      action: 'Align',
      values: { animate: true },
    }),
    operationCommand(api, api.camera.set, {
      category: 'View',
      icon: 'camera',
      alias: ['camera position', 'camera target', 'look at', 'set view', 'mvs camera'],
      action: 'Set',
    }),
    operationCommand(api, api.view.expand, {
      category: 'View',
      icon: 'maximize-2',
      alias: ['expand', 'expand viewport', 'hide panels', 'maximize', 'canvas only'],
      values: { on: 'toggle' },
    }),
    operationCommand(api, api.view.fullscreen, {
      category: 'View',
      icon: 'fullscreen',
      alias: ['fullscreen', 'full screen'],
      values: { on: 'toggle' },
    }),
    ...selectionCommands(api),
    operationCommand(api, api.focus.again, {
      category: 'Selection',
      icon: 'history',
      title: 'Focus again…',
      alias: ['focus again', 'recent focus', 'focus history', 'previous focus', 'back to'],
      action: 'Focus',
    }),
    operationCommand(api, api.history.undo, {
      category: 'Edit',
      icon: 'undo-2',
      alias: ['undo', 'back'],
      shortcut: { key: 'z', mod: true },
    }),
    operationCommand(api, api.history.redo, {
      category: 'Edit',
      icon: 'redo-2',
      alias: ['redo', 'forward'],
      shortcut: { key: 'z', mod: true, shift: true },
    }),
    operationCommand(api, api.history.undoTurn, {
      category: 'Edit',
      icon: 'history',
      alias: ['undo turn', 'undo assistant', 'revert assistant'],
      action: 'Undo',
    }),
  ];
}
