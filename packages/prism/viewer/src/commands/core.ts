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
import { componentCommands } from './component.js';
import { operationCommand } from './operations.js';
import { representationCommands } from './representation.js';
import { selectionCommands } from './selection.js';
import { trajectoryCommands } from './trajectory.js';
import { volumeCommands } from './volume.js';

// The core's commands (operations.md §3), over the API bound to the person; the features bring theirs.

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
    operationCommand(api, api.structure.get, {
      category: 'Data',
      icon: 'list',
      alias: ['what is loaded', 'scene', 'outline'],
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
    operationCommand(api, api.structure.preset, {
      category: 'Data',
      icon: 'sparkles',
      title: 'Apply a representation preset…',
      alias: ['preset', 'representation preset', 'polymer cartoon', 'illustrative', 'auto'],
      action: 'Apply',
    }),
    operationCommand(api, api.structure.model, {
      category: 'Data',
      icon: 'layers',
      alias: ['model', 'next model', 'nmr model', 'switch model'],
      action: 'Show',
    }),
    operationCommand(api, api.structure.assembly, {
      category: 'Data',
      icon: 'box',
      alias: ['assembly', 'biological assembly', 'unit cell', 'deposited coordinates'],
      action: 'Show',
    }),
    operationCommand(api, api.structure.visibility, {
      category: 'Data',
      icon: 'eye',
      alias: ['show structure', 'hide structure', 'toggle structure'],
      action: 'Set',
    }),
    operationCommand(api, api.structure.remove, {
      category: 'Data',
      icon: 'trash-2',
      alias: ['delete structure', 'remove structure', 'close structure'],
      action: 'Remove',
    }),
    operationCommand(api, api.styles.apply, {
      category: 'Style',
      icon: 'swatch-book',
      alias: ['quick style', 'look', 'style', 'default style', 'illustrative style'],
      action: 'Apply',
    }),
    operationCommand(api, api.annotations.apply, {
      category: 'Style',
      icon: 'tag',
      title: 'Show an annotation…',
      alias: ['annotation', 'b-factor', 'secondary structure', 'color by annotation'],
      action: 'Show',
    }),
    operationCommand(api, api.camera.zoom, {
      category: 'View',
      icon: 'zoom-in',
      alias: ['zoom in', 'zoom out', 'closer', 'farther'],
      action: 'Zoom',
    }),
    operationCommand(api, api.camera.orient, {
      category: 'View',
      icon: 'rotate-3d',
      alias: ['orient', 'principal axes', 'align view'],
      action: 'Orient',
    }),
    operationCommand(api, api.camera.resetAxes, {
      category: 'View',
      icon: 'axis-3d',
      alias: ['reset axes', 'align axes', 'cartesian axes', 'xyz'],
      action: 'Align',
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
    operationCommand(api, api.structure.symmetry, {
      category: 'Data',
      icon: 'hexagon',
      alias: ['symmetry', 'symmetry mates', 'crystal', 'unit cell', 'lattice'],
    }),
    operationCommand(api, api.structure.unitCell, {
      category: 'Data',
      icon: 'box',
      alias: ['unit cell', 'crystal cell', 'cell box', 'lattice box'],
      action: 'Show',
    }),
    operationCommand(api, api.structure.transform, {
      category: 'Data',
      icon: 'move-3d',
      alias: ['transform', 'move', 'rotate', 'translate', 'conformation transform', 'matrix'],
      action: 'Move',
    }),
    operationCommand(api, api.structure.transformReset, {
      category: 'Data',
      icon: 'rotate-ccw',
      alias: ['reset transform', 'reset position', 'untransform'],
      action: 'Reset',
    }),
    ...componentCommands(api),
    ...representationCommands(api),
    ...selectionCommands(api),
    operationCommand(api, api.focus.again, {
      category: 'Selection',
      icon: 'history',
      title: 'Focus again…',
      alias: ['focus again', 'recent focus', 'focus history', 'previous focus', 'back to'],
      action: 'Focus',
    }),
    ...volumeCommands(api),
    operationCommand(api, api.volume.streamingView, {
      category: 'Volumes',
      icon: 'volume',
      title: 'Set the streamed density view…',
      alias: ['streaming view', 'density around focus', 'density detail', 'whole cell density', 'streaming detail'],
      action: 'Set',
      // The channels' iso values and colors are the streaming row's inspector's (and the tool's).
      values: { channels: {} },
    }),
    ...trajectoryCommands(api),
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
