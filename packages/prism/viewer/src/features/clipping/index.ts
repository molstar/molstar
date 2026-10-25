/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Ludovic Autin <ludovic.autin@gmail.com>
 */

import { operationCommand } from '../../commands/operations.js';
import type { ViewerExtension, ViewerSkill } from '../../spec.js';
import { ClippingApi } from './operations.js';

declare module '../../api/area.js' {
  // The feature's area (typed as present, as in the default spec: api/area.ts `BuiltInApis`): cuts through
  // representations, read from and written to the state tree (clipping.md §2).
  interface BuiltInApis {
    clipping: ClippingApi;
  }
}

// The `clipping` feature (viewer-spec.md §3.1, clipping.md): the `clipping` area, its commands and skill; its views
// (the Selection tab's Clipping section, Clip… in the selection menu, Clip this… in the scene tree's row menus, the
// gizmos) are a chunk of their own. In the default spec.

export const ClippingSectionId = 'clipping.section.clip';

export const ClippingSkill: ViewerSkill = {
  id: 'clipping',
  title: 'Clipping',
  summary: 'cut through representations with a plane, a sphere, a box, a cylinder or a cone to show their inside',
  guidance: [
    'Skill clipping: a cut removes part of the representations it is on, to show what they hide.',
    '- clipping.add: a shape (plane removes the half its normal points to, facing the viewer at first; sphere, cube,',
    '  cylinder remove their inside; cone its inside) on the representations drawing a target (default: the selection),',
    '  or a node selector (volume > representation) for volumes and shapes; types narrows them (molecular-surface).',
    '- clipping.list gives the cuts by number; clipping.update moves (position), turns (normal or rotation), sizes',
    '  (scale) or inverts one; clipping.apply turns one on or off for a representation (its ref, as the list gives it);',
    '  clipping.keep leaves atoms (a target; the selection by default) whole in a cut, clipping.release cuts them again;',
    '  clipping.show draws a cut in the scene (saved with it); clipping.remove and clipping.clear take them away.',
    '- clipping.settings: variant instance leaves out whole copies whose center is cut away (pixel cuts through), for',
    '  every cut; solidSurface closes cut surfaces with a cap.',
  ].join('\n'),
  tools: [
    {
      id: 'list_cuts',
      operation: 'clipping.list',
      description: 'The cuts: their numbers, shapes, where, and on what.',
    },
    {
      id: 'add_cut',
      operation: 'clipping.add',
      description: 'Cuts representations with a shape to show their inside.',
    },
    { id: 'change_cut', operation: 'clipping.update', description: 'Moves, turns, sizes or inverts a cut.' },
    {
      id: 'apply_cut',
      operation: 'clipping.apply',
      description: 'Turns a cut on or off for one representation (leaves it whole).',
    },
    { id: 'remove_cut', operation: 'clipping.remove', description: 'Removes a cut.' },
    {
      id: 'cut_settings',
      operation: 'clipping.settings',
      description: 'How every cut cuts: through (pixel) or whole copies (instance); solid caps.',
    },
  ],
};

export const ClippingFeature: ViewerExtension = {
  id: 'clipping',
  title: 'Clipping',
  builtIn: true,
  api: { clipping: (ctx) => new ClippingApi(ctx) },
  commands: (api) => {
    const x = api.clipping;
    const category = 'Clipping';
    return [
      operationCommand(api, x.add, {
        category,
        icon: 'crop',
        title: 'Cut…',
        alias: ['clip', 'clipping', 'cut away', 'slice', 'cross section', 'clip plane', 'show inside'],
        action: 'Cut',
      }),
      operationCommand(api, x.update, { category, icon: 'crop', title: 'Change a cut…', action: 'Change' }),
      operationCommand(api, x.apply, {
        category,
        icon: 'scissors',
        title: 'Cut a representation, or leave it whole…',
        alias: ['keep whole', 'unclip', 'clip representation'],
        action: 'Apply',
      }),
      operationCommand(api, x.keep, {
        category,
        icon: 'shield',
        title: 'Keep atoms whole in a cut…',
        alias: ['keep whole', 'exclude from clipping', 'do not clip', 'protect from cut'],
        action: 'Keep',
      }),
      operationCommand(api, x.release, {
        category,
        icon: 'scissors',
        title: 'Cut atoms kept whole again…',
        alias: ['release', 'clip again'],
        action: 'Release',
      }),
      operationCommand(api, x.show, {
        category,
        icon: 'scan',
        title: 'Draw a cut in the scene…',
        alias: ['show cut', 'clip volume', 'draw clipping', 'show clipping'],
        action: 'Draw',
      }),
      operationCommand(api, x.settings, {
        category,
        icon: 'sliders-horizontal',
        title: 'Set how cuts cut…',
        alias: ['clip variant', 'clip instance', 'clip pixel', 'solid surface', 'solid cap', 'cut settings'],
        action: 'Set',
      }),
      operationCommand(api, x.remove, { category, icon: 'trash-2', alias: ['remove clip', 'uncut'], action: 'Remove' }),
      operationCommand(api, x.clear, { category, icon: 'trash-2', alias: ['remove all clipping'], action: 'Remove' }),
    ];
  },
  skills: [ClippingSkill],
  ui: async (ctx) => (await import('./ui.js')).clippingUI(ctx),
};
