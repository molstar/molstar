/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { operationCommand } from '../../commands/operations.js';
import type { ViewerExtension, ViewerSkill } from '../../spec.js';
import { ShapeApi } from './operations.js';

declare module '../../api/area.js' {
  // The feature's area (typed as present, as in the default spec: api/area.ts `BuiltInApis`): the shapes loaded on
  // their own (a mesh file), as the Structure panel's Shapes group shows them.
  interface BuiltInApis {
    shape: ShapeApi;
  }
}

// The `shapes` feature (viewer-spec.md §3.1): the `shape` area and its commands. The Structure panel shows the shapes
// (its Shapes group, `sources/scene.ts`); their look is `look.set`'s. In the default spec.

export const ShapesSkill: ViewerSkill = {
  id: 'shapes',
  title: 'Shapes',
  summary: 'meshes loaded from files (OBJ, PLY): list, show, hide, focus, remove',
  guidance: [
    'Skill shapes: a mesh file loads as a shape with its representations (shape.list gives their refs).',
    '- show_shape shows or hides one (node: its ref) or all; focus_shape points the camera at it; remove_shape removes',
    "  it with its file's data. Its look (material, light, quality) is set_look with node: the shape's ref.",
  ].join('\n'),
  tools: [
    { id: 'list_shapes', operation: 'shape.list', description: 'The shapes loaded, with their representations.' },
    { id: 'show_shape', operation: 'shape.visibility', description: 'Shows or hides a shape, or every shape.' },
    { id: 'focus_shape', operation: 'shape.focus', description: 'Points the camera at a shape.' },
    { id: 'remove_shape', operation: 'shape.remove', description: 'Removes a shape.' },
  ],
};

export const ShapesFeature: ViewerExtension = {
  id: 'shapes',
  title: 'Shapes',
  builtIn: true,
  api: { shape: (ctx) => new ShapeApi(ctx) },
  skills: [ShapesSkill],
  commands: (api) => {
    const x = api.shape;
    const category = 'Shapes';
    return [
      operationCommand(api, x.list, { category, icon: 'list', alias: ['shapes', 'meshes', 'list shapes'] }),
      operationCommand(api, x.visibility, {
        category,
        icon: 'eye',
        alias: ['show shape', 'hide shape', 'show mesh', 'hide mesh'],
        action: 'Set',
      }),
      operationCommand(api, x.focus, {
        category,
        icon: 'focus',
        alias: ['zoom to shape', 'center mesh', 'frame shape'],
        action: 'Focus',
      }),
      operationCommand(api, x.remove, {
        category,
        icon: 'trash-2',
        alias: ['delete shape', 'remove mesh'],
        action: 'Remove',
      }),
    ];
  },
};
