/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { KINFormatProvider } from '@molstar/kinemage-extension/model';
import { operationCommand } from '@molstar/prism-viewer/commands/operations';
import type { ViewerExtension, ViewerSkill } from '@molstar/prism-viewer/spec';
import { KinemageApi } from './operations.js';

declare module '@molstar/prism-viewer/api' {
  // The `kinemage` area: the groups, masters and views of a loaded kinemage (parity.md §6).
  interface ExtensionApis {
    kinemage: KinemageApi;
  }
}

// Kinemage (.kin, MolProbity's output) in the Prism viewer (parity.md §6, README decision 27): the plugin gets the format,
// so `load.files` of a .kin shows it; what the v5 right-panel controls did (groups, masters, views, animate) is
// operations with commands and tools. No view of its own: the controls were buttons and checkboxes over those.

export const KinemageSkill: ViewerSkill = {
  id: 'kinemage',
  title: 'Kinemage',
  summary: 'show and hide the groups and masters of a .kin file, apply its views, step its animation',
  guidance: [
    'Skill kinemage: a .kin file (MolProbity output) opens with load.files and shows as shapes, not as a structure.',
    'get_kinemage lists its groups (with subgroups), masters and views and whether each is shown;',
    'show_kinemage_part shows or hides one by kind and name; apply_kinemage_view moves the camera to a view;',
    'step_kinemage_animation shows the next frame of its animate groups.',
  ].join('\n'),
  tools: [
    {
      id: 'get_kinemage',
      operation: 'kinemage.get',
      description: 'Lists the loaded kinemages: groups, masters, views.',
    },
    {
      id: 'show_kinemage_part',
      operation: 'kinemage.visibility',
      description: 'Shows or hides a kinemage group, subgroup or master by name.',
    },
    {
      id: 'apply_kinemage_view',
      operation: 'kinemage.view',
      description: "Moves the camera to one of a kinemage's views.",
    },
    {
      id: 'step_kinemage_animation',
      operation: 'kinemage.animate',
      description: "Shows the next frame of a kinemage's animate groups.",
    },
  ],
};

export const KinemageExtension: ViewerExtension = {
  id: 'kinemage',
  title: 'Kinemage',
  // The v5 behavior's custom property is never applicable and its controls are React; only the format is needed.
  plugin: { registry: [{ formats: [KINFormatProvider] }], behaviors: [] },
  api: { kinemage: (ctx) => new KinemageApi(ctx) },
  commands: (api) => {
    const k = api.kinemage;
    if (!k) return [];
    return [
      operationCommand(api, k.visibility, {
        category: 'Kinemage',
        icon: 'eye',
        alias: ['kin', 'molprobity', 'master', 'group'],
        action: 'Apply',
      }),
      operationCommand(api, k.view, { category: 'Kinemage', icon: 'camera', alias: ['kin view', 'molprobity view'] }),
      operationCommand(api, k.animate, { category: 'Kinemage', icon: 'play', alias: ['kin animate', 'frames'] }),
    ];
  },
  skills: [KinemageSkill],
};
