/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { CommandDefinition, CommandGroup } from '@molstar/prism-ui/core';
import type { ViewerCommandContext } from '../api/managers.js';
import type { PresetsModel, QuickStylesModel } from '../model/scene.js';

// The Structure tab's Preset and Styles in the command menu (shell.md, the command menu): a group each, a command per
// representation preset and per quick style, so the tab's buttons open the menu inside their group (teaching it as they
// go) and a search anywhere finds them. Each runs what the tab's own lists ran (`structure.preset`, `styles.apply`).

/** The groups the Structure tab's buttons open. */
export const StructureCommandGroups: readonly CommandGroup[] = [
  { id: 'preset', title: 'Preset', icon: 'shapes', order: -2 },
  { id: 'styles', title: 'Styles', icon: 'paintbrush', order: -1 },
];

export function structureCommands(
  presets: PresetsModel,
  styles: QuickStylesModel,
): CommandDefinition<ViewerCommandContext>[] {
  const done = { message: '' };
  return [
    ...presets.groups.flatMap((g) =>
      g.presets.map(
        (p): CommandDefinition<ViewerCommandContext> => ({
          id: `core.preset.${p.id}`,
          title: p.label,
          category: g.label ?? 'Presets',
          group: 'preset',
          icon: 'shapes',
          alias: ['preset', p.label.toLowerCase()],
          run: () => {
            presets.apply(p.id);
            return done;
          },
        }),
      ),
    ),
    ...styles.representations.map(
      (s): CommandDefinition<ViewerCommandContext> => ({
        id: `core.styles.${s.id}`,
        title: s.label,
        category: 'Representation',
        group: 'styles',
        icon: 'paintbrush',
        alias: ['style', 'show as'],
        run: () => {
          styles.applyRepresentation(s.id);
          return done;
        },
      }),
    ),
    ...styles.styles.map(
      (s): CommandDefinition<ViewerCommandContext> => ({
        id: `core.styles.${s.id}`,
        title: s.label,
        category: 'Look',
        group: 'styles',
        icon: 'swatch-book',
        alias: ['style', 'look'],
        run: () => {
          styles.applyStyle(s.id);
          return done;
        },
      }),
    ),
  ];
}
