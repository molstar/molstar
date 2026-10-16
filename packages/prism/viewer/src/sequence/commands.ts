/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { CanvasMenuItem, CommandDefinition, ShellModel } from '@molstar/prism-ui';
import { combineLatest, map } from 'rxjs';
import { colorScheme } from './legend.js';
import type { SequenceViewer } from './model.js';
import { LetterSchemes } from './palette.js';
import type { TracksRail } from './presets.js';

// The sequence's commands (sequence.md §8): the one way to change what the sequence shows. The canvas stack's
// Sequence menu runs them, so does the command menu; `view.sequence` (show or hide) is the shell's own.

/** The overlay card the legend is in (registered by the viewer mode). */
export const LegendOverlay = 'sequence-legend';

const done = (message: string) => ({ message });

const Rails: readonly { id: string; label: string; rail: TracksRail }[] = [
  { id: 'ruler', label: 'Ruler', rail: { navigator: 'ruler' } },
  { id: 'strip', label: 'Strip', rail: { navigator: 'strip' } },
  { id: 'band', label: 'Band', rail: { navigator: 'band' } },
  { id: 'none', label: 'None', rail: { navigator: 'none' } },
];

const schemeId = (name: string | null) => `sequence.colors.${name === null ? 'none' : name.replace(':', '-')}`;
const schemeLabel = (name: string | null) => (name === null ? 'No colors' : (colorScheme(name)?.label ?? name));

export function sequenceCommands(
  viewer: SequenceViewer,
  shell: ShellModel,
  /** The viewer's API (the person's): style, colors and the rail change through it, as steps (sequence.md §7). */
  api?: { run(id: string, values?: Readonly<Record<string, unknown>>): Promise<unknown> },
): CommandDefinition<any>[] {
  const hide = () => {
    if (shell.dock('top').active.value) shell.toggle('top');
  };
  const show = (style: 'compact' | 'tracks') => {
    if (api) void api.run('sequence.style', { style: style === 'compact' ? 'text' : 'tracks' });
    else viewer.setStyle(style);
    shell.show('sequence');
  };
  const setRail = (rail: TracksRail) => {
    if (api)
      void api.run('sequence.rail', { navigator: rail.navigator, ...(rail.height ? { height: rail.height } : {}) });
    else viewer.setRail(rail);
  };
  const color = (letters: string | null) => {
    if (api) void api.run('sequence.colors', { scheme: letters ?? 'none' });
    else viewer.setTheme({ letters });
  };
  return [
    {
      id: 'sequence.none',
      title: 'No sequence',
      alias: ['hide sequence', 'sequence off'],
      category: 'Sequence',
      icon: 'ban',
      run: () => (hide(), done('Sequence hidden')),
    },
    {
      id: 'sequence.text',
      title: 'Sequence: text',
      alias: ['sequence compact', 'sequence letters'],
      category: 'Sequence',
      icon: 'align-left',
      run: () => (show('compact'), done('Sequence as text')),
    },
    {
      id: 'sequence.tracks',
      title: 'Sequence: tracks',
      alias: ['sequence zoom', 'sequence track'],
      category: 'Sequence',
      icon: 'chart-gantt',
      run: () => (show('tracks'), done('Sequence as tracks')),
    },
    ...Rails.map(
      (r): CommandDefinition<any> => ({
        id: `sequence.rail.${r.id}`,
        title: `Sequence tracks rail: ${r.label.toLowerCase()}`,
        category: 'Sequence',
        icon: 'ruler',
        run: () => (setRail(r.rail), show('tracks'), done(`Sequence rail: ${r.label.toLowerCase()}`)),
      }),
    ),
    ...LetterSchemes.map(
      (name): CommandDefinition<any> => ({
        id: schemeId(name),
        title: name === null ? 'Sequence: no colors' : `Sequence colors: ${schemeLabel(name)}`,
        category: 'Sequence',
        icon: 'palette',
        run: () => (color(name), done(`Sequence colors: ${schemeLabel(name).toLowerCase()}`)),
      }),
    ),
    {
      id: 'sequence.legend',
      title: 'Sequence legend…',
      alias: ['what the colors mean', 'sequence colors'],
      category: 'Sequence',
      icon: 'palette',
      run: () => (shell.openOverlay(LegendOverlay), done('')),
    },
  ];
}

/**
 * The Sequence button's menu on the canvas stack: one single choice of three, No sequence, Text, Tracks (the current one
 * in the accent, as the picking levels' list marks its own), then the letters' colors (also a single choice) and the
 * legend. The tracks' rail is the ruler's; the other rails are commands, for code and the command menu. Each item runs a
 * command by id (`open` is the command manager's).
 */
export function sequenceMenu(
  viewer: SequenceViewer,
  shell: ShellModel,
  open: (commandId: string) => void,
): CanvasMenuItem[] {
  const item = (id: string, label: string, run: string, extra: Partial<CanvasMenuItem> = {}): CanvasMenuItem => ({
    id,
    label,
    run: () => open(run),
    ...extra,
  });
  const shown = shell.dock('top').active.pipe(map((t) => !!t));
  const showing = (style: 'compact' | 'tracks') =>
    combineLatest([shown, viewer.style]).pipe(map(([on, s]) => on && s === style));
  return [
    item('none', 'No sequence', 'sequence.none', {
      icon: 'ban',
      choice: true,
      group: 'view',
      pressed: shown.pipe(map((on) => !on)),
    }),
    item('text', 'Text', 'sequence.text', {
      icon: 'align-left',
      choice: true,
      group: 'view',
      pressed: showing('compact'),
    }),
    item('tracks', 'Tracks', 'sequence.tracks', {
      icon: 'chart-gantt',
      choice: true,
      group: 'view',
      pressed: showing('tracks'),
    }),
    ...LetterSchemes.map((name) =>
      item(schemeId(name), schemeLabel(name), schemeId(name), {
        choice: true,
        group: 'colors',
        pressed: viewer.theme.pipe(map((t) => t.letters === name)),
      }),
    ),
    item('legend', 'Legend…', 'sequence.legend', { icon: 'palette', group: 'colors' }),
  ];
}
