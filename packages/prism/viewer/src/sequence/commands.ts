/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import type { CanvasMenuItem, CommandDefinition, ShellModel } from '@molstar/prism-ui';
import { combineLatest, map, of } from 'rxjs';
import { colorScheme } from './legend.js';
import type { SequenceTrack, SequenceViewer } from './model.js';
import { LetterSchemes } from './palette.js';
import type { TracksRail } from './presets.js';

import { derive } from '@molstar/prism-ui/reactive/state';
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

/** The viewer's API as the person: what changes the sequence's state goes through it, as steps (sequence.md §7). */
export interface SequenceApiRunner {
  run(id: string, values?: Readonly<Record<string, unknown>>): Promise<unknown>;
}

/**
 * Hides or shows a track (sequence.md §8): `sequence.track.hide` / `.show` through the API as the person (one step), or
 * on the viewer where there is no API (tests, a bare model).
 */
export function setTrackHidden(
  viewer: SequenceViewer,
  api: SequenceApiRunner | undefined,
  track: string,
  hidden: boolean,
) {
  if (api) return api.run(hidden ? 'sequence.track.hide' : 'sequence.track.show', { track });
  return viewer.setTrackHidden(track, hidden);
}

/** Every hidden track shown: `sequence.track.show` with no track. */
export function showAllTracks(viewer: SequenceViewer, api: SequenceApiRunner | undefined) {
  if (api) return api.run('sequence.track.show', {});
  return viewer.setTracksHidden(
    viewer.tracks.value.filter((t) => t.hidden).map((t) => t.ref),
    false,
  );
}

/** A choice of the tracks that `keep` (by ref, labelled by name), for a command's form. */
function trackChoice(tracks: readonly SequenceTrack[], keep: (t: SequenceTrack) => boolean) {
  const options = tracks.filter(keep).map((t): [string, string] => [t.ref, t.label]);
  return { track: PD.Select(options[0]?.[0] ?? '', options, { label: 'Track' }) };
}

export function sequenceCommands(
  viewer: SequenceViewer,
  shell: ShellModel,
  /** The viewer's API (the person's): style, colors and the rail change through it, as steps (sequence.md §7). */
  api?: SequenceApiRunner,
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
      id: 'sequence.track.hide',
      title: 'Hide a sequence track',
      alias: ['hide track', 'sequence hide track'],
      category: 'Sequence',
      icon: 'eye-off',
      when: () => viewer.tracks.value.some((t) => !t.hidden) || 'No sequence track is shown.',
      params: () => trackChoice(viewer.tracks.value, (t) => !t.hidden),
      run: async (args) => {
        const track = String(args.track ?? '');
        const label = viewer.tracks.value.find((t) => t.ref === track)?.label ?? track;
        await setTrackHidden(viewer, api, track, true);
        return { ...done(`Track hidden: ${label}`), undoable: true };
      },
    },
    {
      id: 'sequence.track.show',
      title: 'Show a hidden sequence track',
      alias: ['show track', 'sequence show track'],
      category: 'Sequence',
      icon: 'eye',
      when: () => viewer.tracks.value.some((t) => t.hidden) || 'No sequence track is hidden.',
      params: () => trackChoice(viewer.tracks.value, (t) => t.hidden),
      run: async (args) => {
        const track = String(args.track ?? '');
        const label = viewer.tracks.value.find((t) => t.ref === track)?.label ?? track;
        await setTrackHidden(viewer, api, track, false);
        return { ...done(`Track shown: ${label}`), undoable: true };
      },
    },
    {
      id: 'sequence.track.show-all',
      title: 'Show all sequence tracks',
      alias: ['show all tracks', 'show hidden tracks'],
      category: 'Sequence',
      icon: 'eye',
      when: () => viewer.tracks.value.some((t) => t.hidden) || 'No sequence track is hidden.',
      run: async () => {
        await showAllTracks(viewer, api);
        return { ...done('All tracks shown'), undoable: true };
      },
    },
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
  /** The API as the person: the Visible tracks submenu's rows hide and show through it (steps). */
  api?: SequenceApiRunner,
): CanvasMenuItem[] {
  const item = (id: string, label: string, run: string, extra: Partial<CanvasMenuItem> = {}): CanvasMenuItem => ({
    id,
    label,
    run: () => open(run),
    ...extra,
  });
  const shown = derive(shell.dock('top').active, (t) => !!t);
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
        pressed: derive(viewer.theme, (t) => t.letters === name),
      }),
    ),
    item('legend', 'Legend…', 'sequence.legend', { icon: 'palette', group: 'colors' }),
    {
      id: 'track-list',
      // Not "Tracks": the style above is called that.
      label: 'Visible tracks',
      icon: 'layers',
      group: 'tracks',
      run: () => {},
      items: derive(viewer.tracks, (tracks) => trackRows(viewer, tracks, api)),
    },
  ];
}

/**
 * The Visible tracks submenu (sequence.md §8): a checkbox per track of the section shown (checked: shown), then Show all, there
 * only while one is hidden.
 */
function trackRows(
  viewer: SequenceViewer,
  tracks: readonly SequenceTrack[],
  api?: SequenceApiRunner,
): CanvasMenuItem[] {
  const rows = tracks.map(
    (t): CanvasMenuItem => ({
      id: `track:${t.ref}`,
      label: t.label,
      group: 'tracks',
      pressed: of(!t.hidden),
      run: () => void setTrackHidden(viewer, api, t.ref, !t.hidden),
    }),
  );
  if (!tracks.some((t) => t.hidden)) return rows;
  return [
    ...rows,
    { id: 'show-all', label: 'Show all', icon: 'eye', group: 'all', run: () => void showAllTracks(viewer, api) },
  ];
}
