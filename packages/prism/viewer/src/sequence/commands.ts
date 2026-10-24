/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import type { CanvasMenuItem, CommandDefinition, CommandGroup, ShellModel } from '@molstar/prism-ui';
import type { SequenceTrackProvider } from './tracks.js';
import { combineLatest, map } from 'rxjs';
import { colorScheme } from './legend.js';
import type { SequenceTrack, SequenceViewer } from './model.js';
import { LetterSchemes } from './palette.js';
import type { TracksRail } from './presets.js';

import { derive } from '@molstar/prism-ui/reactive/state';
// The sequence's commands (sequence.md §8): the one way to change what the sequence shows. The canvas stack's
// Sequence menu runs them, so does the command menu; `view.sequence` (show or hide) is the shell's own. Which
// structure, its tracks and its chains are the structure tree's Sequence node (§8.1); these commands do the same for
// the structure shown.

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

/** The command menu's group of Add track (sequence.md §8.1): what the Sequence node's inspector opens. */
export const AddTrackGroup: CommandGroup = { id: 'sequence-tracks', title: 'Add track', icon: 'plus', order: 5 };

/**
 * A command a track provider (`sequence.track.add`, on the structure the band shows), in Add track's group: there while
 * it applies and is out; its result said as the API says it.
 */
export function addTrackCommands(
  viewer: SequenceViewer,
  providers: readonly SequenceTrackProvider[],
  api?: SequenceApiRunner,
): CommandDefinition<any>[] {
  return providers.map((p) => ({
    id: `sequence.track.add.${p.id}`,
    title: p.label,
    alias: ['add track', `add ${p.label.toLowerCase()} track`],
    category: 'Sequence tracks',
    group: AddTrackGroup.id,
    icon: 'plus',
    when: () => {
      const shown = viewer.structure.value;
      if (!shown) return 'No sequence is shown.';
      const out = viewer.addableTracks({ ref: shown.ref, structure: shown.structure }).some((a) => a.id === p.id);
      if (out) return true;
      const there = viewer.tracks.value.some((t) => t.ref === p.track);
      return there ? 'Already shown' : 'Not for this structure';
    },
    run: async () => {
      const r = await ran(api ? api.run('sequence.track.add', { track: p.id }) : viewer.setTrackIn(p.id, true));
      return { message: r ?? `Track added: ${p.label}`, undoable: true };
    },
  }));
}

/** What a run said when it failed (an API result not ok); undefined: it went. */
async function ran(run: Promise<unknown>): Promise<string | undefined> {
  const r = (await run) as { ok?: boolean; message?: string } | boolean | undefined;
  return typeof r === 'object' && r?.ok === false ? (r.message ?? 'It did not go') : undefined;
}

/** The chain filter's presets (sequence.md §8): every chain, the polymers only, every copy apart. */
const ChainPresets = [
  { id: 'all', title: 'Sequence: every chain', values: {}, done: 'Every chain' },
  { id: 'polymer', title: 'Sequence: polymer chains only', values: { polymer: true }, done: 'Polymer chains only' },
  { id: 'each', title: 'Sequence: every copy apart', values: { each: true }, done: 'Every copy apart' },
] as const;

export function sequenceCommands(
  viewer: SequenceViewer,
  shell: ShellModel,
  /** The viewer's API (the person's): style, colors and the rail change through it, as steps (sequence.md §7). */
  api?: SequenceApiRunner,
  /** Opens the structure tree at the Sequence node of the structure shown (sequence.md §8.1); none: no tree. */
  reveal?: () => void,
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
  const addable = () => {
    const shown = viewer.structure.value;
    return shown ? viewer.addableTracks({ ref: shown.ref, structure: shown.structure }) : [];
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
      id: 'sequence.track.add',
      title: 'Add a sequence track',
      alias: ['add track', 'sequence add track'],
      category: 'Sequence',
      icon: 'plus',
      when: () => addable().length > 0 || 'No other track for this structure.',
      params: () => {
        const options = addable().map((a): [string, string] => [a.id, a.label]);
        return { track: PD.Select(options[0]?.[0] ?? '', options, { label: 'Track' }) };
      },
      run: async (args) => {
        const track = String(args.track ?? '');
        const label = addable().find((a) => a.id === track)?.label ?? track;
        const r = await ran(api ? api.run('sequence.track.add', { track }) : viewer.setTrackIn(track, true));
        return { ...done(r ?? `Track added: ${label}`), undoable: true };
      },
    },
    {
      id: 'sequence.track.remove',
      title: 'Remove a sequence track',
      alias: ['remove track', 'sequence remove track'],
      category: 'Sequence',
      icon: 'trash-2',
      when: () => viewer.tracks.value.some((t) => t.removable) || 'No sequence track can be removed.',
      params: () => trackChoice(viewer.tracks.value, (t) => t.removable),
      run: async (args) => {
        const track = String(args.track ?? '');
        const label = viewer.tracks.value.find((t) => t.ref === track)?.label ?? track;
        const r = await ran(api ? api.run('sequence.track.remove', { track }) : viewer.setTrackIn(track, false));
        return { ...done(r ?? `Track removed: ${label}`), undoable: true };
      },
    },
    ...ChainPresets.map(
      (c): CommandDefinition<any> => ({
        id: `sequence.chains.${c.id}`,
        title: c.title,
        category: 'Sequence',
        icon: 'link',
        when: () => !!viewer.structure.value || 'No sequence is shown.',
        run: async () => {
          const r = await ran(api ? api.run('sequence.chains', c.values) : viewer.setChains(c.values));
          return { ...done(r ?? c.done), undoable: true };
        },
      }),
    ),
    ...(reveal
      ? [
          {
            id: 'sequence.reveal',
            title: 'Sequence: edit in structure tree',
            alias: ['sequence tracks', 'sequence chains', 'sequence structure'],
            category: 'Sequence',
            icon: 'list-tree',
            when: () => !!viewer.structure.value || 'No sequence is shown.',
            run: () => (reveal(), done('')),
          } satisfies CommandDefinition<any>,
        ]
      : []),
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
 * in the accent, as the picking levels' list marks its own), then Edit in structure tree, where its structure, tracks,
 * chains, letter colors and legend are (sequence.md §8.1; nothing else here).
 * The tracks' rail is the ruler's; the other rails are commands, for code and the command menu. Each item runs a command
 * by id (`open` is the command manager's).
 */
export function sequenceMenu(
  viewer: SequenceViewer,
  shell: ShellModel,
  open: (commandId: string) => void,
  /** The structure tree is there to edit in (`sequence.reveal` registered). */
  revealable = false,
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
    ...(revealable
      ? [item('reveal', 'Edit in structure tree', 'sequence.reveal', { icon: 'list-tree', group: 'tree' })]
      : []),
  ];
}
