/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Structure } from '@molstar/model/model/structure';
import { observableState, type ReadonlyState } from '@molstar/prism-ui/core';
import type { ViewerApi } from '../api/area.js';
import type { SceneSequence, SceneSequenceRow } from '../model/scene.js';
import { describeFilter, sequenceChains } from '../sequence/chain-filter.js';
import type { SequenceViewer, SequenceViewers } from '../sequence/model.js';
import type { SequenceSource } from './sequence.js';
import type { SourceScope } from './state.js';
import { LetterSchemes } from '../sequence/palette.js';

/** The letter schemes' names for people (as the legend names them; no Seq* module here). */
const SchemeNames: Readonly<Record<string, string>> = {
  'prism:residues': 'Prism residues',
  zappo: 'Zappo',
  taylor: 'Taylor',
  hydrophobicity: 'Hydrophobicity',
};

// The structure tree's side of the sequence (sequence.md §8.1): each structure's Sequence row read from the sequence
// model (its tracks, what can be added, its chains), and every change through the API as the person, so each is a
// step (or view state: what the band follows). Before the sequence's runtime is in, a row has no tracks yet.

export function createSceneSequence(
  scope: SourceScope,
  source: SequenceSource,
  viewers: ReadonlyState<SequenceViewers | undefined>,
  person: ViewerApi,
  /** Opens the legend card (the shell's); none: the row's menu has no Legend. */
  legend?: () => void,
): SceneSequence {
  const version = observableState(0);
  // Once a burst (an Undo, a switch of structure tells several of these): one pass of the tree a task.
  let queued = false;
  const bump = () => {
    if (queued) return;
    queued = true;
    queueMicrotask(() => {
      queued = false;
      version.set(version.value + 1);
    });
  };
  // What a row reads: what the band follows and shows, and the model's tracks, edits and providers.
  let modelSubs: { unsubscribe(): void }[] = [];
  let bandSub: { unsubscribe(): void } | undefined;
  const follow = (v: SequenceViewer | undefined) => {
    for (const s of modelSubs) s.unsubscribe();
    modelSubs = v
      ? [
          v.tracks.subscribe(bump),
          v.edits.subscribe(bump),
          v.providers.subscribe(bump),
          v.theme.subscribe(bump),
          v.style.subscribe(bump),
          // The letters' swatches follow the ink (the legend is made again with it).
          v.legend.subscribe(bump),
        ]
      : [];
    bump();
  };
  scope.add(viewers.subscribe((vs) => follow(vs?.get())));
  // Which structure it follows, not each frame of it (a trajectory's playback makes a new Structure a frame).
  let followedRef: string | undefined;
  scope.add(
    source.structure.subscribe((h) => {
      if (h?.ref === followedRef) return;
      followedRef = h?.ref;
      bump();
    }),
  );
  scope.add(source.followed.subscribe(bump));
  scope.add(source.shows.subscribe(bump));
  scope.add(
    source.band.subscribe((b) => {
      bandSub?.unsubscribe();
      bandSub = b?.shown.subscribe(bump);
      bump();
    }),
  );
  scope.add(() => {
    for (const s of modelSubs) s.unsubscribe();
    bandSub?.unsubscribe();
  });

  const viewer = () => viewers.value?.get();
  const shownRef = () => source.structure.value?.ref;
  const bandOn = () => source.band.value?.shown.value ?? true;
  const run = (id: string, values: Record<string, unknown>) => void person.run(id, values);

  const row = (structure: string, data: Structure): SceneSequenceRow | undefined => {
    if (data.polymerResidueCount === 0) return undefined;
    const v = viewer();
    const isShown = shownRef() === structure;
    const shown = isShown && bandOn();
    const chosen = isShown && (source.followed.value === structure || source.shows.value === structure);
    const of = { ref: structure, structure: data };
    const chains = sequenceChains(data);
    const filter = v?.chainsOf(structure);
    const tracks = (v?.tracksOf(of) ?? []).map((t) => ({
      ...t,
      // The colorings act on the structure the band shows.
      actions: isShown && v ? v.trackActions(t.ref) : [],
    }));
    // How many tracks (the style is the band's, said in its menu), whether the rule picked it, its chains.
    const hidden = tracks.filter((t) => t.hidden).length;
    // No count before the sequence's runtime is in: its tracks are not known yet (not none).
    const count = v
      ? `${tracks.length} track${tracks.length === 1 ? '' : 's'}${hidden ? `, ${hidden} hidden` : ''}`
      : undefined;
    const description = [
      shown ? undefined : 'Not shown',
      count,
      shown && !chosen ? 'auto' : undefined,
      describeFilter(filter, chains),
    ]
      .filter(Boolean)
      .join(' · ');
    return {
      shown,
      chosen,
      description,
      tracks,
      addable: v?.addableTracks(of).map((a) => ({ id: a.id, label: a.label })) ?? [],
      style: v?.style.value === 'compact' ? 'text' : 'tracks',
      colors: v ? v.theme.value.letters : 'prism:residues',
      chains,
      filter,
    };
  };

  return {
    version: version.asReadonly,
    row,
    toggle: (structure) => {
      const band = source.band.value;
      if (shownRef() === structure) {
        // The band shows it: the eye hides the band (No sequence); hidden, shows it again.
        band?.show(!bandOn());
        return;
      }
      run('sequence.follow', { structure });
      band?.show(true);
    },
    followAutomatically: () => run('sequence.follow', {}),
    setTrackHidden: (structure, _data, track, hidden) =>
      run(hidden ? 'sequence.track.hide' : 'sequence.track.show', { track, structure }),
    addTrack: (structure, _data, provider) => run('sequence.track.add', { track: provider, structure }),
    removeTrack: (structure, _data, track) => run('sequence.track.remove', { track, structure }),
    setChains: (structure, _data, filter) =>
      run('sequence.chains', {
        structure,
        ...(filter?.polymer ? { polymer: true } : {}),
        ...(filter?.each ? { each: true } : {}),
        ...(filter?.chains !== undefined ? { chains: filter.chains.length ? filter.chains.join(',') : 'none' } : {}),
      }),
    ...(legend ? { legend } : {}),
    setStyle: (style) => run('sequence.style', { style }),
    setColors: (scheme) => run('sequence.colors', { scheme: scheme ?? 'none' }),
    schemes: () => {
      const v = viewer();
      return LetterSchemes.map((id) => {
        const swatch = v?.swatchOf(id);
        return { id, label: id === null ? 'No colors' : (SchemeNames[id] ?? id), ...(swatch ? { swatch } : {}) };
      });
    },
    openAddTrack: () => source.band.value?.openAddTrack?.(),
  };
}
