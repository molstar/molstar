/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { mapState, type ReadonlyState } from '@molstar/prism-ui/core';
import { ParamsFormModel } from '@molstar/prism-ui/forms/params-form';
import type { SnapshotStepperModel } from '@molstar/prism-ui/model/playback';
import type { Subscription } from 'rxjs';
import type { ViewerApi } from '../../api/area.js';
import type { SnapshotEntry, SnapshotsModel } from '../../model/snapshots.js';
import type { SnapshotApi } from './operations.js';
import type { SnapshotInfo, SnapshotsState } from './source.js';

// The view models of the Snapshots tab and the playback card, over the snapshot area: they read its state (the source
// in this folder) and act through the person's API, never the plugin. DOM-free.

const Transition = PD.MappedStatic(
  'animate',
  {
    animate: PD.Group(
      {
        durationInMs: PD.Numeric(250, { min: 100, max: 5000, step: 50 }, { label: 'Duration (ms)' }),
        easing: PD.Select('linear', [
          ['linear', 'Linear'],
          ['ease-in-out', 'Ease in-out'],
          ['ease-out', 'Ease out'],
        ]),
      },
      { isFlat: true },
    ),
    instant: PD.Group({}, { isFlat: true }),
  },
  {
    label: 'Camera',
    options: [
      ['animate', 'Animate'],
      ['instant', 'Instant'],
    ],
  },
);

/** What a new snapshot includes: Mol*'s snapshot options (`PluginState.SnapshotParams`), as the tab shows them. */
export const SaveOptions = {
  durationInMs: PD.Numeric(1500, { min: 100, max: 15000, step: 100 }, { label: 'Show for (ms)' }),
  data: PD.Boolean(true, { label: 'Data' }),
  camera: PD.Boolean(true, { label: 'Camera' }),
  cameraTransition: Transition,
  canvas3d: PD.Boolean(true, { label: 'Canvas settings' }),
  componentManager: PD.Boolean(true, { label: 'Component options' }),
  structureSelection: PD.Boolean(false, { label: 'Selection' }),
  animation: PD.Boolean(true, { label: 'Animation' }),
  startAnimation: PD.Boolean(false, { label: 'Start animation' }),
  behavior: PD.Boolean(false, { label: 'Behaviors' }),
  interactivity: PD.Boolean(true, { label: 'Interactivity' }),
  image: PD.Boolean(false, { label: 'Thumbnail', description: 'Save a small image with each snapshot.' }),
};

/** A snapshot's own playback settings. */
export const EntryOptions = {
  durationInMs: PD.Numeric(1500, { min: 100, max: 15000, step: 100 }, { label: 'Show for (ms)' }),
};

const DefaultShowFor = EntryOptions.durationInMs.defaultValue;

/** The name the list and the card show: an unnamed snapshot is its position. */
const titleOf = (e: SnapshotInfo, index: number) => e.name || `Snapshot ${index + 1}`;

const clock = (ms: number) => new Date(ms).toLocaleTimeString([], { hour: 'numeric', minute: '2-digit' });

const toEntry = (e: SnapshotInfo, index: number): SnapshotEntry => ({
  id: e.id,
  name: titleOf(e, index),
  ...(e.key ? { key: e.key } : {}),
  ...(e.description ? { description: e.description } : {}),
  timestamp: clock(e.timestamp),
  ...(e.image ? { image: e.image } : {}),
});

const sameEntries = (a: readonly SnapshotEntry[], b: readonly SnapshotEntry[]) =>
  a.length === b.length &&
  a.every((x, i) => {
    const y = b[i];
    return (
      x.id === y.id &&
      x.name === y.name &&
      x.key === y.key &&
      x.description === y.description &&
      x.timestamp === y.timestamp &&
      x.image === y.image
    );
  });

export interface SnapshotsPanelModel extends SnapshotsModel {
  dispose(): void;
}

/** What `session.save` gives (the export-session feature): the file's bytes, saved by the page. */
type Saved = { blob: Blob; filename: string };

export function createSnapshotsModel(
  api: ViewerApi,
  area: SnapshotApi,
  /** Hands a saved file to the page. */
  save: (blob: Blob, filename: string) => void,
): SnapshotsPanelModel {
  const state = area.state;
  const entries = mapState(state, (s) => s.entries.map(toEntry), sameEntries);
  const current = mapState(state, (s) => s.current, Object.is);
  const forms = new Map<string, ParamsFormModel<typeof EntryOptions>>();
  const timers = new Set<ReturnType<typeof setTimeout>>();
  const subs: Subscription[] = [];
  const entryOf = (id: string) => state.value.entries.find((e) => e.id === id);
  const saveOptions = new ParamsFormModel(SaveOptions);

  return {
    entries,
    current,
    saveOptions,
    entryOptions(id) {
      if (!entryOf(id)) return undefined;
      let form = forms.get(id);
      if (!form) {
        let timer: ReturnType<typeof setTimeout> | undefined;
        form = new ParamsFormModel(
          EntryOptions,
          { durationInMs: entryOf(id)?.showFor ?? DefaultShowFor },
          {
            onChange: (values) => {
              // A slider drags through many values: the edit is made when it rests.
              if (timer) {
                clearTimeout(timer);
                timers.delete(timer);
              }
              timer = setTimeout(() => {
                if (timer) timers.delete(timer);
                void api.snapshot?.update({ entry: id, showFor: values.durationInMs });
              }, 150);
              timers.add(timer);
            },
          },
        );
        forms.set(id, form);
        const shown = form;
        // Follows the list (a session opened, another edit): loaded without reporting an edit.
        subs.push(
          state.subscribe((s) => {
            const e = s.entries.find((x) => x.id === id);
            if (e) shown.load({ durationInMs: e.showFor ?? DefaultShowFor });
          }),
        );
      }
      return form;
    },
    add: (name) => void api.snapshot?.add({ ...(name ? { name } : {}), options: { ...saveOptions.values.value } }),
    apply: (id) => void api.snapshot?.apply({ entry: id }),
    update: (id, patch) => {
      const e = entryOf(id);
      if (!e) return;
      // The operation keeps what is empty, unless `clear` says to empty it: the whole edit goes in.
      void api.snapshot?.update({
        entry: id,
        name: patch.name ?? e.name,
        description: patch.description ?? e.description ?? '',
        key: patch.key ?? e.key ?? '',
        clear: true,
      });
    },
    replace: (id) => void api.snapshot?.replace({ entry: id, options: { ...saveOptions.values.value } }),
    move: (id, by) => {
      const index = state.value.entries.findIndex((e) => e.id === id);
      const to = index + by;
      if (index < 0 || to < 0 || to >= state.value.entries.length) return;
      void api.snapshot?.move({ entry: id, position: to + 1 });
    },
    moveTo: (id, index) => {
      const count = state.value.entries.length;
      void api.snapshot?.move({ entry: id, position: Math.max(0, Math.min(count - 1, index)) + 1 });
    },
    remove: (id) => void api.snapshot?.remove({ entry: id }),
    clear: () => void api.snapshot?.clear(),
    download(kind) {
      // The export-session feature brings `session.save`; without it there is nothing to save.
      if (!api.describe('session.save')) return;
      void api.run('session.save', { kind }).then((r) => {
        const saved = r.ok ? (r.value as Saved | undefined) : undefined;
        if (saved) save(saved.blob, saved.filename);
      });
    },
    open(files) {
      if (files[0]) void api.session.open({ file: files[0] });
    },
    dispose() {
      for (const t of timers) clearTimeout(t);
      timers.clear();
      for (const s of subs.splice(0)) s.unsubscribe();
    },
  };
}

/** The playback card's model: the current snapshot, its place in the list, and stepping and playing. */
export function createSnapshotStepper(api: ViewerApi, area: SnapshotApi): SnapshotStepperModel {
  const state: ReadonlyState<SnapshotsState> = area.state;
  const indexOf = (s: SnapshotsState) => s.entries.findIndex((e) => e.id === s.current);
  const sameTitles = (a: readonly string[], b: readonly string[]) =>
    a.length === b.length && a.every((t, i) => t === b[i]);
  return {
    count: mapState(state, (s) => s.entries.length, Object.is),
    index: mapState(state, indexOf, Object.is),
    title: mapState(
      state,
      (s) => {
        const i = indexOf(s);
        return i < 0 ? undefined : titleOf(s.entries[i], i);
      },
      Object.is,
    ),
    titles: mapState(state, (s) => s.entries.map(titleOf), sameTitles),
    description: mapState(state, (s) => s.entries[indexOf(s)]?.description, Object.is),
    isPlaying: mapState(state, (s) => s.playing, Object.is),
    previous: () => void api.snapshot?.previous(),
    next: () => void api.snapshot?.next(),
    goTo(index) {
      const e = state.value.entries[index];
      if (e) void api.snapshot?.apply({ entry: e.id });
    },
    togglePlay: () => void (state.value.playing ? api.snapshot?.stop() : api.snapshot?.play()),
    openLink(href) {
      // `#key` applies the snapshot with that key; a web link opens as one.
      if (href.startsWith('#')) {
        const e = state.value.entries.find((x) => x.key === href.slice(1));
        if (e) void api.snapshot?.apply({ entry: e.id });
      } else if (/^https?:\/\//i.test(href)) {
        globalThis.open?.(href, '_blank', 'noopener,noreferrer');
      }
    },
  };
}
