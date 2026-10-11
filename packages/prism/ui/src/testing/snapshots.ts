/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { ParamsFormModel } from '../forms/params-form.js';
import type { SnapshotEntry, SnapshotsModel } from '../model/snapshots.js';
import { observableState } from '../reactive/observable-state.js';

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

// Mol*'s PluginState.SnapshotParams: what a new snapshot includes.
const SaveOptions = {
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

// A snapshot's own playback settings.
const EntryOptions = {
  durationInMs: PD.Numeric(1500, { min: 100, max: 15000, step: 100 }, { label: 'Show for (ms)' }),
  cameraTransition: Transition,
};

const Thumbs = [
  'radial-gradient(circle at 40% 55%, #8f6fd6 0 22%, #cfd4e6 23%)',
  'radial-gradient(circle at 55% 45%, #d06a8c 0 18%, #6c7ad6 19% 30%, #d9dbe6 31%)',
  'radial-gradient(circle at 50% 50%, #5ea3d6 0 30%, #e1e3ea 31%)',
];

/** Snapshots for the stories: the 1TQN tour (with thumbnails when `images`), or none. */
export function mockSnapshotsModel(
  post: (text: string, kind?: 'info' | 'success') => void,
  options: { empty?: boolean; images?: boolean; count?: number } = {},
): SnapshotsModel {
  let next = 1;
  const make = (name: string, rest: Partial<SnapshotEntry> = {}): SnapshotEntry => {
    const id = `s${next}`;
    const image = options.images ? Thumbs[(next - 1) % Thumbs.length] : undefined;
    next++;
    return { id, name, image, timestamp: 'today', ...rest };
  };
  const tour: SnapshotEntry[] = options.empty
    ? []
    : options.count
      ? Array.from({ length: options.count }, (_, i) =>
          make(['Overview', 'Active site', 'Heme pocket', 'Substrate channel', 'Surface'][i % 5] + ` ${i + 1}`, {
            description: i % 2 ? `Step ${i + 1} of the tour.` : undefined,
          }),
        )
      : [
          make('Overview', {
            key: 'overview',
            description: 'Cytochrome **P450 3A4** with its heme cofactor.\n\nSee the [active site](#active-site).',
            timestamp: '9:05',
          }),
          make('Active site', {
            key: 'active-site',
            description: '### Active site\n- Heme iron coordinated by `Cys442`',
            timestamp: '9:07',
          }),
          make('Surface', { timestamp: '9:12' }),
        ];
  const entries = observableState<readonly SnapshotEntry[]>(tour);
  const current = observableState<string | undefined>(tour[1]?.id);
  const forms = new Map<string, ParamsFormModel<any>>();
  const set = (id: string, patch: Partial<SnapshotEntry>) =>
    entries.set(entries.value.map((e) => (e.id === id ? { ...e, ...patch } : e)));
  const name = (id: string) => entries.value.find((e) => e.id === id)?.name;
  return {
    entries: entries.asReadonly,
    current: current.asReadonly,
    saveOptions: new ParamsFormModel(SaveOptions),
    entryOptions: (id) => {
      if (!entries.value.some((e) => e.id === id)) return undefined;
      let form = forms.get(id);
      if (!form) forms.set(id, (form = new ParamsFormModel(EntryOptions)));
      return form;
    },
    add: (n) => {
      const entry = make(n ?? `Snapshot ${entries.value.length + 1}`, { timestamp: 'just now' });
      entries.set([...entries.value, entry]);
      current.set(entry.id);
      post(`Saved snapshot “${entry.name}”`, 'success');
    },
    apply: (id) => {
      current.set(id);
      post(`Applied “${name(id)}”`);
    },
    update: (id, patch) => set(id, patch),
    replace: (id) => {
      set(id, { timestamp: 'just now' });
      post(`Replaced “${name(id)}” with the current view`, 'success');
    },
    move: (id, by) => {
      const list = [...entries.value];
      const i = list.findIndex((e) => e.id === id);
      const j = i + by;
      if (i < 0 || j < 0 || j >= list.length) return;
      [list[i], list[j]] = [list[j], list[i]];
      entries.set(list);
    },
    moveTo: (id, index) => {
      const list = entries.value.filter((e) => e.id !== id);
      const entry = entries.value.find((e) => e.id === id);
      if (!entry) return;
      list.splice(Math.max(0, Math.min(list.length, index)), 0, entry);
      entries.set(list);
    },
    remove: (id) => {
      post(`Removed “${name(id)}”`);
      entries.set(entries.value.filter((e) => e.id !== id));
      if (current.value === id) current.set(undefined);
    },
    clear: () => {
      entries.set([]);
      current.set(undefined);
    },
    download: (kind) => post(kind === 'session' ? 'Saved state.molx' : 'Saved state.molj', 'success'),
    open: (files) => post(`Opening ${files.map((f) => f.name).join(', ')}…`),
  };
}
