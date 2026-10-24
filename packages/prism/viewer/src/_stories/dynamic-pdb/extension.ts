/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Model } from '@molstar/model/model/structure';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import type { CommandDefinition } from '@molstar/prism-ui';
import { type ObservableState, observableState } from '@molstar/prism-ui/reactive/observable-state';
import type { ReadonlyState } from '@molstar/prism-ui';
import type { DocNode, Op } from '@molstar/seq-core/format';
import { map, merge, type Observable, startWith } from 'rxjs';
import { ApiArea, ApiResult, type ViewerApi } from '../../api/area.js';
import { loadedStructures } from '../../api/scope.js';
import { streamingOf } from '../../api/volume.js';
import { operationCommand } from '../../commands/operations.js';
import type { SequenceBehavior, SequenceShown } from '../../sequence/behaviors.js';
import type { ViewerExtension } from '../../spec.js';
import {
  conformers,
  type DpdbArtifact,
  type DpdbEntry,
  type DpdbModel,
  type DpdbSource,
  labelChains,
  mutationsTable,
  parseArtifacts,
  parseEntry,
  parseModels,
  type ResidueField,
  residueValues,
  valuesTable,
} from './data.js';

declare module '../../api.js' {
  // The story's own area (api.md §2, extensions add areas).
  interface ExtensionApis {
    'dynamic-pdb': DynamicPdbApi;
  }
}

// Dynamic PDB in the viewer, as an extension a host would write (viewer-spec.md §4), kept with its story: an area that
// loads an entry's models and says what it is (`dynamic-pdb.load`, `dynamic-pdb.entry`), a card with the entry and its
// models' metrics (`ui.tsx`), and the sequence's tracks from the models' per-residue data, each shown or hidden by
// `dynamic-pdb.track`. The card streams the deposited model's density from the volume server (`volume.streaming`:
// 2Fo-Fc and Fo-Fc for X-ray, EMDB's map for EM). It reads through a `DpdbSource` (data.ts): the fixtures in the story.

const Prefix = 'dynamic-pdb';
export const DpdbTracks = {
  mutations: `${Prefix}.mutations-track`,
  bIso: `${Prefix}.b-iso-track`,
  rscc: `${Prefix}.rscc-track`,
  conformers: `${Prefix}.conformers-track`,
  spread: `${Prefix}.spread-track`,
} as const;

const DpdbTrackLabels: Readonly<Record<string, string>> = {
  [DpdbTracks.mutations]: 'Mutations',
  [DpdbTracks.bIso]: 'B-factor',
  [DpdbTracks.rscc]: 'RSCC',
  [DpdbTracks.conformers]: 'Alternate conformers',
  [DpdbTracks.spread]: 'Cα spread (Å)',
};

/** A model of the loaded entry, and the structure it is in the viewer (its state ref). */
export interface DpdbLoadedModel {
  readonly model: DpdbModel;
  readonly artifact?: DpdbArtifact;
  readonly ref?: string;
}

export interface DpdbLoaded {
  readonly entry: DpdbEntry;
  readonly models: readonly DpdbLoadedModel[];
}

/** A track the sequence shows for the entry: the card lists them, each with its eye. */
export interface DpdbTrackInfo {
  readonly id: string;
  readonly label: string;
  readonly color: string;
}

/** What the area, the sequence behavior and the card share (the extension's, one per viewer). */
interface DpdbShared {
  /** The entry the area loaded last. */
  readonly loaded: ObservableState<DpdbLoaded | undefined>;
  /** The tracks the sequence has for the document it shows. */
  readonly tracks: ObservableState<readonly DpdbTrackInfo[]>;
  /** The tracks the person hid (by id): kept across documents. */
  readonly hidden: ObservableState<ReadonlySet<string>>;
}

export class DynamicPdbApi extends ApiArea {
  readonly id = 'dynamic-pdb';

  readonly #shared: DpdbShared;

  constructor(
    ctx: ConstructorParameters<typeof ApiArea>[0],
    readonly source: DpdbSource,
    shared: DpdbShared,
  ) {
    super(ctx);
    this.#shared = shared;
  }

  get loaded(): ReadonlyState<DpdbLoaded | undefined> {
    return this.#shared.loaded.asReadonly;
  }

  get tracks(): ReadonlyState<readonly DpdbTrackInfo[]> {
    return this.#shared.tracks.asReadonly;
  }

  get hiddenTracks(): ReadonlyState<ReadonlySet<string>> {
    return this.#shared.hidden.asReadonly;
  }

  /** Each change of the state tree (a model shown or hidden, its density streamed or stopped): for the views. */
  get changed(): Observable<void> {
    const events = this.ctx.plugin.state.data.events;
    return merge(events.changed, events.cell.stateUpdated).pipe(
      map(() => undefined),
      startWith(undefined),
    );
  }

  #structure(ref: string) {
    return loadedStructures(this.ctx).find((s) => s.cell.transform.ref === ref);
  }

  /** The model of a structure of the entry (by state ref); undefined when it is gone. */
  modelOf(ref: string): Model | undefined {
    return this.#structure(ref)?.cell.obj?.data.model;
  }

  isShown(ref: string): boolean {
    const s = this.#structure(ref);
    return !!s && !s.cell.state.isHidden;
  }

  /** Whether the structure's density streams, and what from (EMDB's map, or X-ray's 2Fo-Fc and Fo-Fc). */
  density(ref: string): { readonly on: boolean; readonly kind?: 'em' | 'x-ray'; readonly id?: string } {
    const s = this.#structure(ref);
    const on = s && streamingOf(this.ctx.plugin, s);
    return on ? { on: true, kind: on.entry.kind === 'em' ? 'em' : 'x-ray', id: on.entry.dataId } : { on: false };
  }

  /** The PDB id of the entry loaded or loading: the scene's rows ask before the load has said which refs it made. */
  #pdb: string | undefined;

  #entryIdOf(ref: string) {
    return this.modelOf(ref)?.entryId.toUpperCase();
  }

  /** Whether a structure is a model of the entry (qFit's is named 3k7f_qFit, the deposited 3K7F). */
  isOfEntry(ref: string): boolean {
    return !!this.#pdb && !!this.#entryIdOf(ref)?.startsWith(this.#pdb);
  }

  /**
   * Whether the volume server has the structure's density: the server knows the PDB's ids, so only the deposited
   * model's; the models share the crystal's frame, so its map is theirs too.
   */
  hasDensity(ref: string): boolean {
    return !!this.#pdb && this.#entryIdOf(ref) === this.#pdb;
  }

  /** The loaded model a structure is (by its Mol* model). */
  loadedOf(model: Model): DpdbLoadedModel | undefined {
    return this.#shared.loaded.value?.models.find((m) => m.ref && this.modelOf(m.ref)?.id === model.id);
  }

  readonly load = this.operation({
    id: 'dynamic-pdb.load',
    title: 'Load a Dynamic PDB entry',
    params: { entry: PD.Text('dpdb_gijqn2no', { description: 'A Dynamic PDB entry id (dpdb_…)' }) },
    modules: ['network'],
    approval: 'network',
    run: async ({ entry: id }, ctx): Promise<ApiResult<DpdbLoaded>> => {
      const entryId = id.trim();
      if (!/^dpdb_\w+$/.test(entryId)) return ApiResult.fail('invalid', `"${id}" is not a Dynamic PDB entry id.`);
      let entry: DpdbEntry;
      let models: DpdbModel[];
      try {
        entry = parseEntry(await this.source.entry(entryId));
        models = parseModels(await this.source.models(entryId));
      } catch (e) {
        return ApiResult.fail('failed', `Could not read ${entryId}: ${(e as Error).message}.`);
      }
      this.#pdb = entry.pdb?.toUpperCase();
      const loaded: DpdbLoadedModel[] = [];
      for (const model of models) {
        const artifacts = parseArtifacts(await this.source.artifacts(model.id));
        const artifact =
          artifacts.find((a) => a.id === model.primaryArtifact) ?? artifacts.find((a) => a.type === 'model');
        if (!artifact) {
          loaded.push({ model });
          continue;
        }
        const before = loadedStructures(ctx).map((s) => s.cell.transform.ref);
        const file = this.source.file(artifact);
        const label = `${entry.pdb ?? entry.id} · ${model.title}`;
        const r = await ctx.api.run('load.data', { url: file.url, format: file.format, label });
        if (!r.ok) return r as ApiResult<never>;
        const ref = loadedStructures(ctx).find((s) => !before.includes(s.cell.transform.ref))?.cell.transform.ref;
        loaded.push({ model, artifact, ...(ref ? { ref } : {}) });
      }
      const value: DpdbLoaded = { entry, models: loaded };
      this.#shared.loaded.set(value);
      const names = loaded.map((m) => m.model.title).join(' and ');
      return ApiResult.ok(
        `Loaded ${entry.id} (${entry.pdb ?? entry.title}): ${names}.`,
        { nodes: loaded.flatMap((m) => (m.ref ? [m.ref] : [])) },
        { value },
      );
    },
  });

  readonly track = this.operation({
    id: 'dynamic-pdb.track',
    title: 'Show or hide a Dynamic PDB track in the sequence',
    params: {
      track: PD.Select<string>(DpdbTracks.bIso, Object.entries(DpdbTrackLabels), { description: 'The track' }),
      shown: PD.Boolean(true),
    },
    run: async ({ track, shown }): Promise<ApiResult<unknown>> => {
      const label = DpdbTrackLabels[track as keyof typeof DpdbTrackLabels] ?? track;
      const hidden = this.#shared.hidden.value;
      if (hidden.has(track) !== shown) return ApiResult.already(`${label} is already ${shown ? 'shown' : 'hidden'}.`);
      const next = new Set(hidden);
      if (shown) next.delete(track);
      else next.add(track);
      this.#shared.hidden.set(next);
      return ApiResult.ok(`${shown ? 'Showed' : 'Hid'} the ${label} track.`);
    },
  });

  readonly entry = this.operation({
    id: 'dynamic-pdb.entry',
    title: 'Read the loaded Dynamic PDB entry',
    params: {},
    read: true,
    run: async (): Promise<ApiResult<unknown>> => {
      const l = this.#shared.loaded.value;
      if (!l) return ApiResult.fail('not-found', 'No Dynamic PDB entry is loaded.', 'Load one with dynamic-pdb.load.');
      const models = l.models.map((m) => ({ title: m.model.title, type: m.model.modelType, metrics: m.model.metrics }));
      const e = l.entry;
      return ApiResult.ok(
        `${e.id}: ${e.title} (${e.method}${e.resolution ? `, ${e.resolution} Å` : ''}); ${models.length} models.`,
        {},
        { value: { id: e.id, pdb: e.pdb, title: e.title, method: e.method, resolution: e.resolution, models } },
      );
    },
  });
}

// --- the sequence's tracks ---------------------------------------------------------------------------------------------

const TrackColors: Readonly<Record<string, string>> = {
  [DpdbTracks.mutations]: '#d0574f',
  [DpdbTracks.bIso]: '#6f86b8',
  [DpdbTracks.rscc]: '#4f9a78',
  [DpdbTracks.conformers]: '#9a6fb8',
  [DpdbTracks.spread]: '#c08a3e',
};
const uniform = (color: string) => ({ theme: 'uniform', params: { color } });

interface TrackNodes {
  readonly id: string;
  readonly data: DocNode;
  readonly track: DocNode;
}

/** A track over an inline annotation in the units' `auth` numbering: its data node and the track with one layer. */
function trackNodes(
  id: string,
  geometry: 'values' | 'ranges',
  value: object,
  layer: Readonly<Record<string, unknown>>,
  hidden: boolean,
  height = 1.5,
): TrackNodes {
  const label = DpdbTrackLabels[id];
  const data = `${id}.data`;
  return {
    id,
    data: { kind: 'annotation', ref: data, params: { geometry, numbering: 'auth', label, value } },
    track: {
      kind: 'track',
      ref: id,
      params: { label, height, hidden },
      children: [{ kind: 'layer', ref: `${id}.layer`, params: { ...layer, data, color: uniform(TrackColors[id]) } }],
    },
  };
}

/** The tracks for the shown structure, when it is a model of the loaded entry; the hidden ones hidden. */
function tracksFor(area: DynamicPdbApi, shown: SequenceShown, hidden: ReadonlySet<string>): TrackNodes[] {
  const loaded = area.loaded.value;
  if (!loaded || !area.loadedOf(shown.structure.model)) return [];
  // The units are the complex's derived refs (`structure/A`: seqspec format.md, nodes and refs).
  const chains = new Map([...labelChains(shown.structure.model)].map(([k, asym]) => [k, `${shown.complex}/${asym}`]));
  // The per-residue data is the deposited model's (the one the API has it for).
  const residues = loaded.models.find((m) => m.model.residueData.length)?.model.residueData ?? [];
  const out: TrackNodes[] = [];
  const track = (
    id: string,
    geometry: 'values' | 'ranges',
    value: object,
    layer: Record<string, unknown>,
    height?: number,
  ) => out.push(trackNodes(id, geometry, value, layer, hidden.has(id), height));
  const variants = mutationsTable(
    chains,
    residues,
    loaded.entry.entities.flatMap((e) => e.mutations),
  );
  if (variants.target.length) track(DpdbTracks.mutations, 'ranges', variants, { type: 'variants' }, 1);
  const field = (f: ResidueField) => valuesTable(chains, residueValues(residues, f));
  if (residues.length) {
    track(DpdbTracks.bIso, 'values', field('bIso'), { type: 'values', encoding: 'line' });
    track(DpdbTracks.rscc, 'values', field('rscc'), { type: 'values', encoding: 'line', domain: [0, 1] });
  }
  // qFit's conformers, against the deposited model's Cα (the same crystal frame).
  const multi = loaded.models.find((m) => m.model.modelType === 'Multiconformer' && m.ref);
  const reference = loaded.models.find((m) => m.model.modelType === 'Deposited' && m.ref);
  const multiModel = multi?.ref && area.modelOf(multi.ref);
  const refModel = reference?.ref && area.modelOf(reference.ref);
  if (multiModel && refModel) {
    const c = conformers(multiModel, refModel);
    // The alternates beyond the first: a residue qFit kept single draws nothing.
    const count = valuesTable(
      chains,
      c.map((x) => ({ ...x, value: x.count - 1 })),
    );
    const spread = valuesTable(
      chains,
      c.map((x) => ({ ...x, value: x.spread })),
    );
    track(DpdbTracks.conformers, 'values', count, { type: 'values', encoding: 'bars', domain: [0, 3] });
    track(DpdbTracks.spread, 'values', spread, { type: 'values', encoding: 'line' });
  }
  return out;
}

const Legend = [
  { label: 'Mutation (from UniProt)', id: DpdbTracks.mutations },
  { label: 'B-factor (Å²)', id: DpdbTracks.bIso },
  { label: 'Real-space correlation', id: DpdbTracks.rscc },
  { label: 'qFit alternate conformers', id: DpdbTracks.conformers },
  { label: 'Cα spread from deposited (Å)', id: DpdbTracks.spread },
];

/**
 * The tracks under the viewer's (sequence.md §3): added once the document and the entry are both in, each after the one
 * before it; the person's eyes (`dynamic-pdb.track`) set their `hidden`, so a hidden track keeps its place.
 */
function sequenceBehavior(area: () => DynamicPdbApi | undefined, shared: DpdbShared): SequenceBehavior {
  return {
    id: `${Prefix}.sequence.tracks`,
    title: 'Dynamic PDB tracks',
    attach(ctx) {
      // The document the tracks are in, and their ids.
      let added: { shown: SequenceShown; ids: readonly string[] } | undefined;
      const add = () => {
        const a = area();
        const shown = ctx.shown.value;
        if (added && added.shown === shown) return;
        const tracks = a && shown ? tracksFor(a, shown, shared.hidden.value) : [];
        added = undefined;
        shared.tracks.set(tracks.map((t) => ({ id: t.id, label: DpdbTrackLabels[t.id], color: TrackColors[t.id] })));
        if (!shown || !tracks.length) {
          ctx.offer(undefined);
          return;
        }
        const ops: Op[] = [];
        let after = shown.after;
        for (const t of tracks) {
          ops.push({ op: 'insert', parent: shown.data, node: t.data });
          ops.push({ op: 'insert', parent: shown.section, after, node: t.track });
          after = t.id;
        }
        added = { shown, ids: tracks.map((t) => t.id) };
        void ctx.patch(shown, ops);
        const ids = new Set(added.ids);
        ctx.offer({
          menus: {
            [DpdbTracks.bIso]: [
              { id: 'color-b-iso', label: 'Color by B-factor', icon: 'paintbrush', theme: 'uncertainty' },
            ],
            [DpdbTracks.conformers]: [
              { id: 'color-occupancy', label: 'Color by occupancy', icon: 'paintbrush', theme: 'occupancy' },
            ],
          },
          legend: [
            {
              title: 'Dynamic PDB',
              items: Legend.filter((l) => ids.has(l.id)).map((l) => ({ label: l.label, color: TrackColors[l.id] })),
            },
          ],
        });
      };
      // The person's eyes: each track's `hidden` set where it is.
      const hide = (hidden: ReadonlySet<string>) => {
        if (!added) return;
        void ctx.patch(
          added.shown,
          added.ids.map((ref): Op => ({ op: 'set', ref, params: { hidden: hidden.has(ref) } })),
        );
      };
      const subs = [
        ctx.shown.subscribe(add),
        shared.loaded.asReadonly.subscribe(add),
        shared.hidden.asReadonly.subscribe(hide),
      ];
      return () => {
        for (const s of subs) s.unsubscribe();
        shared.tracks.set([]);
      };
    },
  };
}

/** The card's id (`ui.tsx`), here so the commands do not load the views. */
export const DpdbCard = `${Prefix}.overlay.entry`;

function commands(api: ViewerApi, area: DynamicPdbApi): CommandDefinition<any>[] {
  return [
    operationCommand(api, area.load, { category: 'Data', icon: 'download', alias: ['dynamic pdb', 'qfit', 'dpdb'] }),
    {
      id: `${Prefix}.show`,
      title: 'Show the Dynamic PDB entry',
      category: 'Data',
      icon: 'layers',
      alias: ['dynamic pdb', 'qfit'],
      async run() {
        const r = await api.view.open({ card: DpdbCard });
        if (ApiResult.isFailure(r)) throw new Error(r.message);
        return { message: 'Opened the Dynamic PDB card.', ops: `view.open({ card: "${DpdbCard}" })` };
      },
    },
  ];
}

/** The tracks hidden at first (their eyes in the card show them): the mutations and qFit's alternates tell the story. */
export const DpdbHiddenTracks: readonly string[] = [DpdbTracks.bIso, DpdbTracks.rscc, DpdbTracks.spread];

/** The extension over a source (the fixtures in the story); `hidden`: the tracks hidden at first. */
export function dynamicPdbExtension(
  source: DpdbSource,
  options: { readonly hidden?: readonly string[] } = {},
): ViewerExtension {
  let area: DynamicPdbApi | undefined;
  const shared: DpdbShared = {
    loaded: observableState<DpdbLoaded | undefined>(undefined),
    tracks: observableState<readonly DpdbTrackInfo[]>([]),
    hidden: observableState<ReadonlySet<string>>(new Set(options.hidden ?? DpdbHiddenTracks)),
  };
  return {
    id: Prefix,
    title: 'Dynamic PDB',
    api: { [Prefix]: (ctx) => (area = new DynamicPdbApi(ctx, source, shared)) },
    sequence: { behaviors: [sequenceBehavior(() => area, shared)] },
    requires: ['structure'],
    commands: (api) => (api[Prefix] ? commands(api, api[Prefix]) : []),
    ui: async (ctx) => (await import('./ui.js')).dynamicPdbUI(ctx),
  };
}
