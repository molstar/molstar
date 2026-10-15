/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { StateObjectSelector } from '@molstar/core/state';
import { Asset } from '@molstar/core/util/assets';
import { Color } from '@molstar/core/util/color';
import { getFileNameInfo } from '@molstar/core/util/file-info';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { StringLike } from '@molstar/core/util/string-like';
import { Volume } from '@molstar/model/model/volume';
import { loadMVSData } from '@molstar/mvs/components/formats';
import { PluginConfig } from '@molstar/plugin/config';
import { DownloadStructure, PdbDownloadProvider } from '@molstar/plugin/state/actions/structure';
import { DownloadDensity } from '@molstar/plugin/state/actions/volume';
import { TopologyFormatCategory } from '@molstar/plugin/state/formats/topology/category';
import { createVolumeRepresentationParams } from '@molstar/plugin/state/helpers/volume-representation-params';
import type { PluginStateObject } from '@molstar/plugin/state/objects';
import { Download, LazyVolume } from '@molstar/plugin/state/transforms/data/fetch';
import { TrajectoryFromModelAndCoordinates } from '@molstar/plugin/state/transforms/structure/hierarchy';
import { VolumeRepresentation3D } from '@molstar/plugin/state/transforms/volume/representation';
import { type ApiContext, ApiResult, memberName, type OperationDefinition } from '../../api/area.js';
import { messageOf } from '../../api/run.js';
import { labelOf, loadedStructures } from '../../api/scope.js';
import { currentStateServer, openStateFrom, stateAddress, transformersOf } from '../../api/state-open.js';

// The `loaders` feature's operations (operations.md §2, the `load` row; README §3.2): Mol*'s own loaders
// (`@molstar/plugin-extension/loaders`: the plugin's state actions and builders) as operations of the core's `load`
// area. That package declares `@molstar/plugin-ui` (React), which the viewer must not depend on (README §3.2 table), so
// the few calls are ported here over `@molstar/plugin` directly, with the same state they make. A load that fails
// returns `failed` and the run's checkpoint restores: nothing half-built. Loading something that is loaded says so
// instead (decision 30, operations.md §1: identity is the entry id, or the URL).

const network = { modules: ['network'] as const, approval: 'network' as const };

/** What a load did, as `value`: the names of what it made. */
export interface Loaded {
  readonly loaded: readonly string[];
}

// --- load.entry ---------------------------------------------------------------------------------------------------

export const EntrySources = ['pdb', 'pdb-ihm', 'pdb-dev', 'afdb', 'emdb', 'model-archive'] as const;
export type EntrySource = (typeof EntrySources)[number];

const SourceOptions: [EntrySource, string][] = [
  ['pdb', 'PDB'],
  ['pdb-ihm', 'PDB-IHM'],
  ['pdb-dev', 'PDB-Dev (deprecated, the same as PDB-IHM)'],
  ['afdb', 'AlphaFold DB'],
  ['emdb', 'EMDB'],
  ['model-archive', 'Model Archive'],
];

export const EntryParams = {
  source: PD.Select<EntrySource>('pdb', SourceOptions, {
    description: 'The database: pdb, pdb-ihm (pdb-dev is accepted), afdb, emdb or model-archive',
  }),
  id: PD.Text('', {
    description:
      'The entry id: 1TQN, 8ZZC (pdb-ihm), Q8W3K0 (afdb, UniProt accession), EMD-1832 (emdb), ma-bak-cepc-0003; several for structures, separated by commas',
  }),
  assembly: PD.Text('', { description: 'Structures: the assembly id to build; empty: the default' }),
  representation: PD.Text('', {
    description: "Structures: a representation preset by id or name; empty: the viewer's default (AFDB: pLDDT)",
  }),
};

const PluginRepresentationDefault = 'preset-structure-representation-ma-quality-assessment-plddt';

/** The ids of a field that takes several (comma or space separated), as the download action reads them. */
export function idsOf(text: string, source: EntrySource): string[] {
  if (source === 'emdb') return text.trim() ? [text.trim()] : [];
  return text
    .split(/[,\s]/)
    .map((id) => id.trim())
    .filter((id) => !!id && (id.length >= 4 || /^[1-9][0-9]*$/.test(id)));
}

const kindOf = (source: EntrySource): Exclude<EntrySource, 'pdb-dev'> => (source === 'pdb-dev' ? 'pdb-ihm' : source);

const displayId = (source: EntrySource, id: string) =>
  ['pdb', 'pdb-ihm', 'pdb-dev', 'emdb'].includes(source) ? id.toUpperCase() : id;

/**
 * The labels the download of an entry carries, by the state actions (`state/actions/structure.ts`, `volume.ts`): the
 * data cell is the one place that still says what a structure was downloaded as, whatever server it came from.
 */
function entryLabels(source: EntrySource, id: string): string[] {
  const upper = id.toUpperCase();
  switch (kindOf(source)) {
    case 'pdb':
      return [
        `PDBe: ${id} (updated cif)`,
        `PDBe: ${id} (cif)`,
        `RCSB PDB: ${id} (bcif)`,
        `RCSB PDB: ${id} (cif)`,
        `PDBj: ${id} (cif)`,
      ].map((l) => l.toLowerCase());
    case 'pdb-ihm': {
      const ihm = (n: string) => (n.startsWith('PDBDEV_') ? n : `PDBDEV_${n.padStart(8, '0')}`);
      return [(upper.match(/^[1-9][A-Z0-9]{3}$/) ? `PDB-IHM: ${upper}` : ihm(upper)).toLowerCase()];
    }
    case 'afdb':
      return [`AlphaFold DB: ${id}`.toLowerCase()];
    case 'model-archive':
      return [`Model Archive: ${id}`.toLowerCase()];
    case 'emdb':
      return [`PDBe EMD Density Server: ${id}`, `RCSB PDB EMD Density Server: ${id}`].map((l) => l.toLowerCase());
  }
}

/** The downloads whose results are still in the scene: their label (lowercase) and URL. */
function liveDownloads(ctx: ApiContext): { label: string; url: string }[] {
  const { cells, tree } = ctx.plugin.state.data;
  const out: { label: string; url: string }[] = [];
  cells.forEach((cell) => {
    if (cell.transform.transformer !== Download) return;
    // A removed structure takes its data cell with it; one that stayed without results is not "loaded".
    if (!tree.children.get(cell.transform.ref).size) return;
    const params = cell.transform.params as { url?: string | Asset.Url; label?: string } | undefined;
    out.push({
      label: (params?.label ?? '').toLowerCase(),
      url: params?.url ? Asset.getUrl(params.url) : '',
    });
  });
  return out;
}

const refsOf = (ctx: ApiContext) => loadedStructures(ctx).map((s) => s.cell.transform.ref);

/** The names of the structures that were not there before. */
function addedNames(ctx: ApiContext, before: readonly string[]): string[] {
  const added = loadedStructures(ctx).filter((s) => !before.includes(s.cell.transform.ref));
  return [...new Set(added.map((s) => labelOf(ctx, s)))];
}

/**
 * Runs plugin work and throws what went wrong. The plugin reports a failed update and carries on (an action rolls back
 * and logs; a transform ends in an error cell), so work that returns is not work that succeeded.
 */
async function guarded<T>(ctx: ApiContext, work: () => Promise<T>): Promise<T> {
  const plugin = ctx.plugin;
  const errorRefs = () => {
    const refs = new Set<string>();
    plugin.state.data.cells.forEach((c) => c.status === 'error' && refs.add(c.transform.ref));
    return refs;
  };
  const known = errorRefs();
  const logged: string[] = [];
  const sub = plugin.events.log.subscribe((e) => {
    if (e.type === 'error') logged.push(String(e.message));
  });
  try {
    const out = await work();
    let reason: string | undefined = logged[0];
    plugin.state.data.cells.forEach((c) => {
      if (c.status === 'error' && !known.has(c.transform.ref)) reason ??= c.errorText;
    });
    if (reason) throw new Error(reason);
    return out;
  } finally {
    sub.unsubscribe();
  }
}

function pluginFailure(what: string, e: unknown, hint = 'Check the address and the network connection.') {
  return ApiResult.fail('failed', `Could not load ${what}: ${messageOf(e)}`, hint);
}

const loadedLine = (names: readonly string[], fallback: string) =>
  names.length ? `Loaded ${names.join(', ')}.` : `Loaded ${fallback}.`;

async function loadEntry(
  ctx: ApiContext,
  source: EntrySource,
  ids: readonly string[],
  o: { assembly: string; representation: string },
) {
  const plugin = ctx.plugin;
  await plugin.initialized;
  const root = plugin.state.data.root.obj!;
  if (kindOf(source) === 'emdb') {
    const provider = plugin.config.get(PluginConfig.Download.DefaultEmdbProvider)!;
    await plugin.runTask(
      plugin.state.data.applyAction(DownloadDensity, {
        source: { name: 'pdb-emd-ds' as const, params: { provider: { id: ids[0], server: provider }, detail: 3 } },
      }),
    );
    return;
  }
  const defaults = DownloadStructure.createDefaultParams(root, plugin);
  const options = { ...defaults.source.params.options } as Record<string, any>;
  if (o.assembly) {
    options.type = { name: 'assembly', params: { id: o.assembly, dynamicBonds: false, mergeBySymmetry: false } };
  }
  if (o.representation) options.representation = o.representation;
  const idText = ids.join(' ');
  let src: any;
  switch (kindOf(source)) {
    case 'pdb': {
      const provider = plugin.config.get(PluginConfig.Download.DefaultPdbProvider)!;
      src = {
        name: 'pdb',
        params: {
          provider: {
            id: idText,
            server: { name: provider, params: PdbDownloadProvider[provider].defaultValue },
          },
          options,
        },
      };
      break;
    }
    case 'pdb-ihm':
      src = { name: 'pdb-ihm', params: { provider: { id: idText, encoding: 'bcif' }, options } };
      break;
    case 'afdb':
      // The pLDDT preset is the model-archive extension's: without it the viewer's default stands.
      if (
        !o.representation &&
        plugin.builders.structure.representation.getPresets().some((p) => p.id === PluginRepresentationDefault)
      ) {
        options.representation = PluginRepresentationDefault;
      }
      src = { name: 'alphafolddb', params: { provider: { id: idText, encoding: 'bcif' }, options } };
      break;
    case 'model-archive':
      src = { name: 'modelarchive', params: { id: idText, options } };
      break;
  }
  await plugin.runTask(plugin.state.data.applyAction(DownloadStructure, { source: src }));
}

const entry: OperationDefinition<typeof EntryParams, Loaded> = {
  id: 'load.entry',
  title: 'Load an entry',
  params: EntryParams,
  ...network,
  exists({ source, id }, ctx) {
    const ids = idsOf(id, source);
    if (!ids.length) return undefined;
    const have = liveDownloads(ctx);
    const loaded = ids.every((i) => entryLabels(source, i).some((l) => have.some((d) => d.label === l)));
    return loaded
      ? `${ids.map((i) => displayId(source, i)).join(', ')} ${ids.length === 1 ? 'is' : 'are'} already loaded.`
      : undefined;
  },
  async run({ source, id, assembly, representation }, ctx) {
    let ids = idsOf(id, source);
    if (!ids.length) {
      return ApiResult.fail(
        'invalid',
        `'${id}' is not an id.`,
        source === 'emdb' ? 'Give an EMDB id such as EMD-1832.' : 'Give a database id such as 1TQN.',
      );
    }
    if (representation) {
      const presets = ctx.plugin.builders.structure.representation.getPresets();
      const found = presets.find(
        (p) => p.id === representation || p.display.name.toLowerCase() === representation.toLowerCase(),
      );
      if (!found) {
        return ApiResult.fail(
          'not-found',
          `There is no representation preset '${representation}'.`,
          `Presets: ${presets.map((p) => p.display.name).join(', ')}.`,
        );
      }
      representation = found.id;
    }
    // Those already loaded are left alone (a partial list loads the rest); `copy` loads them all again.
    if (!ctx.copy) {
      const have = liveDownloads(ctx);
      ids = ids.filter((i) => !entryLabels(source, i).some((l) => have.some((d) => d.label === l)));
    }
    const what = ids.map((i) => displayId(source, i)).join(', ');
    const before = refsOf(ctx);
    const cellsBefore = ctx.plugin.state.data.cells.size;
    try {
      await guarded(ctx, () => loadEntry(ctx, source, ids, { assembly, representation }));
    } catch (e) {
      return pluginFailure(what, e, 'Check the id and the network connection.');
    }
    if (ctx.plugin.state.data.cells.size === cellsBefore) {
      return ApiResult.fail('failed', `Loading ${what} made nothing.`, 'Check the id; the log says more.');
    }
    const names = kindOf(source) === 'emdb' ? [what] : addedNames(ctx, before);
    return ApiResult.ok(
      loadedLine(names, what),
      { nodes: refsOf(ctx).filter((r) => !before.includes(r)) },
      { value: { loaded: names } },
    );
  },
};

// --- load.url -----------------------------------------------------------------------------------------------------

export const UrlParams = {
  url: PD.Text('', { description: 'The address of a structure, volume or other data file' }),
  format: PD.Text('mmcif', {
    description: "The data format by name (mmcif, pdb, gro, ccp4, …), or 'auto' to read it from the address",
  }),
  binary: PD.Boolean(false, { description: 'The file is binary (bcif, ccp4, dcd, …)' }),
};

async function loadUrl(ctx: ApiContext, url: string, format: string, isBinary: boolean) {
  const plugin = ctx.plugin;
  await plugin.initialized;
  await plugin.dataTransaction(
    async () => {
      // 'auto' reads the format from the address, and whether it is binary from the formats that are.
      const info = getFileNameInfo(url);
      const data = await plugin.builders.data.download(
        { url, isBinary: isBinary || (format === 'auto' && plugin.dataFormats.binaryExtensions.has(info.ext)) },
        { state: { isGhost: true } },
      );
      const provider =
        format === 'auto' ? plugin.dataFormats.auto(info, data.cell!.obj!) : plugin.dataFormats.get(format);
      if (!provider) throw new Error(`Could not tell the format of ${url}`);
      const parsed = await provider.parse(plugin, data);
      await provider.visuals?.(plugin, parsed);
    },
    { rethrowErrors: true },
  );
}

const url: OperationDefinition<typeof UrlParams, Loaded> = {
  id: 'load.url',
  title: 'Load from a URL',
  params: UrlParams,
  ...network,
  exists({ url }, ctx) {
    const u = url.trim();
    return u && liveDownloads(ctx).some((d) => d.url === u) ? `${u} is already loaded.` : undefined;
  },
  async run({ url, format, binary }, ctx) {
    const address = url.trim();
    if (!address) return ApiResult.fail('invalid', 'Which address?', 'Give the URL of a data file.');
    format = format.trim() || 'mmcif';
    if (format !== 'auto') {
      // A lazy extension that reads this format loads first.
      await ctx.extensions.ensureFor({ formats: [format] });
      if (!ctx.plugin.dataFormats.has(format)) {
        return ApiResult.fail(
          'invalid',
          `Unknown data format '${format}'.`,
          `Formats: auto, ${ctx.plugin.dataFormats.types.map(([n]) => n).join(', ')}.`,
        );
      }
    }
    const before = refsOf(ctx);
    try {
      await guarded(ctx, () => loadUrl(ctx, address, format, binary));
    } catch (e) {
      return pluginFailure(address, e);
    }
    const names = addedNames(ctx, before);
    return ApiResult.ok(
      loadedLine(names, address),
      { nodes: refsOf(ctx).filter((r) => !before.includes(r)) },
      { value: { loaded: names.length ? names : [address] } },
    );
  },
};

// --- load.mvs -----------------------------------------------------------------------------------------------------

export const MvsParams = {
  url: PD.Text('', { description: 'The address of a MolViewSpec file (.mvsj, .mvsx); empty: use data' }),
  data: PD.Text('', {
    description: "The MolViewSpec itself: JSON for mvsj, 'base64,…' for mvsx; used when there is no url",
  }),
  format: PD.Select<'mvsj' | 'mvsx'>('mvsj', PD.arrayToOptions(['mvsj', 'mvsx'] as const)),
  keepCamera: PD.Boolean(false, { description: "Keep the camera as it is, ignoring the MolViewSpec's" }),
};

const mvs: OperationDefinition<typeof MvsParams, Loaded> = {
  id: 'load.mvs',
  title: 'Load a MolViewSpec',
  params: MvsParams,
  // No `modules: ['network']`: the module is checked per call, as only a url needs it (data is local).
  // It replaces the scene (and an MVS the Snapshots list): destructive wins over network (api/area.ts).
  approval: 'destructive',
  async run({ url, data, format, keepCamera }, ctx) {
    const address = url.trim();
    if (address && data.trim()) {
      return ApiResult.fail('invalid', 'Give a url or data, not both.', 'Neither is loaded: say which one is meant.');
    }
    if (!address && !data.trim()) return ApiResult.fail('invalid', 'Which MolViewSpec?', 'Give a url or its data.');
    if (address && !ctx.modules.has('network')) {
      return ApiResult.fail(
        'unavailable',
        'Load a MolViewSpec needs the network, which this viewer does not have.',
        'Give its data instead of a url.',
      );
    }
    const plugin = ctx.plugin;
    const before = refsOf(ctx);
    try {
      await plugin.initialized;
      const options = { sanityChecks: true, keepCamera };
      await guarded(ctx, async () => {
        if (!address) {
          await loadMVSData(plugin, data.trim(), format, options);
        } else if (format === 'mvsj') {
          const text = await plugin.runTask(plugin.fetch({ url: address, type: 'string' }));
          await loadMVSData(plugin, StringLike.toString(text), 'mvsj', { ...options, sourceUrl: address });
        } else {
          const bytes = await plugin.runTask(plugin.fetch({ url: address, type: 'binary' }));
          await loadMVSData(plugin, bytes, 'mvsx', { ...options, sourceUrl: address });
        }
      });
    } catch (e) {
      return pluginFailure(address || 'the MolViewSpec', e, 'Check the MolViewSpec and the network connection.');
    }
    const names = addedNames(ctx, before);
    return ApiResult.ok(
      loadedLine(names, 'the MolViewSpec'),
      { nodes: refsOf(ctx).filter((r) => !before.includes(r)), view: true, sequence: true },
      { value: { loaded: names } },
    );
  },
};

// --- load.snapshot ------------------------------------------------------------------------------------------------

export const SnapshotParams = {
  url: PD.Text('', { description: 'The address of a session (.molx) or a state (.molj); empty: use id' }),
  type: PD.Select<'auto' | 'molj' | 'molx'>(
    'auto',
    [
      ['auto', "From the address's extension"],
      ['molj', 'State (.molj)'],
      ['molx', 'Session (.molx)'],
    ],
    { description: "The file's kind; by default the address's extension says (.molx, else .molj)" },
  ),
  id: PD.Text('', { description: "A state on the state server (the viewer's `state-server` option) by id" }),
};

/** The kind a url's extension says (`.molx`, `.molj`; a query or fragment is ignored). */
export function snapshotTypeOf(url: string): 'molj' | 'molx' | undefined {
  const path = url.split(/[?#]/)[0];
  const m = /\.(molx|molj)$/i.exec(path);
  return m ? (m[1].toLowerCase() as 'molj' | 'molx') : undefined;
}

const snapshot: OperationDefinition<typeof SnapshotParams, Loaded> = {
  id: 'load.snapshot',
  title: 'Load a session',
  params: SnapshotParams,
  ...network,
  // It replaces the scene (and an MVS the Snapshots list): destructive wins over network (api/area.ts).
  approval: 'destructive',
  async run({ url, type, id }, ctx) {
    const address = url.trim();
    if (address && id.trim()) return ApiResult.fail('invalid', 'Give a url or an id, not both.');
    if (!address && !id.trim()) return ApiResult.fail('invalid', 'Which session?', 'Give a url or an id.');
    const plugin = ctx.plugin;
    try {
      await plugin.initialized;
      await guarded(ctx, async () => {
        if (address) {
          const kind = type === 'auto' ? (snapshotTypeOf(address) ?? 'molj') : type;
          const data = await plugin.runTask(plugin.fetch({ url: address, type: 'binary' }));
          // A lazy extension whose transformers the session uses loads first, as `session.open` does for a file.
          if (kind === 'molj') {
            await ctx.extensions.ensureFor({ transformers: transformersOf(new TextDecoder().decode(data)) });
          }
          await plugin.managers.snapshot.open(new File([data], `state.${kind}`));
        } else {
          const server = currentStateServer(ctx);
          if (!server) throw new Error('the viewer has no state server');
          const key = id.trim();
          await openStateFrom(ctx, stateAddress(server, key), key);
        }
      });
    } catch (e) {
      return pluginFailure(address || `state ${id}`, e, 'Check the address, the id and the network connection.');
    }
    const names = loadedStructures(ctx).map((s) => labelOf(ctx, s));
    return ApiResult.ok(
      address ? 'Loaded the session.' : `Loaded state ${id.trim()}.`,
      { view: true, sequence: true },
      { value: { loaded: names } },
    );
  },
};

// --- load.trajectory ----------------------------------------------------------------------------------------------

export const TrajectoryParams = {
  model: PD.Text('', { description: 'The address of the model (or topology) file' }),
  modelFormat: PD.Text('mmcif', { description: 'Its format: mmcif, pdb, gro, …, or a topology (psf, prmtop, top)' }),
  modelBinary: PD.Boolean(false),
  coordinates: PD.Text('', { description: 'The address of the coordinates file' }),
  coordinatesFormat: PD.Text('xtc', { description: 'Its format: xtc, dcd, trr, nctraj, lammpstrj' }),
  coordinatesBinary: PD.Boolean(true),
  preset: PD.Text('default', { description: "The structure preset: 'default' or 'all-models'" }),
};

async function loadTrajectory(
  ctx: ApiContext,
  p: {
    model: string;
    modelFormat: string;
    modelBinary: boolean;
    coordinates: string;
    coordinatesFormat: string;
    coordinatesBinary: boolean;
    preset: string;
  },
) {
  const plugin = ctx.plugin;
  await plugin.initialized;
  await plugin.dataTransaction(
    async () => {
      const modelData = await plugin.builders.data.download({ url: p.model, isBinary: p.modelBinary });
      const isTopology = plugin.dataFormats.options.some(
        ([n, , c]) => n === p.modelFormat && c === TopologyFormatCategory,
      );
      let model: StateObjectSelector;
      if (isTopology) {
        model = (await plugin.dataFormats.get(p.modelFormat)!.parse(plugin, modelData)).topology;
      } else {
        const trajectory = await plugin.builders.structure.parseTrajectory(modelData, p.modelFormat as 'mmcif');
        model = await plugin.builders.structure.createModel(trajectory);
      }
      const coordData = await plugin.builders.data.download({ url: p.coordinates, isBinary: p.coordinatesBinary });
      const coords = await plugin.dataFormats.get(p.coordinatesFormat)!.parse(plugin, coordData);
      const trajectory = await plugin
        .build()
        .toRoot()
        .apply(
          TrajectoryFromModelAndCoordinates,
          { modelRef: model.ref, coordinatesRef: coords.ref },
          { dependsOn: [model.ref, coords.ref] },
        )
        .commit();
      await plugin.builders.structure.hierarchy.applyPreset(trajectory, p.preset as 'default');
    },
    { rethrowErrors: true },
  );
}

const trajectory: OperationDefinition<typeof TrajectoryParams, Loaded> = {
  id: 'load.trajectory',
  title: 'Load a trajectory',
  params: TrajectoryParams,
  ...network,
  async run(values, ctx) {
    const { model, coordinates } = values;
    if (!model.trim() || !coordinates.trim()) {
      return ApiResult.fail('invalid', 'A trajectory is a model and coordinates.', 'Give both addresses.');
    }
    await ctx.extensions.ensureFor({ formats: [values.modelFormat, values.coordinatesFormat] });
    for (const f of [values.modelFormat, values.coordinatesFormat]) {
      if (!ctx.plugin.dataFormats.has(f)) {
        return ApiResult.fail(
          'invalid',
          `Unknown data format '${f}'.`,
          `Formats: ${ctx.plugin.dataFormats.types.map(([n]) => n).join(', ')}.`,
        );
      }
    }
    const before = refsOf(ctx);
    try {
      await guarded(ctx, () =>
        loadTrajectory(ctx, { ...values, model: model.trim(), coordinates: coordinates.trim() }),
      );
    } catch (e) {
      return pluginFailure(`${model} with ${coordinates}`, e);
    }
    const names = addedNames(ctx, before);
    return ApiResult.ok(
      loadedLine(names, 'the trajectory'),
      { nodes: refsOf(ctx).filter((r) => !before.includes(r)) },
      { value: { loaded: names } },
    );
  },
};

// --- load.volume --------------------------------------------------------------------------------------------------

export interface VolumeIsovalueInfo {
  type: 'absolute' | 'relative';
  value: number;
  color: Color;
  alpha: number;
  volumeIndex: number;
}

export const VolumeParams = {
  url: PD.Text('', { description: 'The address of a volume file, or of a density server cell' }),
  format: PD.Text('ccp4', { description: 'Its format: ccp4, dsn6, cube, dx, dscif, …' }),
  binary: PD.Boolean(true),
  isovalues: PD.ObjectList(
    {
      type: PD.Select<'absolute' | 'relative'>('relative', PD.arrayToOptions(['absolute', 'relative'] as const)),
      value: PD.Numeric(1.5, { step: 0.1 }),
      color: PD.Color(Color(0x3362b2)),
      alpha: PD.Numeric(1, { min: 0, max: 1, step: 0.05 }),
      volumeIndex: PD.Numeric(
        0,
        { min: 0, max: 10, step: 1 },
        { description: "Which of the file's volumes (2Fo-Fc is 0)" },
      ),
    },
    (e) => `${e.type} ${e.value}`,
    {
      defaultValue: [
        { type: 'relative', value: 1.5, color: Color(0x3362b2), alpha: 1, volumeIndex: 0 },
      ] as VolumeIsovalueInfo[],
    },
  ),
  entryId: PD.Text('', { description: 'The entry the volume belongs to (its label)' }),
  lazy: PD.Boolean(false, { description: 'Read only what the view needs (density servers)' }),
};

async function loadVolume(
  ctx: ApiContext,
  p: { url: string; format: string; binary: boolean; isovalues: VolumeIsovalueInfo[]; entryId: string; lazy: boolean },
) {
  const plugin = ctx.plugin;
  await plugin.initialized;
  const entryId = p.entryId || undefined;
  if (p.lazy) {
    const update = plugin.build();
    update.toRoot().apply(LazyVolume, {
      url: p.url,
      format: p.format,
      entryId,
      isBinary: p.binary,
      isovalues: p.isovalues.map((v) => ({ ...v })),
    });
    await update.commit();
    return;
  }
  await plugin.dataTransaction(
    async () => {
      const data = await plugin.builders.data.download(
        { url: p.url, isBinary: p.binary },
        { state: { isGhost: true } },
      );
      const parsed = await plugin.dataFormats.get(p.format)!.parse(plugin, data, { entryId });
      const first = (parsed.volume || parsed.volumes?.[0]) as StateObjectSelector<PluginStateObject.Volume.Data>;
      if (!first?.isOk) throw new Error('Failed to parse any volume.');
      const repr = plugin.build();
      for (const iso of p.isovalues) {
        const volume: StateObjectSelector<PluginStateObject.Volume.Data> =
          parsed.volumes?.[iso.volumeIndex ?? 0] ?? parsed.volume;
        repr.to(volume).apply(
          VolumeRepresentation3D,
          createVolumeRepresentationParams(plugin, first.data!, {
            type: 'isosurface',
            typeParams: {
              alpha: iso.alpha ?? 1,
              isoValue: Volume.adjustedIsoValue(volume.cell!.obj!.data, iso.value, iso.type),
            },
            color: 'uniform',
            colorParams: { value: iso.color },
          }),
        );
      }
      await repr.commit();
    },
    { rethrowErrors: true },
  );
}

const volume: OperationDefinition<typeof VolumeParams, Loaded> = {
  id: 'load.volume',
  title: 'Load a volume',
  params: VolumeParams,
  ...network,
  exists({ url }, ctx) {
    const u = url.trim();
    return u && liveDownloads(ctx).some((d) => d.url === u) ? `${u} is already loaded.` : undefined;
  },
  async run(values, ctx) {
    const address = values.url.trim();
    if (!address) return ApiResult.fail('invalid', 'Which address?', 'Give the URL of a volume file.');
    await ctx.extensions.ensureFor({ formats: [values.format] });
    if (!ctx.plugin.dataFormats.has(values.format)) {
      return ApiResult.fail(
        'invalid',
        `Unknown density format '${values.format}'.`,
        `Formats: ${ctx.plugin.dataFormats.types.map(([n]) => n).join(', ')}.`,
      );
    }
    const cellsBefore = ctx.plugin.state.data.cells.size;
    try {
      await guarded(ctx, () => loadVolume(ctx, { ...values, url: address }));
    } catch (e) {
      return pluginFailure(address, e);
    }
    if (ctx.plugin.state.data.cells.size === cellsBefore) {
      return ApiResult.fail('failed', `Loading ${address} made nothing.`, 'Check the format; the log says more.');
    }
    const label = values.entryId || address;
    return ApiResult.ok(`Loaded the volume ${label}.`, { nodes: [] }, { value: { loaded: [label] } });
  },
};

/** The operations the feature adds to the core's `load` area, under their members' names (`load.entry` is `entry`). */
export const loadersOperations = (): { area: string; name: string; operation: OperationDefinition<any, any> }[] =>
  [entry, url, mvs, snapshot, trajectory, volume].map((operation) => ({
    area: 'load',
    name: memberName(operation.id),
    operation: operation as OperationDefinition<any, any>,
  }));
