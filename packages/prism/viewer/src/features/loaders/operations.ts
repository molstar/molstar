/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type StateObjectCell, StateObjectSelector } from '@molstar/core/state';
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
import type { DataFormatProvider } from '@molstar/plugin/state/formats/provider';
import { TopologyFormatCategory } from '@molstar/plugin/state/formats/topology/category';
import { TrajectoryFormatCategory } from '@molstar/plugin/state/formats/trajectory/category';
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
  // The known data services only (the database's own address).
  addresses: () => [],
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
  ...StructureOptionParams(),
};

/**
 * How a structure file is built (v5's `loadAllModelsOrAssemblyFromUrl`): every model, and an assembly. Structure
 * formats only; the defaults are the format's own load (the first model, its default assembly).
 */
function StructureOptionParams() {
  return {
    models: PD.Select<'first' | 'all'>(
      'first',
      [
        ['first', 'The first model'],
        ['all', 'Every model (one model: as first)'],
      ],
      {
        description: "Structures: 'all' loads every model of an ensemble or a trajectory file, each its own structure",
      },
    ),
    assembly: PD.Text('', { description: "Structures: the assembly id to build ('1', …); empty: the default" }),
  };
}

interface StructureOptions {
  readonly models: 'first' | 'all';
  readonly assembly: string;
}

const hasStructureOptions = (o: StructureOptions) => o.models === 'all' || !!o.assembly.trim();

/**
 * Parses the data and shows it: the format's own visuals, or, with structure options, the trajectory built by the
 * hierarchy presets as v5's loaders do (`all-models` falls back to the default for a single model).
 */
async function parseAndShow(
  ctx: ApiContext,
  data: StateObjectSelector<PluginStateObject.Data.String | PluginStateObject.Data.Binary>,
  provider: DataFormatProvider,
  options: StructureOptions,
) {
  const plugin = ctx.plugin;
  if (!hasStructureOptions(options)) {
    const parsed = await provider.parse(plugin, data);
    await provider.visuals?.(plugin, parsed);
    return;
  }
  const trajectory = await plugin.builders.structure.parseTrajectory(data, provider as never);
  const assembly = options.assembly.trim();
  const structure = assembly ? { name: 'assembly' as const, params: { id: assembly } } : undefined;
  if (options.models === 'all') {
    await plugin.builders.structure.hierarchy.applyPreset(trajectory, 'all-models', {
      useDefaultIfSingleModel: true,
      ...(structure ? { structure } : {}),
    } as never);
  } else {
    await plugin.builders.structure.hierarchy.applyPreset(trajectory, 'default', { structure });
  }
}

/** Why the structure options do not apply to a format, or undefined. */
function structureOptionsProblem(provider: DataFormatProvider | undefined, options: StructureOptions) {
  if (!provider || !hasStructureOptions(options) || provider.category === TrajectoryFormatCategory) return undefined;
  return ApiResult.fail(
    'invalid',
    `Models and assemblies are for structures; ${provider.label} is not one.`,
    "Leave out 'models' and 'assembly'.",
  );
}

async function loadUrl(ctx: ApiContext, url: string, format: string, isBinary: boolean, options: StructureOptions) {
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
      const problem = structureOptionsProblem(provider, options);
      if (problem) throw new Error(problem.message);
      await parseAndShow(ctx, data, provider, options);
    },
    { rethrowErrors: true },
  );
}

/** A format by name, its lazy extension loaded first; a failure names the formats. */
async function formatProblem(ctx: ApiContext, format: string) {
  if (format === 'auto') return undefined;
  await ctx.extensions.ensureFor({ formats: [format] });
  if (ctx.plugin.dataFormats.has(format)) return undefined;
  return ApiResult.fail(
    'invalid',
    `Unknown data format '${format}'.`,
    `Formats: auto, ${ctx.plugin.dataFormats.types.map(([n]) => n).join(', ')}.`,
  );
}

const url: OperationDefinition<typeof UrlParams, Loaded> = {
  id: 'load.url',
  title: 'Load from a URL',
  params: UrlParams,
  addresses: ({ url }) => [url],
  ...network,
  exists({ url }, ctx) {
    const u = url.trim();
    return u && liveDownloads(ctx).some((d) => d.url === u) ? `${u} is already loaded.` : undefined;
  },
  async run({ url, format, binary, models, assembly }, ctx) {
    const address = url.trim();
    if (!address) return ApiResult.fail('invalid', 'Which address?', 'Give the URL of a data file.');
    format = format.trim() || 'mmcif';
    const unknown = await formatProblem(ctx, format);
    if (unknown) return unknown;
    const options = { models, assembly };
    const known = structureOptionsProblem(format === 'auto' ? undefined : ctx.plugin.dataFormats.get(format), options);
    if (known) return known;
    const before = refsOf(ctx);
    try {
      await guarded(ctx, () => loadUrl(ctx, address, format, binary, options));
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

// --- load.data ----------------------------------------------------------------------------------------------------

export const DataParams = {
  data: PD.Text('', {
    multiline: true,
    description: "The file's contents: text (a PDB or mmCIF file, SDF, …), or 'base64,…' for a binary format",
  }),
  format: PD.Text('mmcif', {
    description:
      "The data format by name (mmcif, pdb, sdf, bcif, ccp4, …), or 'auto' to read it from the label's extension",
  }),
  label: PD.Text('', { description: "A name for the data (a file name: '1abc.pdb'); empty: the format's" }),
  /** Scripts: the bytes of a binary file, instead of `data` (not for the command's form or the agent). */
  bytes: PD.Value<Uint8Array | ArrayBuffer | number[] | undefined>(undefined, { isHidden: true }),
  ...StructureOptionParams(),
};

const Base64 = /^base64,/i;

/** The data as the plugin's raw data takes it: text, or bytes (given, or decoded from `base64,…`). */
function rawOf(data: string, bytes: Uint8Array | ArrayBuffer | number[] | undefined): string | Uint8Array | undefined {
  if (bytes !== undefined) {
    if (bytes instanceof Uint8Array) return bytes;
    return new Uint8Array(bytes as ArrayBuffer);
  }
  if (!data) return undefined;
  if (!Base64.test(data)) return data;
  const text = atob(data.replace(Base64, '').replace(/\s+/g, ''));
  const out = new Uint8Array(text.length);
  for (let i = 0; i < text.length; i++) out[i] = text.charCodeAt(i);
  return out;
}

const rawData: OperationDefinition<typeof DataParams, Loaded> = {
  id: 'load.data',
  title: 'Load data',
  params: DataParams,
  // Local: nothing is fetched. Each call loads (there is no identity to say "already loaded" by).
  async run({ data, format, label, bytes, models, assembly }, ctx) {
    if (data.trim() && bytes !== undefined) {
      return ApiResult.fail('invalid', 'Give data or bytes, not both.', 'Neither is loaded: say which one is meant.');
    }
    let raw: string | Uint8Array | undefined;
    try {
      raw = rawOf(data.trim() ? data : '', bytes);
    } catch {
      return ApiResult.fail('invalid', 'The base64 data does not decode.', "Give 'base64,' and the file's bytes.");
    }
    if (raw === undefined || (typeof raw !== 'string' && raw.length === 0)) {
      return ApiResult.fail('invalid', 'Which data?', "Give the file's contents as data.");
    }
    format = format.trim() || 'mmcif';
    if (format === 'auto' && !getFileNameInfo(label.trim()).ext) {
      return ApiResult.fail('invalid', 'Which format?', "Give the format, or a label with the file's extension.");
    }
    const unknown = await formatProblem(ctx, format);
    if (unknown) return unknown;
    const options = { models, assembly };
    const known = structureOptionsProblem(format === 'auto' ? undefined : ctx.plugin.dataFormats.get(format), options);
    if (known) return known;
    const plugin = ctx.plugin;
    const name = label.trim();
    const before = refsOf(ctx);
    const cellsBefore = plugin.state.data.cells.size;
    try {
      await plugin.initialized;
      await guarded(ctx, () =>
        plugin.dataTransaction(
          async () => {
            const cell = await plugin.builders.data.rawData(
              { data: raw as string | Uint8Array<ArrayBuffer>, ...(name ? { label: name } : {}) },
              { state: { isGhost: true } },
            );
            const provider =
              format === 'auto'
                ? plugin.dataFormats.auto(getFileNameInfo(name), cell.cell!.obj!)
                : plugin.dataFormats.get(format);
            if (!provider) throw new Error(`could not tell the format of ${name}`);
            const problem = structureOptionsProblem(provider, options);
            if (problem) throw new Error(problem.message);
            await parseAndShow(ctx, cell, provider, options);
          },
          { rethrowErrors: true },
        ),
      );
    } catch (e) {
      return ApiResult.fail('failed', `Could not load the data: ${messageOf(e)}`, 'Check the data and its format.');
    }
    if (plugin.state.data.cells.size === cellsBefore) {
      return ApiResult.fail('failed', 'Loading the data made nothing.', 'Check the format; the log says more.');
    }
    const names = addedNames(ctx, before);
    const what = name || 'the data';
    return ApiResult.ok(
      loadedLine(names, what),
      { nodes: refsOf(ctx).filter((r) => !before.includes(r)) },
      { value: { loaded: names.length ? names : [what] } },
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
  // v5's `loadMvs*` options (`MVSLoadOptions`, parity review M6).
  append: PD.Boolean(false, {
    description: "Add the MolViewSpec's snapshots to the Snapshots list instead of replacing the list",
  }),
  keepCameraOrientation: PD.Boolean(false, {
    description: "Keep the camera's direction, moving only to the MolViewSpec's target (keepCamera wins)",
  }),
  snapshot: PD.Numeric(0, { min: 0, step: 1 }, { description: 'The snapshot to show first, by index (0: the first)' }),
  extensions: PD.Boolean(true, { description: "Apply Mol*'s MVS extensions (not part of the MolViewSpec standard)" }),
};

const mvs: OperationDefinition<typeof MvsParams, Loaded> = {
  id: 'load.mvs',
  title: 'Load a MolViewSpec',
  params: MvsParams,
  // No `modules: ['network']`: the module is checked per call, as only a url needs it (data is local).
  // It replaces the scene (and an MVS the Snapshots list): destructive wins over network (api/area.ts).
  approval: 'destructive',
  async run({ url, data, format, keepCamera, append, keepCameraOrientation, snapshot, extensions }, ctx) {
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
      const options = {
        sanityChecks: true,
        keepCamera,
        keepCameraOrientation,
        appendSnapshots: append,
        defaultSnapshotIndex: Math.max(0, Math.floor(snapshot)),
        ...(extensions ? {} : { extensions: [] }),
      };
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

// --- load.svs -----------------------------------------------------------------------------------------------------

export const SvsParams = {
  url: PD.Text('', {
    description: 'The address of a SeqViewSpec document (.svsj, .svs; JSON or MessagePack); empty: use data',
  }),
  data: PD.Text('', {
    multiline: true,
    description: 'The SeqViewSpec document itself (JSON); used when there is no url',
  }),
};

/** SeqViewSpec structure formats as the plugin's loaders name them (BinaryCIF is mmCIF read as binary). */
const SvsFormats: Readonly<Record<string, { format: string; binary: boolean }>> = {
  mmcif: { format: 'mmcif', binary: false },
  bcif: { format: 'mmcif', binary: true },
  pdb: { format: 'pdb', binary: false },
  sdf: { format: 'sdf', binary: false },
  mol2: { format: 'mol2', binary: false },
};

/** The entry-id services a document may name its structure by, as `load.entry`'s sources. */
const SvsServices: Readonly<Record<string, EntrySource>> = { pdbe: 'pdb', afdb: 'afdb' };

/** The state refs of the data cells above a cell (the cell itself first). */
function dataAncestors(ctx: ApiContext, ref: string): StateObjectCell[] {
  const { cells } = ctx.plugin.state.data;
  const out: StateObjectCell[] = [];
  for (let cell = cells.get(ref); cell; ) {
    out.push(cell);
    const parent = cell.transform.parent;
    if (parent === cell.transform.ref) break;
    cell = cells.get(parent);
  }
  return out;
}

/** The loaded structure a document's source names (sequence.md §6): by its address, entry id, or file name. */
function matchStructure(ctx: ApiContext, source: Readonly<Record<string, unknown>>, uri: string | undefined) {
  const entryId = typeof source.id === 'string' ? source.id.toLowerCase() : undefined;
  const hostFile = typeof source.host_file === 'string' ? source.host_file : undefined;
  const fileId = hostFile?.split('/').pop()?.split('.')[0]?.toLowerCase();
  return loadedStructures(ctx).find((s) => {
    const entry = s.cell.obj?.data.model.entryId?.toLowerCase();
    if (entryId && entry === entryId) return true;
    if (hostFile && fileId && entry === fileId) return true;
    return dataAncestors(ctx, s.cell.transform.ref).some((cell) => {
      const params = cell.transform.params as { url?: string | Asset.Url; file?: Asset.File } | undefined;
      if (uri && cell.transform.transformer === Download && params?.url && Asset.getUrl(params.url) === uri)
        return true;
      return !!hostFile && (params?.file?.name === hostFile || cell.obj?.label === hostFile);
    });
  });
}

/** An address resolved against the document's (a relative URI in a document from an address). */
function resolveUri(uri: string, base: string): string {
  try {
    return base ? new URL(uri, base).href : uri;
  } catch {
    return uri;
  }
}

const svs: OperationDefinition<typeof SvsParams, Loaded> = {
  id: 'load.svs',
  title: 'Load a SeqViewSpec',
  params: SvsParams,
  // No `addresses`: the document names a structure that is fetched in the same run, which a check of the document's
  // address alone cannot see, so the agent always asks (agent.md §4).
  // No `modules: ['network']`: data is local, and a structure that is loaded already needs nothing.
  approval: 'network',
  async run({ url, data }, ctx) {
    const address = url.trim();
    if (address && data.trim()) {
      return ApiResult.fail('invalid', 'Give a url or data, not both.', 'Neither is loaded: say which one is meant.');
    }
    if (!address && !data.trim()) return ApiResult.fail('invalid', 'Which SeqViewSpec?', 'Give a url or its data.');
    const sequence = ctx.api.sequence;
    if (!sequence) {
      return ApiResult.fail('unavailable', 'This viewer has no sequence view to show a SeqViewSpec in.');
    }
    if (address && !ctx.modules.has('network')) {
      return ApiResult.fail(
        'unavailable',
        'Load a SeqViewSpec needs the network, which this viewer does not have.',
        'Give its data instead of a url.',
      );
    }
    if (/\.svsx$/i.test(address.split(/[?#]/)[0])) {
      return ApiResult.fail('invalid', 'SeqViewSpec archives (.svsx) are not read yet.', 'Give its .svsj document.');
    }
    // The document: JSON or MessagePack from an address, JSON as data; its envelope checked (Seq* checks the rest).
    // Seq*'s readers come with the sequence's runtime (its own chunk).
    const { fromJSON, readDocument, svsStructureSource } = await sequence.runtime();
    let read: ReturnType<typeof readDocument>;
    try {
      read = address
        ? readDocument(await ctx.plugin.runTask(ctx.plugin.fetch({ url: address, type: 'binary' })))
        : fromJSON(data.trim());
    } catch (e) {
      return pluginFailure(address, e);
    }
    const doc = read.doc;
    if (!doc) {
      return ApiResult.fail(
        'invalid',
        `Not a SeqViewSpec: ${read.issues[0]?.message ?? 'the document cannot be read'}.`,
        'Give a SeqViewSpec document (format "seq_view_spec").',
      );
    }
    const target = svsStructureSource(doc);
    if (!target) {
      return ApiResult.fail(
        'invalid',
        'The SeqViewSpec names no structure (a molstar:structure complex with a source).',
        'The viewer shows a SeqViewSpec over a structure it loads.',
      );
    }
    // Its first molstar:structure complex: the loaded structure it names, else loaded now (one step with the SVS).
    const src = target.source;
    const service = typeof src.service === 'string' ? SvsServices[src.service] : undefined;
    const uri = typeof src.uri === 'string' ? resolveUri(src.uri, address) : undefined;
    const before = refsOf(ctx);
    let structure = matchStructure(ctx, src, uri);
    if (!structure) {
      const format = SvsFormats[String(src.format ?? '')] ?? { format: 'auto', binary: false };
      let loaded: ApiResult<Loaded>;
      if (service && typeof src.id === 'string') {
        loaded = (await ctx.api.run('load.entry', { source: service, id: src.id })) as ApiResult<Loaded>;
      } else if (uri) {
        loaded = (await ctx.api.run('load.url', { url: uri, ...format })) as ApiResult<Loaded>;
      } else if (typeof src.value === 'string') {
        const text = format.format === 'auto' ? 'mmcif' : format.format;
        loaded = (await ctx.api.run('load.data', { data: src.value, format: text })) as ApiResult<Loaded>;
      } else {
        const name = typeof src.host_file === 'string' ? src.host_file : 'its structure';
        return ApiResult.fail(
          'not-found',
          `The SeqViewSpec's structure (${name}) is not loaded.`,
          'Open the file first, or with the document.',
        );
      }
      if (!loaded.ok) return loaded;
      structure = loadedStructures(ctx).find((s) => !before.includes(s.cell.transform.ref));
      if (!structure) return ApiResult.fail('failed', "Loading the SeqViewSpec's structure made none.");
    }
    const shown = await sequence.showDocument({
      document: doc,
      structure: structure.cell.transform.ref,
      ...(address ? { url: address } : {}),
    });
    if (!shown.ok) {
      const issue = shown.issues.find((i) => i.severity === 'error') ?? shown.issues[0];
      return ApiResult.fail(
        'invalid',
        `The sequence cannot show the SeqViewSpec: ${issue?.message ?? 'it is not valid'}.`,
        'Check the document.',
      );
    }
    const names = addedNames(ctx, before);
    const title = doc.metadata?.title ?? address ?? 'the SeqViewSpec';
    return ApiResult.ok(
      `${names.length ? `Loaded ${names.join(', ')}; the` : 'The'} sequence shows ${title} over ${labelOf(ctx, structure)}.`,
      { nodes: refsOf(ctx).filter((r) => !before.includes(r)), sequence: true },
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
  // An id is the configured state server's.
  addresses: ({ url }) => (url ? [url] : []),
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
  addresses: ({ model, coordinates }) => [model, coordinates],
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
  addresses: ({ url }) => [url],
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
  [entry, url, rawData, mvs, svs, snapshot, trajectory, volume].map((operation) => ({
    area: 'load',
    name: memberName(operation.id),
    operation: operation as OperationDefinition<any, any>,
  }));
