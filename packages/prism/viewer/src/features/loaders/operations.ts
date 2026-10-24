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
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { Volume } from '@molstar/model/model/volume';
import { loadMVSData } from '@molstar/mvs/components/formats';
import { PluginConfig } from '@molstar/plugin/config';
import { DownloadStructure, PdbDownloadProvider } from '@molstar/plugin/state/actions/structure';
import { DownloadDensity } from '@molstar/plugin/state/actions/volume';
import type { DataFormatProvider } from '@molstar/plugin/state/formats/provider';
import { CoordinatesFormatCategory } from '@molstar/plugin/state/formats/coordinates/category';
import { ParseXsc } from '@molstar/plugin/state/formats/coordinates/xsc';
import { TopologyFormatCategory } from '@molstar/plugin/state/formats/topology/category';
import { TrajectoryFormatCategory } from '@molstar/plugin/state/formats/trajectory/category';
import { VolumeFormatCategory } from '@molstar/plugin/state/formats/volume/category';
import { createVolumeRepresentationParams } from '@molstar/plugin/state/helpers/volume-representation-params';
import { PluginStateObject } from '@molstar/plugin/state/objects';
import { Download, LazyVolume } from '@molstar/plugin/state/transforms/data/fetch';
import {
  CoordinatesWithBoxes,
  TrajectoryFromModelAndCoordinates,
} from '@molstar/plugin/state/transforms/structure/hierarchy';
import { VolumeRepresentation3D } from '@molstar/plugin/state/transforms/volume/representation';
import { type ApiContext, ApiResult, memberName, type OperationDefinition } from '../../api/area.js';
import {
  asBytes,
  asText,
  contentOf,
  extensionOf,
  fileNameOf,
  type Input,
  InputParams,
  inputModuleProblem,
  inputOf,
  isZip,
} from '../../api/input.js';
import { OptionalParam } from '../../api/params.js';
import { messageOf } from '../../api/run.js';
import { labelOf, loadedStructures } from '../../api/scope.js';
import { type Opened, openSession, SessionOpenParams, sessionAddresses } from '../../api/session.js';
import type { SequenceArchive } from '../../sequence/model.js';
import { chooseFormat, type FormatChoice, type FormatScope, formatNameOf, providerFor } from '../../api/formats.js';

// The `loaders` feature's operations (api.md §2, the `load` row; api.md §5): Mol*'s own loaders
// (`@molstar/plugin-extension/loaders`: the plugin's state actions and builders) as operations of the core's `load`
// area. That package declares `@molstar/plugin-ui` (React), which the viewer must not depend on (api.md §5 table), so
// the few calls are ported here over `@molstar/plugin` directly, with the same state they make. A load that fails
// returns `failed` and the run's checkpoint restores: nothing half-built. Loading something that is loaded says so
// instead (decision 30, api.md §1: identity is the entry id, or the URL).

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
  preset: PD.Text('', {
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
    let out: T;
    try {
      out = await work();
    } catch (e) {
      // A failed download reverts and the work then fails on what is missing; the log has the reason (the fetch).
      if (logged[0]) throw new Error(logged[0]);
      throw e;
    }
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

/** The data object of a download or read. A failed download reverts, so there may be no cell at all. */
function objOf<O>(selector: { cell?: { obj?: O } } | undefined, what: string): O {
  const obj = selector?.cell?.obj;
  if (!obj) throw new Error(`could not read ${what}`);
  return obj;
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

/** "1TQN is already loaded." */
const entriesLoaded = (source: EntrySource, ids: readonly string[]) =>
  `${ids.map((i) => displayId(source, i)).join(', ')} ${ids.length === 1 ? 'is' : 'are'} already loaded.`;

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
    return loaded ? entriesLoaded(source, ids) : undefined;
  },
  async run({ source, id, assembly, preset }, ctx) {
    let representation = preset;
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
      const asked = ids;
      ids = ids.filter((i) => !entryLabels(source, i).some((l) => have.some((d) => d.label === l)));
      // `exists` is asked before the run waits its turn (`api/run.ts`): a load queued behind the same one finds it here.
      if (!ids.length) return ApiResult.already(entriesLoaded(source, asked));
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

// --- load.url, load.data ------------------------------------------------------------------------------------------

/** `binary`, absent: from the extension or the format (api.md §1, "Inputs"). */
const BinaryParam = () =>
  OptionalParam(
    PD.Boolean(false, {
      description:
        'Whether an address or a file is binary; absent: from its extension or the format (bcif, ccp4, dcd, … are)',
    }),
  );

export const UrlParams = {
  url: PD.Text('', { description: 'The address of a structure, volume or other data file' }),
  format: PD.Text('auto', {
    description: "The data format by name (mmcif, bcif, pdb, gro, ccp4, …); 'auto': from the address's extension",
  }),
  binary: BinaryParam(),
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
  data: StateObjectSelector<
    PluginStateObject.Data.String | PluginStateObject.Data.Binary | PluginStateObject.Data.Response
  >,
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

/**
 * The input as a data cell of the plugin's state (a ghost: the parsed nodes are what the scene shows): a download, a
 * file read, or raw data. A download is what `load.url`'s identity reads (`liveDownloads`).
 */
/** The format chosen again with the downloaded content (the format asked for, or the one by the name, otherwise) */
async function withContent(
  ctx: ApiContext,
  choice: FormatChoice,
  o: { format: string; binary: boolean | undefined },
  name: string,
  content: unknown,
): Promise<FormatChoice> {
  // a response is sniffed in the form it holds (its bytes, or its text once a reader decoded it), converting nothing;
  // a text past the string limit is not sniffed, as a String cell's was not
  if (content instanceof PluginStateObject.Data.ResponseBody) content = content.heldBytes ?? content.textNow();
  if (typeof content !== 'string' && !(content instanceof Uint8Array)) return choice;
  const again = await chooseFormat(ctx, { format: o.format, name, binary: choice.binary, content });
  return ApiResult.isFailure(again) ? choice : again;
}

async function dataCellOf(
  ctx: ApiContext,
  input: Input,
  binary: boolean,
  label: string,
): Promise<
  StateObjectSelector<PluginStateObject.Data.String | PluginStateObject.Data.Binary | PluginStateObject.Data.Response>
> {
  const plugin = ctx.plugin;
  const ghost = { state: { isGhost: true } };
  const named = label ? { label } : {};
  switch (input.kind) {
    case 'url':
      return plugin.builders.data.download({ url: input.url, isBinary: binary, ...named }, ghost);
    case 'file':
      return (await plugin.builders.data.readFile({ file: Asset.File(input.file), isBinary: binary, ...named }, ghost))
        .data;
    case 'data':
      return plugin.builders.data.rawData({ data: input.data as string | Uint8Array<ArrayBuffer>, ...named }, ghost);
  }
}

/** "<url> is already loaded.", or undefined. */
function urlLoaded(url: string | undefined, ctx: ApiContext): string | undefined {
  const u = url?.trim();
  return u && liveDownloads(ctx).some((d) => d.url === u) ? `${u} is already loaded.` : undefined;
}

/** The addresses a call reaches (the agent's policy, agent.md §4): its url, if it has one. */
const urlAddresses = ({ url }: { url?: string }): string[] => (url?.trim() ? [url.trim()] : []);

/**
 * `load.url` and `load.data`: the input as a structure (or any data format Mol* reads), its format and whether it is
 * binary told when not given (api.md §1, "Inputs"), with the structure options. One step; a failure leaves
 * nothing.
 */
async function loadContent(
  ctx: ApiContext,
  title: string,
  input: Input,
  o: { format: string; binary: boolean | undefined; label: string; models: 'first' | 'all'; assembly: string },
): Promise<ApiResult<Loaded>> {
  const unavailable = inputModuleProblem(input, ctx, title);
  if (unavailable) return unavailable;
  // Asked again in its turn: one queued behind the same load finds it loaded (`load.entry`).
  const loaded = input.kind === 'url' && !ctx.copy ? urlLoaded(input.url, ctx) : undefined;
  if (loaded) return ApiResult.already(loaded);
  // Data is binary when it is bytes (or base64); a label names it (its extension may tell the format). A download is
  // named by its file (what a particle list or a map is called after), not its whole address.
  const label = input.kind === 'url' ? input.name : o.label.trim();
  const name = label || input.name;
  const choice = await chooseFormat(ctx, {
    format: o.format,
    name,
    binary: input.kind === 'data' ? typeof input.data !== 'string' : o.binary,
    ...(input.kind === 'data' ? { content: input.data } : {}),
  });
  if (ApiResult.isFailure(choice)) return choice;
  const options = { models: o.models, assembly: o.assembly };
  const known = structureOptionsProblem(choice.provider, options);
  if (known) return known;
  const plugin = ctx.plugin;
  const what = input.kind === 'url' ? input.url : name || 'the data';
  const before = refsOf(ctx);
  const cellsBefore = plugin.state.data.cells.size;
  try {
    await plugin.initialized;
    await guarded(ctx, () =>
      plugin.dataTransaction(
        async () => {
          const cell = await dataCellOf(ctx, input, choice.binary, label);
          const obj = objOf(cell, what);
          // A download is told by its content too once it is here: a lazy feature that claims it by a sniff (an
          // mmCIF holding particles) loads, and the format is chosen again with it (viewer-spec.md §3).
          const chosen = input.kind === 'url' ? await withContent(ctx, choice, o, name, obj.data) : choice;
          const provider = providerFor(ctx, chosen, name, obj);
          if (!provider) throw new Error(`could not tell the format of ${what}`);
          const problem = structureOptionsProblem(provider, options);
          if (problem) throw new Error(problem.message);
          await parseAndShow(ctx, cell, provider, options);
        },
        { rethrowErrors: true },
      ),
    );
  } catch (e) {
    if (input.kind === 'url') return pluginFailure(input.url, e);
    return ApiResult.fail('failed', `Could not load ${what}: ${messageOf(e)}`, 'Check the data and its format.');
  }
  if (plugin.state.data.cells.size === cellsBefore) {
    return ApiResult.fail('failed', `Loading ${what} made nothing.`, 'Check the format; the log says more.');
  }
  const names = addedNames(ctx, before);
  return ApiResult.ok(
    loadedLine(names, what),
    { nodes: refsOf(ctx).filter((r) => !before.includes(r)) },
    { value: { loaded: names.length ? names : [what] } },
  );
}

const url: OperationDefinition<typeof UrlParams, Loaded> = {
  id: 'load.url',
  title: 'Load from a URL',
  params: UrlParams,
  addresses: ({ url }) => [url],
  ...network,
  exists: ({ url }, ctx) => urlLoaded(url, ctx),
  async run({ url, format, binary, models, assembly }, ctx) {
    const address = url.trim();
    if (!address) return ApiResult.fail('invalid', 'Which address?', 'Give the URL of a data file.');
    const input: Input = { kind: 'url', url: address, name: fileNameOf(address) };
    return loadContent(ctx, 'Load from a URL', input, { format, binary, label: '', models, assembly });
  },
};

export const DataParams = {
  ...InputParams({
    url: 'The address of a data file',
    data: "The file's contents: text (a PDB or mmCIF file, SDF, …), bytes, or 'base64,…' for a binary format",
    file: 'A structure, volume or other data file',
    accept: '',
  }),
  format: PD.Text('auto', {
    description:
      "The data format by name (mmcif, bcif, pdb, sdf, ccp4, …); 'auto': from the name's extension, else the data",
  }),
  binary: BinaryParam(),
  label: PD.Text('', {
    description: "A name for the data (a file name: '1abc.pdb'); empty: the file's or the format's",
  }),
  ...StructureOptionParams(),
};

const rawData: OperationDefinition<typeof DataParams, Loaded> = {
  id: 'load.data',
  title: 'Load data',
  params: DataParams,
  // Data and a file are local (each call loads: no identity to say "already loaded" by); a url is a download, asked
  // for as `load.url`'s is. The modules are the input's, checked per call.
  approval: 'network',
  addresses: urlAddresses,
  exists: ({ url }, ctx) => urlLoaded(url, ctx),
  async run(values, ctx) {
    const input = inputOf(values, { label: values.label });
    if (!input || ApiResult.isFailure(input)) return input ?? ApiResult.fail('invalid', 'Load what?');
    return loadContent(ctx, 'Load data', input, values);
  },
};

// --- load.mvs -----------------------------------------------------------------------------------------------------

export const MvsParams = {
  ...InputParams({
    url: 'The address of a MolViewSpec file (.mvsj, .mvsx)',
    data: "The MolViewSpec itself: its JSON (mvsj), or an archive's bytes (mvsx; 'base64,…' as text)",
    file: 'A MolViewSpec file (.mvsj, .mvsx)',
    accept: '.mvsj,.mvsx,.json',
  }),
  format: PD.Select<'auto' | 'mvsj' | 'mvsx'>('auto', PD.arrayToOptions(['auto', 'mvsj', 'mvsx'] as const), {
    description: "'auto': from the name's extension (.mvsj, .mvsx), else the content (a zip is mvsx)",
  }),
  keepCamera: PD.Boolean(false, { description: "Keep the camera as it is, ignoring the MolViewSpec's" }),
  // v5's `loadMvs*` options (`MVSLoadOptions`, parity.md).
  append: PD.Boolean(false, {
    description: "Add the MolViewSpec's snapshots to the Snapshots list instead of replacing the list",
  }),
  keepCameraOrientation: PD.Boolean(false, {
    description: "Keep the camera's direction, moving only to the MolViewSpec's target (keepCamera wins)",
  }),
  snapshot: PD.Numeric(0, { min: 0, step: 1 }, { description: 'The snapshot to show first, by index (0: the first)' }),
  extensions: PD.Boolean(true, { description: "Apply Mol*'s MVS extensions (not part of the MolViewSpec standard)" }),
};

/** The kind of a MolViewSpec: given, else the name's extension, else the content (a zip is an archive). */
function mvsFormatOf(format: 'auto' | 'mvsj' | 'mvsx', name: string, content: string | Uint8Array): 'mvsj' | 'mvsx' {
  if (format !== 'auto') return format;
  const ext = extensionOf(name);
  if (ext === 'mvsx') return 'mvsx';
  if (ext === 'mvsj' || ext === 'json') return 'mvsj';
  return isZip(content) ? 'mvsx' : 'mvsj';
}

const mvs: OperationDefinition<typeof MvsParams, Loaded> = {
  id: 'load.mvs',
  title: 'Load a MolViewSpec',
  params: MvsParams,
  // No `modules`: they are the input's, checked per call (a url needs the network; data nothing).
  // It replaces the scene (and an MVS the Snapshots list): destructive wins over network (api/area.ts).
  approval: 'destructive',
  addresses: urlAddresses,
  async run(values, ctx) {
    const { format, keepCamera, append, keepCameraOrientation, snapshot, extensions } = values;
    const input = inputOf(values);
    if (!input || ApiResult.isFailure(input)) return input ?? ApiResult.fail('invalid', 'Load what?');
    const unavailable = inputModuleProblem(input, ctx, 'Load a MolViewSpec');
    if (unavailable) return unavailable;
    const address = input.kind === 'url' ? input.url : '';
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
        ...(address ? { sourceUrl: address } : {}),
      };
      await guarded(ctx, async () => {
        const content = await contentOf(ctx, input);
        const kind = mvsFormatOf(format, input.name, content);
        const data = kind === 'mvsj' ? asText(content) : (asBytes(content) as Uint8Array<ArrayBuffer>);
        await loadMVSData(plugin, data, kind, options);
      });
    } catch (e) {
      return pluginFailure(
        address || input.name || 'the MolViewSpec',
        e,
        'Check the MolViewSpec and the network connection.',
      );
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

export const SvsParams = InputParams({
  url: 'The address of a SeqViewSpec document (.svsj, .svs; JSON or MessagePack) or archive (.svsx)',
  data: "The SeqViewSpec document itself (JSON), or a document's or an archive's bytes ('base64,…' as text)",
  file: 'A SeqViewSpec file: an archive (.svsx) or a document',
  accept: '.svsx,.svsj,.svs',
});

/** An `.svsx` (seqspec:sources#archives): by its name, else a zip's signature. */
const isArchive = (name: string, bytes: Uint8Array) => extensionOf(name) === 'svsx' || isZip(bytes);

/** A file name from a path or address (its last part). */
const baseName = (path: string) => path.split(/[?#]/)[0].split('/').pop() || path;

/** SeqViewSpec structure formats, which the loaders take by these names (`bcif` too); others are read as `auto`. */
const SvsFormats: readonly string[] = ['mmcif', 'bcif', 'pdb', 'sdf', 'mol2'];

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

/** An AlphaFold DB model's entry id (`AF-Q8W3K0-F1`): its accession (lowercase), as a document names it. */
const afdbAccession = (entryId: string | undefined) => (entryId ? /^af-(.+)-f\d+$/.exec(entryId)?.[1] : undefined);

/**
 * The loaded structure a document's source names (sequence.md §6): by its address, entry id (an AlphaFold DB model's
 * by its accession), the label of the entry's download (`load.entry`'s identity), or file name.
 */
function matchStructure(ctx: ApiContext, source: Readonly<Record<string, unknown>>, uri: string | undefined) {
  const entryId = typeof source.id === 'string' ? source.id.toLowerCase() : undefined;
  const service = typeof source.service === 'string' ? SvsServices[source.service] : undefined;
  const labels = service && typeof source.id === 'string' ? entryLabels(service, source.id) : [];
  const hostFile = typeof source.host_file === 'string' ? source.host_file : undefined;
  const fileId = hostFile?.split('/').pop()?.split('.')[0]?.toLowerCase();
  return loadedStructures(ctx).find((s) => {
    const entry = s.cell.obj?.data.model.entryId?.toLowerCase();
    if (entryId && (entry === entryId || afdbAccession(entry) === entryId)) return true;
    if (hostFile && fileId && entry === fileId) return true;
    return dataAncestors(ctx, s.cell.transform.ref).some((cell) => {
      const params = cell.transform.params as
        | { url?: string | Asset.Url; file?: Asset.File; label?: string }
        | undefined;
      if (cell.transform.transformer === Download) {
        if (uri && params?.url && Asset.getUrl(params.url) === uri) return true;
        if (labels.includes((params?.label ?? '').toLowerCase())) return true;
      }
      return !!hostFile && (params?.file?.name === hostFile || cell.obj?.label === hostFile);
    });
  });
}

/** A file the archive holds at a relative URI (against its main file's directory), as a `File` of its name. */
function archivedFile(archive: SequenceArchive, uri: string): File | undefined {
  if (/^[a-z][a-z0-9+.-]*:/i.test(uri) || uri.startsWith('/')) return undefined;
  const parts: string[] = archive.base ? archive.base.split('/') : [];
  for (const p of uri.split('/')) {
    if (p === '..') {
      if (!parts.length) return undefined;
      parts.pop();
    } else if (p && p !== '.') parts.push(p);
  }
  const bytes = archive.files.get(parts.join('/'));
  return bytes ? new File([bytes as Uint8Array<ArrayBuffer>], baseName(uri)) : undefined;
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
  // No `modules`: they are the input's, checked per call (data is local, and a structure that is loaded already needs
  // nothing).
  approval: 'network',
  async run(values, ctx) {
    const input = inputOf(values);
    if (!input || ApiResult.isFailure(input)) return input ?? ApiResult.fail('invalid', 'Load what?');
    const sequence = ctx.api.sequence;
    if (!sequence) {
      return ApiResult.fail('unavailable', 'This viewer has no sequence view to show a SeqViewSpec in.');
    }
    const unavailable = inputModuleProblem(input, ctx, 'Load a SeqViewSpec');
    if (unavailable) return unavailable;
    const address = input.kind === 'url' ? input.url : '';
    // The document: JSON or MessagePack from an address or a file, JSON as data, or an archive's (`.svsx`: the document
    // and the files it reads); its envelope checked (Seq* checks the rest). Seq*'s readers come with the sequence's
    // runtime (its own chunk).
    const { fromJSON, readDocument, readSvsx, svsStructureSource } = await sequence.runtime();
    let read: ReturnType<typeof readDocument>;
    let archive: SequenceArchive | undefined;
    try {
      if (input.kind === 'data' && typeof input.data === 'string') read = fromJSON(input.data.trim());
      else {
        const bytes = asBytes(await contentOf(ctx, input));
        if (isArchive(input.name, bytes)) {
          const x = await readSvsx(bytes);
          if (!x.doc) {
            return ApiResult.fail(
              'invalid',
              `Not a SeqViewSpec archive: ${x.issues[0]?.message ?? 'it cannot be read'}.`,
              'An .svsx is a zip with index.svsm or index.svsj and the files it reads.',
            );
          }
          read = { doc: x.doc, issues: x.issues };
          archive = { files: x.files, base: x.base };
        } else read = readDocument(bytes);
      }
    } catch (e) {
      return pluginFailure(address || input.name || 'the SeqViewSpec', e);
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
    // A structure in the archive (a relative URI it holds) opens as a file of that name.
    const archived = archive && typeof src.uri === 'string' ? archivedFile(archive, src.uri) : undefined;
    const uri = typeof src.uri === 'string' && !archived ? resolveUri(src.uri, address) : undefined;
    const before = refsOf(ctx);
    let structure = matchStructure(ctx, src, uri);
    if (!structure) {
      // Its format by name, else told from the address, the name or the data (`auto`).
      const format = SvsFormats.includes(String(src.format)) ? { format: String(src.format) } : {};
      let loaded: ApiResult<Loaded>;
      if (archived) {
        loaded = (await ctx.api.run('load.files', { files: [archived] })) as ApiResult<Loaded>;
      } else if (service && typeof src.id === 'string') {
        loaded = (await ctx.api.run('load.entry', { source: service, id: src.id })) as ApiResult<Loaded>;
      } else if (uri) {
        loaded = (await ctx.api.run('load.url', { url: uri, ...format })) as ApiResult<Loaded>;
      } else if (typeof src.value === 'string') {
        loaded = (await ctx.api.run('load.data', { data: src.value, ...format })) as ApiResult<Loaded>;
      } else {
        const name = typeof src.host_file === 'string' ? src.host_file : 'its structure';
        return ApiResult.fail(
          'not-found',
          `The SeqViewSpec's structure (${name}) is not loaded.`,
          'Open the file first, or with the document.',
        );
      }
      if (!loaded.ok) return loaded;
      // Already loaded (the nested load said so, by its own identity): the structure is the one that matches it.
      structure =
        loadedStructures(ctx).find((s) => !before.includes(s.cell.transform.ref)) ??
        (loaded.unchanged ? matchStructure(ctx, src, uri) : undefined);
      if (!structure) return ApiResult.fail('failed', "Loading the SeqViewSpec's structure made none.");
    }
    const shown = await sequence.showDocument(
      { document: doc, structure: structure.cell.transform.ref, ...(address ? { url: address } : {}) },
      archive,
    );
    if (!shown.ok) {
      const issue = shown.issues.find((i) => i.severity === 'error') ?? shown.issues[0];
      return ApiResult.fail(
        'invalid',
        `The sequence cannot show the SeqViewSpec: ${issue?.message ?? 'it is not valid'}.`,
        'Check the document.',
      );
    }
    const names = addedNames(ctx, before);
    const title = doc.metadata?.title ?? (address || input.name || 'the SeqViewSpec');
    return ApiResult.ok(
      `${names.length ? `Loaded ${names.join(', ')}; the` : 'The'} sequence shows ${title} over ${labelOf(ctx, structure)}.`,
      { nodes: refsOf(ctx).filter((r) => !before.includes(r)), sequence: true },
      { value: { loaded: names } },
    );
  },
};

// --- load.snapshot ------------------------------------------------------------------------------------------------

/**
 * `session.open` under the loaders' name (api.md §1, "Inputs"): the same params and the same run (`openSession`), for
 * the app's `snapshot-url` and `snapshot-id` URL options and v5's scripts. New code opens sessions with `session.open`.
 */
const snapshot: OperationDefinition<typeof SessionOpenParams, Opened> = {
  id: 'load.snapshot',
  title: 'Load a session (as session.open)',
  params: SessionOpenParams,
  // As `session.open`: the modules are the input's, checked per call; it replaces the scene (destructive).
  approval: 'destructive',
  addresses: sessionAddresses,
  run: openSession,
};

// --- load.trajectory ----------------------------------------------------------------------------------------------

export const TrajectoryParams = {
  model: PD.Text('', { description: 'The address of the model (or topology) file' }),
  modelFormat: PD.Text('auto', {
    description: "Its format: mmcif, pdb, gro, …, or a topology (psf, prmtop, top); 'auto': from its extension",
  }),
  modelBinary: BinaryParam(),
  coordinates: PD.Text('', { description: 'The address of the coordinates file' }),
  coordinatesFormat: PD.Text('auto', {
    description: "Its format: xtc, dcd, trr, nctraj, lammpstrj; 'auto': from its extension",
  }),
  coordinatesBinary: BinaryParam(),
  boxes: PD.Text('', {
    description:
      "The address of a NAMD extended system file (.xsc, .xst): each frame's periodic box with its origin, which a DCD lacks; empty: none",
  }),
  preset: PD.Text('default', { description: "The structure preset: 'default' or 'all-models'" }),
};

const ModelFormats: FormatScope = {
  categories: [TrajectoryFormatCategory, TopologyFormatCategory],
  what: 'a model or topology',
};
const CoordinatesFormats: FormatScope = { categories: [CoordinatesFormatCategory], what: 'a coordinates' };

async function loadTrajectory(
  ctx: ApiContext,
  p: {
    model: string;
    modelFormat: FormatChoice;
    coordinates: string;
    coordinatesFormat: FormatChoice;
    boxes?: string;
    preset: string;
  },
) {
  const plugin = ctx.plugin;
  await plugin.initialized;
  await plugin.dataTransaction(
    async () => {
      // Named by their files (not their whole addresses): the structure's row says protein.pdb.
      const modelName = fileNameOf(p.model) || p.model;
      const modelData = await plugin.builders.data.download({
        url: p.model,
        isBinary: p.modelFormat.binary,
        label: modelName,
      });
      const modelProvider = providerFor(ctx, p.modelFormat, modelName, objOf(modelData, modelName));
      if (!modelProvider) throw new Error(`could not tell the format of ${modelName}`);
      let model: StateObjectSelector;
      if (modelProvider.category === TopologyFormatCategory) {
        model = (await modelProvider.parse(plugin, modelData)).topology;
      } else {
        const trajectory = await plugin.builders.structure.parseTrajectory(modelData, modelProvider as never);
        model = await plugin.builders.structure.createModel(trajectory);
      }
      const coordinatesName = fileNameOf(p.coordinates) || p.coordinates;
      const coordData = await plugin.builders.data.download({
        url: p.coordinates,
        isBinary: p.coordinatesFormat.binary,
        label: coordinatesName,
      });
      const coordinatesProvider = providerFor(
        ctx,
        p.coordinatesFormat,
        coordinatesName,
        objOf(coordData, coordinatesName),
      );
      if (!coordinatesProvider) throw new Error(`could not tell the format of ${coordinatesName}`);
      let coords: StateObjectSelector = await coordinatesProvider.parse(plugin, coordData);
      // the frames' boxes from a NAMD extended system file: their origins and orientation
      if (p.boxes) {
        const boxesName = fileNameOf(p.boxes) || p.boxes;
        const boxesData = await plugin.builders.data.download({ url: p.boxes, isBinary: false, label: boxesName });
        const boxes = await plugin.build().to(boxesData).apply(ParseXsc).commit();
        coords = await plugin
          .build()
          .to(coords)
          .apply(CoordinatesWithBoxes, { boxesRef: boxes.ref }, { dependsOn: [boxes.ref] })
          .commit();
      }
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
  addresses: ({ model, coordinates, boxes }) => [model, coordinates, ...(boxes?.trim() ? [boxes.trim()] : [])],
  ...network,
  async run(values, ctx) {
    const model = values.model.trim();
    const coordinates = values.coordinates.trim();
    if (!model || !coordinates) {
      return ApiResult.fail('invalid', 'A trajectory is a model and coordinates.', 'Give both addresses.');
    }
    const modelFormat = await chooseFormat(ctx, {
      format: values.modelFormat,
      name: fileNameOf(model),
      binary: values.modelBinary,
      scope: ModelFormats,
    });
    if (ApiResult.isFailure(modelFormat)) return modelFormat;
    const coordinatesFormat = await chooseFormat(ctx, {
      format: values.coordinatesFormat,
      name: fileNameOf(coordinates),
      binary: values.coordinatesBinary,
      scope: CoordinatesFormats,
    });
    if (ApiResult.isFailure(coordinatesFormat)) return coordinatesFormat;
    const before = refsOf(ctx);
    try {
      await guarded(ctx, () =>
        loadTrajectory(ctx, {
          model,
          modelFormat,
          coordinates,
          coordinatesFormat,
          boxes: values.boxes?.trim() || undefined,
          preset: values.preset,
        }),
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
  ...InputParams({
    url: 'The address of a volume file, or of a density server cell',
    data: "The volume file's contents: bytes, 'base64,…' for a binary format, or text (cube, dx)",
    file: 'A volume file (.ccp4, .mrc, .map, .dsn6, .cube, .dx, …)',
    accept: '.ccp4,.mrc,.map,.dsn6,.brix,.cube,.cub,.dx,.dxbin,.bcif,.cif',
  }),
  format: PD.Text('auto', {
    description:
      "Its format: ccp4 (mrc, map), dsn6, cube, dx, dscif, …; 'auto': from the name's extension, else the data",
  }),
  binary: BinaryParam(),
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
  lazy: PD.Boolean(false, { description: 'Read only what the view needs (density servers; an address only)' }),
};

const VolumeFormats: FormatScope = { categories: [VolumeFormatCategory], what: 'a volume' };

async function loadVolume(
  ctx: ApiContext,
  input: Input,
  choice: FormatChoice,
  p: { isovalues: VolumeIsovalueInfo[]; entryId: string; lazy: boolean },
) {
  const plugin = ctx.plugin;
  await plugin.initialized;
  const entryId = p.entryId || undefined;
  if (p.lazy && input.kind === 'url') {
    // The format by name: the one given, else the extension's.
    const format = formatNameOf(choice);
    const update = plugin.build();
    update.toRoot().apply(LazyVolume, {
      url: input.url,
      format,
      entryId,
      isBinary: choice.binary,
      isovalues: p.isovalues.map((v) => ({ ...v })),
    });
    await update.commit();
    return;
  }
  await plugin.dataTransaction(
    async () => {
      const data = await dataCellOf(ctx, input, choice.binary, '');
      const provider = providerFor(ctx, choice, input.name, objOf(data, input.name || 'the data'));
      if (!provider) throw new Error(`could not tell the format of ${input.name || 'the data'}`);
      const parsed = await provider.parse(plugin, data, { entryId });
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
  // Data and a file are local; a url is a download. The modules are the input's, checked per call.
  approval: 'network',
  addresses: urlAddresses,
  exists: ({ url }, ctx) => urlLoaded(url, ctx),
  async run(values, ctx) {
    const input = inputOf(values);
    if (!input || ApiResult.isFailure(input)) return input ?? ApiResult.fail('invalid', 'Load what?');
    const unavailable = inputModuleProblem(input, ctx, 'Load a volume');
    if (unavailable) return unavailable;
    if (values.lazy && input.kind !== 'url') {
      return ApiResult.fail('invalid', 'A lazy volume is read from an address.', "Give its url, or leave out 'lazy'.");
    }
    const loaded = input.kind === 'url' && !ctx.copy ? urlLoaded(input.url, ctx) : undefined;
    if (loaded) return ApiResult.already(loaded);
    const choice = await chooseFormat(ctx, {
      format: values.format,
      name: input.name,
      binary: input.kind === 'data' ? typeof input.data !== 'string' : values.binary,
      ...(input.kind === 'data' ? { content: input.data } : {}),
      scope: VolumeFormats,
    });
    if (ApiResult.isFailure(choice)) return choice;
    const what = input.kind === 'url' ? input.url : input.name || 'the data';
    const cellsBefore = ctx.plugin.state.data.cells.size;
    try {
      await guarded(ctx, () => loadVolume(ctx, input, choice, values));
    } catch (e) {
      return pluginFailure(what, e, input.kind === 'url' ? undefined : 'Check the data and its format.');
    }
    if (ctx.plugin.state.data.cells.size === cellsBefore) {
      return ApiResult.fail('failed', `Loading ${what} made nothing.`, 'Check the format; the log says more.');
    }
    const label = values.entryId || what;
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
