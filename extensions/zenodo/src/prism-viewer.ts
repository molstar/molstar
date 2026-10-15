/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { getFileNameInfo } from '@molstar/core/util/file-info';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import type { PluginContext } from '@molstar/plugin/context';
import { CoordinatesFormatCategory } from '@molstar/plugin/state/formats/coordinates/category';
import { TopologyFormatCategory } from '@molstar/plugin/state/formats/topology/category';
import { TrajectoryFormatCategory } from '@molstar/plugin/state/formats/trajectory/category';
import { VolumeFormatCategory } from '@molstar/plugin/state/formats/volume/category';
import { ApiArea, type ApiContext, type ApiFailure, ApiResult } from '@molstar/prism-viewer/api/area';
import { messageOf } from '@molstar/prism-viewer/api/run';
import { operationCommand } from '@molstar/prism-viewer/commands/operations';
import type { ViewerExtension } from '@molstar/prism-viewer/spec';
import {
  ZenodoCommands,
  ZenodoOperations,
  ZenodoSection,
  ZenodoSkill,
  zenodoFileAddress,
  zenodoRecordAddress,
  zenodoRecordId,
} from './prism-viewer-stub.js';

declare module '@molstar/prism-viewer/api' {
  // The `zenodo` area: the extension's own (operations.md §2, extensions add areas).
  interface ExtensionApis {
    zenodo: ZenodoApi;
  }
}

// Zenodo in the Prism viewer (plan 5.2, parity.md §6: Home's import slot; lazy): v5's Zenodo Import as two operations.
// `zenodo.list` reads a record's files from Zenodo's API and says which the viewer reads (by the plugin's data formats,
// as v5 did); `zenodo.load` loads one through the loaders feature's `load.url`, `load.volume` or `load.trajectory`, as
// part of its own step (one Undo). Both reach only zenodo.org, which is what their `addresses` say: the agent's policy
// asks unless the person allowed that domain (agent.md §4).

/** What the viewer reads a file as, by its extension. */
export type ZenodoFileKind = 'structure' | 'topology' | 'coordinates' | 'volume';

export interface ZenodoFileInfo {
  readonly name: string;
  readonly size: number;
  readonly url: string;
  /** Absent: the viewer does not read it. */
  readonly kind?: ZenodoFileKind;
  readonly format?: string;
  readonly binary?: boolean;
}

export interface ZenodoRecordInfo {
  readonly id: string;
  readonly title: string;
  readonly doi?: string;
  readonly files: readonly ZenodoFileInfo[];
}

/** The part of Zenodo's record JSON that is read (the API's current shape and the older one: `type`, `bucket`). */
interface RecordJson {
  id?: number | string;
  doi?: string;
  metadata?: { title?: string; doi?: string };
  files?: { key?: string; filename?: string; size?: number; type?: string; links?: Record<string, string> }[];
}

const RecordParam = PD.Text('', {
  description:
    'The record: its id (1234567), DOI (10.5281/zenodo.1234567) or link (https://zenodo.org/records/1234567)',
});

export const ListParams = { record: RecordParam };

export const LoadParams = {
  record: RecordParam,
  file: PD.Text('', { description: 'The file to load, by its name in the record (zenodo.list names them)' }),
  coordinates: PD.Text('', {
    description: 'A coordinates file of the same record (xtc, dcd, …): loads a trajectory with file as its model',
  }),
};

const KindsByCategory: [string, ZenodoFileKind][] = [
  [TrajectoryFormatCategory, 'structure'],
  [VolumeFormatCategory, 'volume'],
  [TopologyFormatCategory, 'topology'],
  [CoordinatesFormatCategory, 'coordinates'],
];

/** What the plugin reads, by extension: the first format of each kind (as v5's Zenodo import). */
function formatsByExtension(
  plugin: PluginContext,
): Map<string, { kind: ZenodoFileKind; format: string; binary: boolean }> {
  const out = new Map<string, { kind: ZenodoFileKind; format: string; binary: boolean }>();
  for (const [category, kind] of KindsByCategory) {
    for (const { name, provider } of plugin.dataFormats.list) {
      if (provider.category !== category) continue;
      for (const e of provider.binaryExtensions ?? [])
        if (!out.has(e)) out.set(e, { kind, format: name, binary: true });
      for (const e of provider.stringExtensions ?? [])
        if (!out.has(e)) out.set(e, { kind, format: name, binary: false });
    }
  }
  return out;
}

const underZenodo = (url: string) => {
  try {
    const u = new URL(url);
    return u.protocol === 'https:' && (u.hostname === 'zenodo.org' || u.hostname.endsWith('.zenodo.org'));
  } catch {
    return false;
  }
};

function recordInfo(plugin: PluginContext, id: string, json: RecordJson): ZenodoRecordInfo {
  const formats = formatsByExtension(plugin);
  const files = (json.files ?? []).flatMap((f): ZenodoFileInfo[] => {
    const name = f.key ?? f.filename;
    if (!name) return [];
    // The record's own link when it is Zenodo's; otherwise the API's address for the file (what `addresses` says).
    const link = f.links?.content ?? f.links?.self;
    const url = link && underZenodo(link) ? link : zenodoFileAddress(id, name);
    const ext = getFileNameInfo(name).ext || (f.type ?? '').toLowerCase();
    const known = formats.get(ext);
    return [{ name, size: f.size ?? 0, url, ...(known ?? {}) }];
  });
  const doi = json.doi ?? json.metadata?.doi;
  return { id, title: json.metadata?.title ?? `Record ${id}`, ...(doi ? { doi } : {}), files };
}

const notAnId = (record: string) =>
  ApiResult.fail(
    'invalid',
    record.trim() ? `'${record}' is not a Zenodo record.` : 'Which Zenodo record?',
    'Give its id (1234567), its DOI (10.5281/zenodo.1234567) or its link.',
  );

const sizeOf = (n: number) =>
  n >= 1e9 ? `${(n / 1e9).toFixed(1)} GB` : n >= 1e6 ? `${(n / 1e6).toFixed(1)} MB` : `${Math.ceil(n / 1e3)} kB`;

export class ZenodoApi extends ApiArea {
  readonly id = 'zenodo';

  /** Records read in this session, by id: a load right after a list does not ask Zenodo again. */
  readonly #records = new Map<string, ZenodoRecordInfo>();

  readonly list = this.operation({
    ...ZenodoOperations.list,
    params: ListParams,
    run: async ({ record }, ctx): Promise<ApiResult<ZenodoRecordInfo>> => {
      const id = zenodoRecordId(record);
      if (!id) return notAnId(record);
      const r = await this.#record(ctx, id);
      if (ApiResult.isFailure(r)) return r;
      const readable = r.files.filter((f) => f.kind);
      const lines = r.files.map(
        (f) => `${f.name} (${sizeOf(f.size)}${f.kind ? `, ${f.kind}` : ', not read by the viewer'})`,
      );
      return ApiResult.ok(
        `Record ${id}: ${r.title}. ${r.files.length} file${r.files.length === 1 ? '' : 's'}, ${readable.length} the viewer reads${lines.length ? `:\n${lines.join('\n')}` : '.'}`,
        {},
        { value: r },
      );
    },
  });

  readonly load = this.operation({
    ...ZenodoOperations.load,
    params: LoadParams,
    run: async ({ record, file, coordinates }, ctx) => {
      const id = zenodoRecordId(record);
      if (!id) return notAnId(record);
      if (!file.trim()) return ApiResult.fail('invalid', 'Which file?', 'Give its name; zenodo.list names them.');
      const r = await this.#record(ctx, id);
      if (ApiResult.isFailure(r)) return r;
      const find = (name: string) => {
        const n = name.trim();
        return r.files.find((f) => f.name === n) ?? r.files.find((f) => f.name.toLowerCase() === n.toLowerCase());
      };
      const missing = (name: string) =>
        ApiResult.fail(
          'not-found',
          `Record ${id} has no file '${name}'.`,
          `Its files: ${r.files.map((f) => f.name).join(', ') || 'none'}.`,
        );
      const main = find(file);
      if (!main) return missing(file);
      if (!main.kind || !main.format) {
        return ApiResult.fail(
          'invalid',
          `${main.name} is not a format the viewer reads.`,
          'Pick a structure, volume, topology or coordinates file (zenodo.list says which).',
        );
      }
      let result: ApiResult;
      if (coordinates.trim()) {
        const coords = find(coordinates);
        if (!coords) return missing(coordinates);
        if (coords.kind !== 'coordinates' || (main.kind !== 'structure' && main.kind !== 'topology')) {
          return ApiResult.fail(
            'invalid',
            `${main.name} and ${coords.name} are not a model and its coordinates.`,
            'Give a structure or topology as file, and a coordinates file (xtc, dcd, …) as coordinates.',
          );
        }
        result = await runLoader(ctx, 'load.trajectory', {
          model: main.url,
          modelFormat: main.format,
          modelBinary: main.binary,
          coordinates: coords.url,
          coordinatesFormat: coords.format,
          coordinatesBinary: coords.binary,
        });
      } else if (main.kind === 'structure') {
        result = await runLoader(ctx, 'load.url', { url: main.url, format: main.format, binary: main.binary });
      } else if (main.kind === 'volume') {
        result = await runLoader(ctx, 'load.volume', { url: main.url, format: main.format, binary: main.binary });
      } else {
        return ApiResult.fail(
          'invalid',
          `${main.name} is ${main.kind === 'topology' ? 'a topology' : 'coordinates'} only.`,
          main.kind === 'topology'
            ? 'Load it with a coordinates file of the record (coordinates).'
            : 'Give the model or topology as file, and this as coordinates.',
        );
      }
      if (!result.ok || result.unchanged) return result;
      return ApiResult.ok(`Loaded ${main.name} from Zenodo record ${id}.`, result.changed, { value: result.value });
    },
  });

  async #record(ctx: ApiContext, id: string): Promise<ZenodoRecordInfo | ApiFailure> {
    const known = this.#records.get(id);
    if (known) return known;
    let json: RecordJson;
    try {
      json = (await ctx.plugin.runTask(ctx.plugin.fetch({ url: zenodoRecordAddress(id), type: 'json' }))) as RecordJson;
    } catch (e) {
      const message = messageOf(e);
      if (/\b(404|410)\b/.test(message)) {
        return ApiResult.fail('not-found', `There is no Zenodo record ${id}.`, 'Check the id or DOI.');
      }
      return ApiResult.fail(
        'failed',
        `Could not read Zenodo record ${id}: ${message}`,
        'Check the network connection.',
      );
    }
    if (!json || typeof json !== 'object' || !Array.isArray(json.files)) {
      return ApiResult.fail(
        'failed',
        `Zenodo record ${id} lists no files.`,
        'Its files may be restricted; open the record on zenodo.org.',
      );
    }
    const info = recordInfo(ctx.plugin, id, json);
    this.#records.set(id, info);
    return info;
  }
}

/** Runs a loader of the loaders feature as part of this run (its step); says so when the viewer has none. */
async function runLoader(ctx: ApiContext, id: string, values: Record<string, unknown>): Promise<ApiResult> {
  if (!ctx.api.describe(id)) {
    return ApiResult.fail('unavailable', `This viewer cannot load files from an address (no '${id}').`);
  }
  return ctx.api.run(id, values);
}

export const ZenodoExtension: ViewerExtension = {
  id: 'zenodo',
  title: 'Zenodo',
  api: { zenodo: (ctx) => new ZenodoApi(ctx) },
  commands: (api) => {
    const z = api.zenodo;
    if (!z) return [];
    return [
      operationCommand(api, z.list, {
        category: 'Load',
        icon: 'cloud',
        alias: ['zenodo', 'zenodo record'],
        action: 'List',
      }),
      operationCommand(api, z.load, {
        category: 'Load',
        icon: 'cloud',
        alias: ['zenodo', 'zenodo file'],
        action: 'Load',
      }),
      {
        id: ZenodoCommands.open,
        title: 'Show the Zenodo import',
        category: 'Load',
        icon: 'cloud',
        alias: ['zenodo import'],
        // Opens Home at the section: a view, not a change.
        async run() {
          const r = await api.view.open({ tab: 'home', section: ZenodoSection });
          if (ApiResult.isFailure(r)) throw new Error(r.message);
          return { message: 'Opened the Zenodo import.', ops: 'view.open({ tab: "home" })' };
        },
      },
    ];
  },
  skills: [ZenodoSkill],
  ui: async (ctx) => (await import('./prism-viewer-ui.js')).zenodoUI(ctx),
};
