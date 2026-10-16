/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { ApiResult, ViewerApi } from '../../api/area.js';

// The load parameters of the URL (app.md §4, the table; v5's `apps/viewer` names kept): read once, when the viewer is
// ready, and run as the `load.*` operations with the `url` origin, in v5's order. DOM-free: the app passes
// `location.search`.

/** One load the URL asks for: an operation and its values. */
export interface UrlLoad {
  readonly operation: string;
  readonly values: Record<string, unknown>;
}

/**
 * How the structure loads build what they load (Prism's own; v5 has them as `loadAllModelsOrAssemblyFromUrl`, not in
 * its URL): `all-models=1` loads every model of `structure-url` and `url`; `assembly=<id>` builds that assembly for
 * those and for the `pdb`, `pdb-ihm` and `model-archive` entries. They change loads, so they are not loads themselves:
 * the app knows them (`KnownUrlParameters`), and they alone do not stop the light restore.
 */
export const LoadUrlOptions: readonly string[] = ['all-models', 'assembly'];

/** The load parameters (the ones this module reads), for the app's "unknown parameter" log. */
export const LoadUrlParameters: readonly string[] = [
  'snapshot-id',
  'snapshot-url',
  'snapshot-url-type',
  'structure-url',
  'structure-url-format',
  'structure-url-is-binary',
  'mvs-url',
  'mvs-data',
  'mvs-format',
  'pdb',
  'pdb-ihm',
  'pdb-dev',
  'emdb',
  'afdb',
  'model-archive',
  'url',
  'url-format',
  'url-is-binary',
  'svs-url',
];

/** The parameters of a query string: names lowercased, the first of a name wins; values decoded as v5 does (`+` stays). */
export function parseSearch(search: string): Map<string, string> {
  const out = new Map<string, string>();
  for (const part of search.replace(/^\?/, '').split('&')) {
    if (!part) continue;
    const eq = part.indexOf('=');
    const name = decode(eq < 0 ? part : part.slice(0, eq)).toLowerCase();
    if (!out.has(name)) out.set(name, eq < 0 ? '' : decode(part.slice(eq + 1)));
  }
  return out;
}

function decode(s: string): string {
  try {
    return decodeURIComponent(s);
  } catch {
    return s;
  }
}

/** The loads a query string asks for, in order. Empty values are absent, as in v5. */
export function loadsFromUrl(search: string): UrlLoad[] {
  const q = parseSearch(search);
  const get = (name: string) => q.get(name)?.trim() || undefined;
  const flag = (name: string) => ['1', 'true'].includes(get(name)?.toLowerCase() ?? '');
  // A binary flag left out is told from the format and the address (operations.md §1, "Inputs").
  const binary = (name: string) => (get(name) === undefined ? {} : { binary: flag(name) });
  const format = (name: string) => (get(name) ? { format: get(name)!.toLowerCase() } : {});
  const loads: UrlLoad[] = [];
  const assembly = get('assembly');
  const structureOptions = {
    ...(flag('all-models') ? { models: 'all' } : {}),
    ...(assembly ? { assembly } : {}),
  };

  const snapshotId = get('snapshot-id');
  if (snapshotId) loads.push({ operation: 'load.snapshot', values: { id: snapshotId } });

  const snapshotUrl = get('snapshot-url');
  if (snapshotUrl) {
    // `load.snapshot` is `session.open`. Without `snapshot-url-type` the operation follows the address's extension (.molx
    // or .molj), else the content.
    const type = get('snapshot-url-type')?.toLowerCase();
    loads.push({ operation: 'load.snapshot', values: { url: snapshotUrl, ...(type ? { type } : {}) } });
  }

  const structureUrl = get('structure-url');
  if (structureUrl) {
    // v5's default format is mmCIF (text, or binary at a .bcif address).
    loads.push({
      operation: 'load.url',
      values: {
        url: structureUrl,
        format: get('structure-url-format')?.toLowerCase() ?? 'mmcif',
        ...binary('structure-url-is-binary'),
        ...structureOptions,
      },
    });
  }

  // Both is an error, not a guess: the operation says so (v5 logs it and ignores both). Without `mvs-format` the kind is
  // read from the address, else the content (base64 data is an archive's).
  const mvsUrl = get('mvs-url');
  const mvsData = get('mvs-data');
  if (mvsUrl || mvsData) {
    loads.push({
      operation: 'load.mvs',
      values: {
        ...(mvsUrl ? { url: mvsUrl } : {}),
        ...(mvsData ? { data: mvsData } : {}),
        ...format('mvs-format'),
      },
    });
  }

  for (const [name, source] of [
    ['pdb', 'pdb'],
    ['pdb-ihm', 'pdb-ihm'],
    ['pdb-dev', 'pdb-ihm'], // deprecated, the same
    ['emdb', 'emdb'],
    ['afdb', 'afdb'],
    ['model-archive', 'model-archive'],
  ] as const) {
    const id = get(name);
    const structure = source === 'pdb' || source === 'pdb-ihm' || source === 'model-archive';
    if (id) {
      loads.push({ operation: 'load.entry', values: { source, id, ...(structure && assembly ? { assembly } : {}) } });
    }
  }

  const url = get('url');
  if (url) {
    // v5 needs `url-format`; without it the format is read from the address (`auto`, the operation's default).
    loads.push({
      operation: 'load.url',
      values: { url, ...format('url-format'), ...binary('url-is-binary'), ...structureOptions },
    });
  }

  // Prism's own: a SeqViewSpec in the sequence, last, so the structure it names may be one loaded above.
  const svsUrl = get('svs-url');
  if (svsUrl) loads.push({ operation: 'load.svs', values: { url: svsUrl } });
  return loads;
}

/** What one URL load did. */
export interface UrlLoadResult extends UrlLoad {
  readonly result: ApiResult;
}

/**
 * Runs the loads of a query string one after the other, as the `url` origin (the history and the run log say who asked).
 * A failed one is in the results, and does not stop the rest.
 */
export async function runUrlLoads(api: ViewerApi, search: string): Promise<UrlLoadResult[]> {
  const asUrl = api.as({ kind: 'url' });
  const results: UrlLoadResult[] = [];
  for (const load of loadsFromUrl(search)) {
    results.push({ ...load, result: await asUrl.run(load.operation, load.values) });
  }
  return results;
}
