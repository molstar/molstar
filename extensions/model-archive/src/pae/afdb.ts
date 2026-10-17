/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { readPae } from '@molstar/prism-viewer/sequence/molstar';

// AlphaFold DB's PAE (parity.md §6, "PAE plot"): which entry a structure is, where its PAE JSON is, and the JSON as the
// card's matrix. No plugin here, so it is tested on plain values. The PAE file is AlphaFold DB's, apart from the model.

/** An AlphaFold DB entry: `AF-P04637-F1`, `AF-P04637-2-F1` for an isoform (UniProt accession, fragment), and the model file's version when it said. */
export interface AfdbEntry {
  /** `AF-P04637-F1`. */
  readonly id: string;
  readonly accession: string;
  readonly fragment: number;
  readonly version?: number;
}

// A UniProt accession, with an isoform suffix (`P04637-2`) when the entry is an isoform's.
const IdPattern = /AF-([A-Z0-9]+(?:-\d+)?)-F(\d+)/i;

/**
 * The entry a structure is, from the words it carries (its entry id, label, file name, the address it was loaded from,
 * in order of trust); undefined when none names one. A model file's `_v6` suffix is its version.
 */
export function afdbEntryOf(texts: readonly (string | undefined)[]): AfdbEntry | undefined {
  for (const t of texts) {
    const m = t && IdPattern.exec(t);
    if (!m) continue;
    const accession = m[1].toUpperCase();
    const fragment = Number(m[2]);
    const v = /_v(\d+)(?:[^0-9]|$)/i.exec(t.slice(m.index + m[0].length));
    return {
      id: `AF-${accession}-F${fragment}`,
      accession,
      fragment,
      ...(v ? { version: Number(v[1]) } : {}),
    };
  }
  return undefined;
}

export const AfdbHost = 'alphafold.ebi.ac.uk';

/** The API's listing of an accession's entries. */
export const afdbApiUrl = (accession: string) => `https://${AfdbHost}/api/prediction/${encodeURIComponent(accession)}`;

/** The conventional file address of an entry's PAE. */
export const afdbPaeUrl = (entry: AfdbEntry, version: number) =>
  `https://${AfdbHost}/files/${entry.id}-predicted_aligned_error_v${version}.json`;

/** The versions the conventional address is tried with when the API does not say, after the model's own: the newest known. */
const FallbackVersions = [6, 4];

/** Why a PAE could not be got (`Error.message` is for the person). */
export class AfdbError extends Error {}

type Fetch = (
  url: string,
  init?: { method?: string; signal?: AbortSignal },
) => Promise<Pick<Response, 'ok' | 'status' | 'json'>>;

/** Where a PAE is, and what to tell when it is not the model's own version's. */
export interface AfdbPaeSource {
  readonly url: string;
  /** Set when the file is of another version than the model file's. */
  readonly note?: string;
}

const versionOfUrl = (url: string) => {
  const v = /_v(\d+)\.json(?:$|[?#])/i.exec(url);
  return v ? Number(v[1]) : undefined;
};

/**
 * Where the entry's PAE JSON is. When the model file names its version, that version's file: the API's `paeDocUrl`
 * (the latest version's) when that is it, else the conventional address, probed. Otherwise `paeDocUrl` of the API's
 * entry for this id and fragment, else the conventional address in the versions the API or the known ones name. A file
 * of another version than the model's comes with a note saying so. Throws `AfdbError`.
 */
export async function resolveAfdbPaeUrl(
  entry: AfdbEntry,
  fetcher: Fetch = (url, init) => fetch(url, init),
  signal?: AbortSignal,
): Promise<AfdbPaeSource> {
  let unreachable: string | undefined;
  let latest: number | undefined;
  let docUrl: string | undefined;
  const source = (url: string): AfdbPaeSource => {
    const v = versionOfUrl(url);
    return entry.version && v && v !== entry.version
      ? {
          url,
          note: `The model file is version ${entry.version}, but this PAE is version ${v}: its errors may belong to another prediction.`,
        }
      : { url };
  };
  const probe = async (url: string): Promise<boolean> => {
    try {
      return (await fetcher(url, { method: 'HEAD', signal })).ok;
    } catch (e) {
      if (signal?.aborted) throw e;
      unreachable ??= e instanceof Error ? e.message : String(e);
      return false;
    }
  };
  try {
    const r = await fetcher(afdbApiUrl(entry.accession), { signal });
    if (r.ok) {
      const list = (await r.json()) as unknown;
      const found = (Array.isArray(list) ? list : []).find(
        (e): e is Record<string, unknown> =>
          !!e && typeof e === 'object' && String((e as Record<string, unknown>).entryId).toUpperCase() === entry.id,
      );
      const url = found?.paeDocUrl;
      if (typeof url === 'string' && isAfdbUrl(url)) docUrl = url;
      if (typeof found?.latestVersion === 'number') latest = found.latestVersion;
    }
  } catch (e) {
    if (signal?.aborted) throw e;
    unreachable = e instanceof Error ? e.message : String(e);
  }
  if (docUrl && (!entry.version || versionOfUrl(docUrl) === entry.version)) return { url: docUrl };
  // The model names a version the API's address is not: its own file first, the API's (the latest) only without it.
  if (entry.version) {
    const own = afdbPaeUrl(entry, entry.version);
    if (await probe(own)) return { url: own };
  }
  if (docUrl) return source(docUrl);
  for (const v of new Set([latest, ...FallbackVersions])) {
    if (!v || v === entry.version) continue;
    const url = afdbPaeUrl(entry, v);
    if (await probe(url)) return source(url);
  }
  throw new AfdbError(
    unreachable ? `Could not reach AlphaFold DB (${unreachable}).` : `AlphaFold DB has no PAE for ${entry.id}.`,
  );
}

/** An address on AlphaFold DB (https): what the API may point the fetch at. */
function isAfdbUrl(url: string): boolean {
  try {
    const u = new URL(url);
    return u.protocol === 'https:' && (u.hostname === AfdbHost || u.hostname.endsWith('.alphafold.ebi.ac.uk'));
  } catch {
    return false;
  }
}

/**
 * AlphaFold DB's PAE JSON (as parsed) as `size` × `size` row-major numbers (row i: residue i aligned on, column j:
 * residue j scored); a message when it is not one. The reader is Seq*'s `afdb_pae` (it keeps the matrix by scored
 * residue, so it is transposed back for the card's).
 */
export function paeValuesFromAfdb(json: unknown): { size: number; values: Float32Array; max: number } | string {
  const m = readPae(json);
  if (typeof m === 'string') return m;
  const { size } = m;
  const values = new Float32Array(size * size);
  let max = 0;
  for (let i = 0; i < size; i++) {
    for (let j = 0; j < size; j++) {
      const v = m.values[j * size + i];
      values[i * size + j] = v;
      if (v > max) max = v;
    }
  }
  return { size, values, max };
}
