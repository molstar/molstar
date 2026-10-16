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

/** An AlphaFold DB entry: `AF-P04637-F1` (UniProt accession, fragment), and the model file's version when it said. */
export interface AfdbEntry {
  /** `AF-P04637-F1`. */
  readonly id: string;
  readonly accession: string;
  readonly fragment: number;
  readonly version?: number;
}

const IdPattern = /AF-([A-Z0-9]+)-F(\d+)/i;

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

/** The versions the conventional address is tried with when the API does not say: the model's own, then the newest known. */
const FallbackVersions = [6, 4];

/** Why a PAE could not be got (`Error.message` is for the person). */
export class AfdbError extends Error {}

type Fetch = (
  url: string,
  init?: { method?: string; signal?: AbortSignal },
) => Promise<Pick<Response, 'ok' | 'status' | 'json'>>;

/**
 * Where the entry's PAE JSON is: `paeDocUrl` of the API's entry for this id and fragment, else the conventional file
 * address (probed, in the versions the model file or the API named, then the known ones). Throws `AfdbError`.
 */
export async function resolveAfdbPaeUrl(
  entry: AfdbEntry,
  fetcher: Fetch = (url, init) => fetch(url, init),
  signal?: AbortSignal,
): Promise<string> {
  let unreachable: string | undefined;
  let latest: number | undefined;
  try {
    const r = await fetcher(afdbApiUrl(entry.accession), { signal });
    if (r.ok) {
      const list = (await r.json()) as unknown;
      const found = (Array.isArray(list) ? list : []).find(
        (e): e is Record<string, unknown> =>
          !!e && typeof e === 'object' && String((e as Record<string, unknown>).entryId).toUpperCase() === entry.id,
      );
      const url = found?.paeDocUrl;
      if (typeof url === 'string' && isAfdbUrl(url)) return url;
      if (typeof found?.latestVersion === 'number') latest = found.latestVersion;
    }
  } catch (e) {
    if (signal?.aborted) throw e;
    unreachable = e instanceof Error ? e.message : String(e);
  }
  const versions = [...new Set([entry.version, latest, ...FallbackVersions].filter((v): v is number => !!v))];
  for (const v of versions) {
    const url = afdbPaeUrl(entry, v);
    try {
      if ((await fetcher(url, { method: 'HEAD', signal })).ok) return url;
    } catch (e) {
      if (signal?.aborted) throw e;
      unreachable ??= e instanceof Error ? e.message : String(e);
    }
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
