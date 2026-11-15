/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Asset, type AssetManager } from '@molstar/core/util/assets';
import { ParticlePackSource } from '@molstar/model/formats/particles/pack-stream';

// Where a particle pack is read by range (`ParticlePackSource`): a file by `Blob.slice`, a URL by HTTP range requests.
// A URL whose server ignores ranges (a whole body for a range) none: read whole, as other data.

/**
 * a pack's source of `asset`: a file's (the file held: the asset's, or `assets`' where it holds it), a URL's where its
 * server answers ranges; else none
 */
export async function getParticlePackSource(
  asset: Asset,
  assets?: AssetManager,
): Promise<ParticlePackSource | undefined> {
  if (Asset.isUrl(asset)) {
    const source = await getUrlRangeSource(asset.url, asset.headers);
    // a few requests at a time, ranges side by side merged (a view's targets: thousands of small ranges)
    return source && ParticlePackSource.batched(source);
  }
  const file = Asset.isFile(asset) ? getPackFile(asset, assets) : undefined;
  return file && ParticlePackSource.batched(ParticlePackSource.fromBlob(file));
}

/**
 * The file of a pack's file asset: the one `assets` holds (a session's, its asset read back; the asset's own file
 * not kept by the state), else the asset's own
 */
export function getPackFile(asset: Asset.File, assets?: AssetManager): File | undefined {
  const held = assets?.get(asset)?.file;
  if (held) return held;
  return asset.file instanceof File ? asset.file : undefined;
}

/**
 * `asset`'s file held by `assets` (not read: a pack read by range, never whole), so that a session saved keeps it;
 * whether it was held here (none before it): its holder deletes it as it goes
 */
export function holdParticlePackFile(asset: Asset.File, assets: AssetManager): boolean {
  if (assets.has(asset)) return false;
  const file = getPackFile(asset);
  if (!file) return false;
  assets.set(asset, file);
  return true;
}

/** the total of a `Content-Range` (`bytes 0-63/731477500`), none where unknown (`*`) */
function getRangeTotal(header: string | null) {
  const total = header ? /\/(\d+)\s*$/.exec(header)?.[1] : undefined;
  return total ? Number(total) : undefined;
}

/**
 * The source of `url` by range requests: its first bytes asked as a range, a partial answer (206, its total) a
 * source, else none (the answer's body let go: a server ignoring ranges sends it whole)
 */
export async function getUrlRangeSource(
  url: string,
  headers?: Record<string, string>,
): Promise<ParticlePackSource | undefined> {
  const probe = new AbortController();
  let response: Response;
  try {
    response = await fetch(url, { headers: { ...headers, Range: 'bytes=0-63' }, signal: probe.signal });
  } catch {
    return undefined;
  }
  const total = response.status === 206 ? getRangeTotal(response.headers.get('content-range')) : undefined;
  if (total === undefined) {
    probe.abort();
    return undefined;
  }
  await response.arrayBuffer();
  return {
    byteLength: total,
    read: async (offset, length) => {
      // none asked (an empty block): no request, `bytes=X-(X-1)` a range a server refuses (416)
      if (length <= 0) return new Uint8Array(0);
      const r = await fetch(url, { headers: { ...headers, Range: `bytes=${offset}-${offset + length - 1}` } });
      if (r.status !== 206) throw new Error(`particle pack range of '${url}' not answered (${r.status})`);
      const bytes = new Uint8Array(await r.arrayBuffer());
      if (bytes.byteLength !== length) throw new Error(`particle pack range of '${url}' answered short`);
      return bytes;
    },
  };
}
