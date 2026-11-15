/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Asset } from '@molstar/core/util/assets';
import { ParticlePackSource } from '@molstar/model/formats/particles/pack-stream';

// Where a particle pack is read by range (`ParticlePackSource`): a file by `Blob.slice`, a URL by HTTP range requests.
// A URL whose server ignores ranges (a whole body for a range) none: read whole, as other data.

/** a pack's source of `asset`: a file's (the file held), a URL's where its server answers ranges; else none */
export async function getParticlePackSource(asset: Asset): Promise<ParticlePackSource | undefined> {
  if (Asset.isUrl(asset)) {
    const source = await getUrlRangeSource(asset.url, asset.headers);
    // a few requests at a time, ranges side by side merged (a view's targets: thousands of small ranges)
    return source && ParticlePackSource.batched(source);
  }
  if (Asset.isFile(asset) && asset.file) return ParticlePackSource.batched(ParticlePackSource.fromBlob(asset.file));
  return undefined;
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
