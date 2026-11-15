/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Asset, type AssetManager } from '@molstar/core/util/assets';
import {
  ParticlePackHeaderBytes,
  parseParticlePackDocument,
  readParticlePackHeader,
} from '@molstar/model/formats/particles/pack';
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
  /** told why a URL is not read by range (read whole then) */
  warn?: (message: string) => void,
): Promise<ParticlePackSource | undefined> {
  if (Asset.isUrl(asset)) {
    const source = await getUrlRangeSource(asset.url, asset.headers, warn);
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

/** the size of `url` by a HEAD request (`Content-Length`, CORS-safelisted: read on another origin too), none unknown */
async function getHeadLength(url: string, headers?: Record<string, string>) {
  try {
    const head = await fetch(url, { method: 'HEAD', headers });
    const length = head.ok ? Number(head.headers.get('content-length')) : Number.NaN;
    return Number.isFinite(length) && length > 0 ? length : undefined;
  } catch {
    return undefined;
  }
}

/**
 * The size of the pack starting with `head` (its header's bytes) from the pack itself: its document's end or its last
 * block's, whichever is further; the document's bytes kept (the reader's first read of them). None where `head` is no
 * pack's or its document unread
 */
async function getPackLength(
  head: Uint8Array,
  read: (offset: number, length: number) => Promise<Uint8Array>,
): Promise<{ byteLength: number; document: { offset: number; bytes: Uint8Array } } | undefined> {
  try {
    const { jsonOffset, jsonLength } = readParticlePackHeader(head);
    const bytes = await read(jsonOffset, jsonLength);
    let byteLength = jsonOffset + jsonLength;
    // its blocks checked by the reader (within the size given there: their own ends here)
    for (const b of parseParticlePackDocument(bytes, Infinity).blocks) {
      byteLength = Math.max(byteLength, b.byteOffset + b.byteLength);
    }
    return { byteLength, document: { offset: jsonOffset, bytes } };
  } catch {
    return undefined;
  }
}

/**
 * The source of `url` by range requests: its first bytes asked as a range, a partial answer (206) a source of its
 * total, else none (the answer's body let go: a server ignoring ranges sends it whole). The total the answer's
 * `Content-Range`; on another origin that header is read only where the server exposes it
 * (`Access-Control-Expose-Headers: Content-Range`), else by a HEAD request's `Content-Length`, else from the pack
 * itself (its header and document: a presigned URL, signed for GET alone, its bucket exposing no headers). Where none
 * gives it (or the server ignores ranges) none, `warn` told why: the pack then read whole
 */
export async function getUrlRangeSource(
  url: string,
  headers?: Record<string, string>,
  warn?: (message: string) => void,
): Promise<ParticlePackSource | undefined> {
  const probe = new AbortController();
  let response: Response;
  try {
    response = await fetch(url, {
      headers: { ...headers, Range: `bytes=0-${ParticlePackHeaderBytes - 1}` },
      signal: probe.signal,
    });
  } catch {
    return undefined;
  }
  if (response.status !== 206) {
    probe.abort();
    if (response.ok) warn?.(`'${url}' answers no range requests (status ${response.status}): read whole`);
    return undefined;
  }

  const readRange = async (offset: number, length: number, signal?: AbortSignal) => {
    // none asked (an empty block): no request, `bytes=X-(X-1)` a range a server refuses (416)
    if (length <= 0) return new Uint8Array(0);
    // aborted: the request dropped (a tile no longer wanted)
    const r = await fetch(url, {
      headers: { ...headers, Range: `bytes=${offset}-${offset + length - 1}` },
      signal,
    });
    if (r.status !== 206) throw new Error(`particle pack range of '${url}' not answered (${r.status})`);
    const bytes = new Uint8Array(await r.arrayBuffer());
    if (bytes.byteLength !== length) throw new Error(`particle pack range of '${url}' answered short`);
    return bytes;
  };

  const head = new Uint8Array(await response.arrayBuffer());
  let total = getRangeTotal(response.headers.get('content-range')) ?? (await getHeadLength(url, headers));
  /** the document's bytes where read for the size: the reader's first read of them answered with them */
  let kept: { offset: number; bytes: Uint8Array } | undefined;
  if (total === undefined) {
    const known = await getPackLength(head, readRange);
    if (!known) {
      warn?.(
        `'${url}' answers range requests but its size is unknown: its server exposes no 'Content-Range' ` +
          "('Access-Control-Expose-Headers: Content-Range' on another origin) nor a HEAD 'Content-Length', and its " +
          'start is no particle pack: read whole',
      );
      return undefined;
    }
    total = known.byteLength;
    kept = known.document;
  }

  return {
    byteLength: total,
    read: async (offset, length, signal) => {
      if (kept && offset === kept.offset && length === kept.bytes.byteLength) {
        const { bytes } = kept;
        kept = undefined;
        return bytes;
      }
      return readRange(offset, length, signal);
    },
  };
}
