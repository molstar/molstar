/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { type RuntimeContext, Task } from '@molstar/core/task';
import type { ReadonlyUint8Array } from '@molstar/core/util/type-helpers';
import { unzip } from '@molstar/core/util/zip/zip';
import { utf8Read } from '@molstar/core/util/utf8';
import { ReaderResult as Result } from '../result.js';
import type { SogArchive, SogFile, SogImage, SogImageBytes, SogImageDecoder, SogMeta } from './schema.js';
import { decodeWebpInBrowser } from './webp-decoder.js';

async function parseArchive(ctx: RuntimeContext, data: ReadonlyUint8Array): Promise<Result<SogArchive>> {
  // the archive's buffer as it is where the bytes are all of it (no copy)
  const whole = data.byteOffset === 0 && data.byteLength === data.buffer.byteLength;
  const buffer = whole ? (data.buffer as ArrayBuffer) : (data as Uint8Array).slice().buffer;
  const entries = await unzip(ctx, buffer);
  const file = (name: string) => {
    const e = entries[name];
    return e instanceof Uint8Array ? e : undefined;
  };
  const metaBytes = file('meta.json');
  if (!metaBytes) return Result.error('not a SOG file (no meta.json)');
  const meta = JSON.parse(utf8Read(metaBytes, 0, metaBytes.length)) as SogMeta;
  if (meta.version === undefined) return Result.error('SOG version 1 (no version in meta.json) is not supported');
  if (meta.version !== 2 && meta.version !== 3) return Result.error(`SOG version ${meta.version} is not supported`);

  const image = (name: string | undefined): SogImageBytes => {
    const bytes = name !== undefined ? file(name) : undefined;
    if (!bytes) throw new Error(`SOG image '${name}' missing`);
    return { name: name!, bytes };
  };
  return Result.success<SogArchive>({
    meta,
    version: meta.version,
    count: meta.count,
    meansLow: image(meta.means.files[0]),
    meansHigh: image(meta.means.files[1]),
    scales: image(meta.scales.files[0]),
    quats: image(meta.quats.files[0]),
    sh0: image(meta.sh0.files[0]),
    shN: meta.shN ? { centroids: image(meta.shN.files[0]), labels: image(meta.shN.files[1]) } : undefined,
  });
}

/**
 * Decodes the image `image` of `archive` by `decode`, checked to hold the splats (the palette apart: its size its own)
 */
export async function decodeSogImage(
  archive: SogArchive,
  image: SogImageBytes,
  decode: SogImageDecoder = decodeWebpInBrowser,
): Promise<SogImage> {
  const decoded = await decode(image.bytes, image.name);
  if (decoded.width * decoded.height < archive.count && image !== archive.shN?.centroids) {
    throw new Error(`SOG image '${image.name}' of ${decoded.width}x${decoded.height} for ${archive.count} splats`);
  }
  return decoded;
}

/** The images of `archive` decoded by `decode` */
export async function decodeSogArchive(
  archive: SogArchive,
  decode: SogImageDecoder = decodeWebpInBrowser,
): Promise<SogFile> {
  const image = (i: SogImageBytes) => decodeSogImage(archive, i, decode);
  const { shN } = archive;
  const [meansLow, meansHigh, scales, quats, sh0, centroids, labels] = await Promise.all([
    image(archive.meansLow),
    image(archive.meansHigh),
    image(archive.scales),
    image(archive.quats),
    image(archive.sh0),
    shN ? image(shN.centroids) : undefined,
    shN ? image(shN.labels) : undefined,
  ]);
  return {
    meta: archive.meta,
    version: archive.version,
    count: archive.count,
    meansLow,
    meansHigh,
    scales,
    quats,
    sh0,
    shN: centroids && labels ? { centroids, labels } : undefined,
  };
}

/**
 * Reads PlayCanvas's SOG archive (a zip of `meta.json` and lossless WebP images, versions 2 and 3): the meta and the
 * images as stored (views of the archive's entries)
 */
export function parseSogArchive(data: ReadonlyUint8Array) {
  return Task.create<Result<SogArchive>>('Parse SOG archive', async (ctx) => {
    try {
      return await parseArchive(ctx, data);
    } catch (e) {
      return Result.error(e instanceof Error ? e.message : String(e));
    }
  });
}

/**
 * Parses PlayCanvas's SOG (a zip of `meta.json` and lossless WebP images, versions 2 and 3): the images decoded by
 * `decode` (by default the browser's, `decodeWebpInBrowser`: none in Node)
 */
export function parseSog(data: ReadonlyUint8Array, decode: SogImageDecoder = decodeWebpInBrowser) {
  return Task.create<Result<SogFile>>('Parse SOG', async (ctx) => {
    try {
      const archive = await parseArchive(ctx, data);
      if (archive.isError) return archive;
      return Result.success(await decodeSogArchive(archive.result, decode));
    } catch (e) {
      return Result.error(e instanceof Error ? e.message : String(e));
    }
  });
}
