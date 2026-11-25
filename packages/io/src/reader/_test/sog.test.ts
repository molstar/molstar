/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { Task } from '@molstar/core/task';
import { zip } from '@molstar/core/util/zip/zip';
import { decodeSogImage, parseSog, parseSogArchive } from '../sog/parser.js';
import type { SogImage, SogImageDecoder } from '../sog/schema.js';

const Files = ['means_l.webp', 'means_u.webp', 'scales.webp', 'quats.webp', 'sh0.webp'];
const ShFiles = ['shN_centroids.webp', 'shN_labels.webp'];

/** a SOG of `count` splats with its meta (`version` undefined: none), each image's bytes its name's */
async function sog(meta: object) {
  const entries: { [k: string]: Uint8Array<ArrayBuffer> } = {
    'meta.json': new TextEncoder().encode(JSON.stringify(meta)),
  };
  for (const f of [...Files, ...ShFiles]) entries[f] = new TextEncoder().encode(f);
  // the archive as the plugin's bytes are: a `Uint8Array` (`zip` gives its buffer)
  return new Uint8Array(await Task.create('zip', (ctx) => zip(ctx, entries)).run());
}

function meta(version: number | undefined, count: number, sh = false) {
  return {
    ...(version !== undefined ? { version } : {}),
    count,
    means: { mins: [0, 0, 0], maxs: [1, 1, 1], files: Files.slice(0, 2) },
    scales: { codebook: [], files: [Files[2]] },
    quats: { files: [Files[3]] },
    sh0: { codebook: [], files: [Files[4]] },
    ...(sh ? { shN: { codebook: [], files: ShFiles } } : {}),
  };
}

/** a decoder of images of `width` x `height`, each filled with its name's length, the names decoded recorded */
function decoder(width: number, height: number, decoded: string[]): SogImageDecoder {
  return async (bytes, name) => {
    decoded.push(name);
    expect(new TextDecoder().decode(bytes)).toBe(name);
    const image: SogImage = { width, height, data: new Uint8Array(width * height * 4).fill(name.length) };
    return image;
  };
}

describe('sog', () => {
  it('reads meta.json and decodes the images it names, the palette where it has one', async () => {
    for (const version of [2, 3]) {
      const decoded: string[] = [];
      const parsed = await parseSog(await sog(meta(version, 6, true)), decoder(3, 2, decoded)).run();
      if (parsed.isError) throw new Error(parsed.message);
      const f = parsed.result;
      expect([f.version, f.count]).toEqual([version, 6]);
      expect(decoded.sort()).toEqual([...Files, ...ShFiles].sort());
      expect(f.meansLow.data[0]).toBe('means_l.webp'.length);
      expect(f.sh0.width * f.sh0.height).toBe(6);
      expect(f.shN?.labels.data[0]).toBe('shN_labels.webp'.length);
    }
    const plain = await parseSog(await sog(meta(3, 6)), decoder(3, 2, [])).run();
    expect(plain.isError ? undefined : plain.result.shN).toBeUndefined();
  });

  it('rejects version 1 (no version), other versions, images smaller than the splats, no meta.json', async () => {
    for (const [m, error] of [
      [meta(undefined, 4), /version 1/],
      [meta(4, 4), /version 4/],
      [meta(3, 7), /for 7 splats/],
    ] as const) {
      const parsed = await parseSog(await sog(m), decoder(3, 2, [])).run();
      expect(parsed.isError && parsed.message).toMatch(error);
    }
    const noMeta = new Uint8Array(await Task.create('zip', (ctx) => zip(ctx, { 'a.webp': new Uint8Array(4) })).run());
    expect((await parseSog(noMeta, decoder(1, 1, [])).run()).isError).toBe(true);
  });

  // the graphics upload the images as stored, decoding the means alone
  it("reads an archive of a Node Buffer at an offset in a larger one (a pool's view)", async () => {
    const bytes = await sog(meta(3, 6, true));
    const pool = new Uint8Array(bytes.length + 16).fill(0xee);
    pool.set(bytes, 3);
    const parsed = await parseSogArchive(Buffer.from(pool.buffer, 3, bytes.length)).run();
    if (parsed.isError) throw new Error(parsed.message);
    expect([parsed.result.version, parsed.result.count]).toEqual([3, 6]);
    expect(new TextDecoder().decode(parsed.result.quats.bytes)).toBe('quats.webp');
  });

  it('reads the archive without decoding: the images as stored, decoded one by one on request', async () => {
    const parsed = await parseSogArchive(await sog(meta(3, 6, true))).run();
    if (parsed.isError) throw new Error(parsed.message);
    const archive = parsed.result;
    expect([archive.version, archive.count]).toEqual([3, 6]);
    expect(new TextDecoder().decode(archive.quats.bytes)).toBe('quats.webp');
    expect(archive.shN?.centroids.name).toBe('shN_centroids.webp');
    const decoded: string[] = [];
    const means = await decodeSogImage(archive, archive.meansLow, decoder(3, 2, decoded));
    expect(decoded).toEqual(['means_l.webp']);
    expect(means.data[0]).toBe('means_l.webp'.length);
    // smaller than the splats: refused, the palette apart (its size its own)
    await expect(decodeSogImage(archive, archive.scales, decoder(2, 2, []))).rejects.toThrow(/for 6 splats/);
    await expect(decodeSogImage(archive, archive.shN!.centroids, decoder(1, 1, []))).resolves.toBeDefined();
  });
});
