/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Paul Pillot <paul.pillot@tandemai.com>
 */

import fs from 'node:fs';
import path from 'node:path';
import readline from 'node:readline';
import { Readable } from 'node:stream';
import type { ReadableStream } from 'node:stream/web';
import { pipeline } from 'node:stream/promises';
import { createGunzip } from 'node:zlib';
import { CIF } from '@molstar/io/reader/cif';
import type { CifBlock } from '@molstar/io/reader/cif/data-model';

const ComponentsUrl = 'https://files.wwpdb.org/pub/pdb/data/monomers/components.cif.gz';

/** The CCD, decompressed once into the cache, read one component at a time through a byte index (plan §3.1). */
export interface Ccd {
  /** component IDs, in the file's order */
  ids: string[];
  read(id: string): Promise<CifBlock>;
}

export async function openCcd(cacheDir: string): Promise<Ccd> {
  const file = path.join(cacheDir, 'components.cif');
  const indexFile = path.join(cacheDir, 'components.index.json');
  if (!fs.existsSync(file)) await download(file);
  if (!fs.existsSync(indexFile)) fs.writeFileSync(indexFile, JSON.stringify(await buildIndex(file)));
  const index: Record<string, [number, number]> = JSON.parse(fs.readFileSync(indexFile, 'utf8'));
  const fd = fs.openSync(file, 'r');
  return {
    ids: Object.keys(index),
    async read(id) {
      const range = index[id];
      if (!range) throw new Error(`${id} is not in the CCD`);
      const buffer = Buffer.alloc(range[1] - range[0]);
      fs.readSync(fd, buffer, 0, buffer.length, range[0]);
      const parsed = await CIF.parseText(buffer.toString('utf8')).run();
      if (parsed.isError) throw new Error(`${id}: ${parsed.message}`);
      return parsed.result.blocks[0];
    },
  };
}

async function download(file: string) {
  console.log(`downloading ${ComponentsUrl}`);
  const response = await fetch(ComponentsUrl);
  if (!response.ok || !response.body) throw new Error(`${ComponentsUrl}: ${response.status}`);
  fs.mkdirSync(path.dirname(file), { recursive: true });
  const partial = `${file}.partial`;
  await pipeline(Readable.fromWeb(response.body as ReadableStream), createGunzip(), fs.createWriteStream(partial));
  fs.renameSync(partial, file);
}

/** Byte range of each `data_` block. */
async function buildIndex(file: string) {
  const index: Record<string, [number, number]> = {};
  const lines = readline.createInterface({ input: fs.createReadStream(file), crlfDelay: Infinity });
  let offset = 0;
  let current: string | undefined;
  let start = 0;
  for await (const line of lines) {
    if (line.startsWith('data_')) {
      if (current !== undefined) index[current] = [start, offset];
      current = line.slice(5).trim();
      start = offset;
    }
    offset += Buffer.byteLength(line) + 1;
  }
  if (current !== undefined) index[current] = [start, offset];
  return index;
}
