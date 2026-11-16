#!/usr/bin/env node
/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as argparse from 'argparse';
import * as fs from 'fs';
import { isParticlePack } from '@molstar/model/formats/particles/pack';
import { mmcifToPack, packToMmcif, packToPack, parseMmcif } from '@molstar/particle-pack-cli/convert';

interface Args {
  src: string;
  out: string;
  assembly?: string;
  variant?: 'auto' | 'cellpack' | 'petworld' | 'standard';
  tiles?: number;
  levels?: boolean;
}

const parser = new argparse.ArgumentParser({
  add_help: true,
  description:
    'Convert an mmCIF particle assembly (cellPACK, PetWorld; .cif or .bcif) to a particle pack (.mspk), a pack back to mmCIF (the output extension: .cif or .bcif), or a pack to a pack (.mspk: the current version, tile-ordered with --tiles)',
});
parser.add_argument('src', { help: 'Source file path' });
parser.add_argument('out', { help: 'Output file path' });
parser.add_argument('-a', '--assembly', { help: 'The assembly id (the first by default)', required: false });
parser.add_argument('--variant', {
  help: 'The mmCIF layout (auto by default)',
  choices: ['auto', 'cellpack', 'petworld', 'standard'],
  required: false,
});
parser.add_argument('--tiles', {
  help: 'Tile-ordered for a tile set of this many tiles an axis (a power of two, e.g. 16; or coarser, 8, 4): made without a pass over the particles',
  type: 'int',
  required: false,
});

parser.add_argument('--levels', {
  help: 'With --tiles, coarse levels: far tiles drawn as samples of their particles, read as a few ranges (about a seventh more particles written)',
  action: 'store_true',
  required: false,
});

/** A read or a write at most this long: one of 2 GiB or more fails (`readFileSync` and `writeFileSync` too) */
const Chunk = 1 << 28;

/** `path`'s bytes, read in chunks: a pack of the yeast cell's grids is several GB */
function readFile(path: string): Uint8Array {
  const fd = fs.openSync(path, 'r');
  try {
    const out = new Uint8Array(fs.fstatSync(fd).size);
    for (let o = 0; o < out.byteLength; ) {
      const n = fs.readSync(fd, out, o, Math.min(Chunk, out.byteLength - o), o);
      if (n <= 0) throw new Error(`'${path}' read short`);
      o += n;
    }
    return out;
  } finally {
    fs.closeSync(fd);
  }
}

/** `data` written to `path` in chunks; a failed write leaves no file (not one of 0 bytes) */
function writeFile(path: string, data: Uint8Array | string) {
  const bytes = typeof data === 'string' ? Buffer.from(data) : data;
  const fd = fs.openSync(path, 'w');
  try {
    for (let o = 0; o < bytes.byteLength; ) o += fs.writeSync(fd, bytes, o, Math.min(Chunk, bytes.byteLength - o));
    fs.closeSync(fd);
  } catch (e) {
    fs.closeSync(fd);
    fs.rmSync(path, { force: true });
    throw e;
  }
}

async function run(args: Args) {
  const data = readFile(args.src);
  const start = performance.now();
  if (isParticlePack(data) && /\.mspk$/i.test(args.out)) {
    writeFile(args.out, packToPack(data, { tiles: args.tiles, levels: args.levels }));
  } else if (isParticlePack(data)) {
    const binary = /\.bcif$/i.test(args.out);
    const { data: out, unmatchedRadii } = packToMmcif(data, binary);
    writeFile(args.out, out);
    if (unmatchedRadii > 0) console.warn(`${unmatchedRadii} sites of a radius no element has, written as carbon`);
  } else {
    const cif = await parseMmcif(data, /\.bcif$/i.test(args.src));
    writeFile(
      args.out,

      await mmcifToPack(cif, {
        assemblyId: args.assembly,
        variant: args.variant,
        tiles: args.tiles,
        levels: args.levels,
      }),
    );
  }
  const size = fs.statSync(args.out).size;
  console.log(`${args.out}: ${(size / 1e6).toFixed(1)} MB in ${((performance.now() - start) / 1000).toFixed(1)} s`);
}

run(parser.parse_args() as Args).catch((e) => {
  console.error(e instanceof Error ? e.message : e);
  process.exit(1);
});
