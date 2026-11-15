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

async function run(args: Args) {
  const data = new Uint8Array(fs.readFileSync(args.src));
  const start = performance.now();
  if (isParticlePack(data) && /\.mspk$/i.test(args.out)) {
    fs.writeFileSync(args.out, packToPack(data, { tiles: args.tiles, levels: args.levels }));
  } else if (isParticlePack(data)) {
    const binary = /\.bcif$/i.test(args.out);
    const { data: out, unmatchedRadii } = packToMmcif(data, binary);
    fs.writeFileSync(args.out, out);
    if (unmatchedRadii > 0) console.warn(`${unmatchedRadii} sites of a radius no element has, written as carbon`);
  } else {
    const cif = await parseMmcif(data, /\.bcif$/i.test(args.src));
    fs.writeFileSync(
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
