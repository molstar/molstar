/**
 * Copyright (c) 2024-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Ludovic Autin <ludovic.autin@gmail.com>
 */

import { Task, type RuntimeContext, chunkedSubtask } from '@molstar/core/task';
import {
  TextByteTokenizer,
  TextByteTokenBuilder as TokenBuilder,
  type TextByteTokens,
} from '../../common/text/text-byte-tokenizer.js';
import { ReaderResult as Result } from '../../result.js';
import { TextByteTokenColumnProvider as TokenColumn } from '../../common/text/column/text-byte-token.js';
import { Column } from '@molstar/core/data/db';
import type { LammpsFrame, LammpsTrajectoryFile } from '../schema.js';
import type { StringLike } from '@molstar/core/util/string-like';
import type { ReadonlyUint8Array } from '@molstar/core/util/type-helpers';

const { readLine } = TextByteTokenizer;

function State(tokenizer: TextByteTokenizer, runtimeCtx: RuntimeContext) {
  return {
    tokenizer,
    runtimeCtx,
  };
}
type State = ReturnType<typeof State>;

async function handleAtoms(state: State, count: number, parts: string[]): Promise<LammpsFrame> {
  const { tokenizer } = state;
  const columnIndexMap = Object.fromEntries(parts.map((colName, index) => [colName, index]));
  // declare column x, y, and z by check first caracter to 'x' or 'y' or 'z'
  // x,y,z = unscaled atom coordinates
  // xs,ys,zs = scaled atom coordinates this need the boundary box
  // xu,yu,zu = unwrapped atom coordinates
  // xsu,ysu,zsu = scaled unwrapped atom coordinates
  // ix,iy,iz = box image that the atom is in: read when given, with wrapped coordinates they give the unwrapped ones
  // how should we handle the different scenario ?
  const xCol = parts.findIndex((p) => p[0] === 'x');
  const yCol = parts.findIndex((p) => p[0] === 'y');
  const zCol = parts.findIndex((p) => p[0] === 'z');
  // retrieve the atom type colum for x only
  const atomMode = parts[xCol]; // x,xs,xu,xsu
  const atomId = TokenBuilder.create(tokenizer.data, count * 2);
  const moleculeType = TokenBuilder.create(tokenizer.data, count * 2);
  const atomType = TokenBuilder.create(tokenizer.data, count * 2);
  const x = TokenBuilder.create(tokenizer.data, count * 2);
  const y = TokenBuilder.create(tokenizer.data, count * 2);
  const z = TokenBuilder.create(tokenizer.data, count * 2);
  const hasImages = 'ix' in columnIndexMap && 'iy' in columnIndexMap && 'iz' in columnIndexMap;
  const ix = TokenBuilder.create(tokenizer.data, hasImages ? count * 2 : 0);
  const iy = TokenBuilder.create(tokenizer.data, hasImages ? count * 2 : 0);
  const iz = TokenBuilder.create(tokenizer.data, hasImages ? count * 2 : 0);

  const n = parts.length;
  // the tokens each column goes to, by index: a lookup per token, not a switch over the column map; the first match
  // wins, as in a switch (`x` is the first column starting with x, and so on)
  const targets: (TextByteTokens | undefined)[] = new Array(n).fill(undefined);
  const assign = (col: number | undefined, tokens: TextByteTokens) => {
    if (col !== undefined && col >= 0 && col < n && !targets[col]) targets[col] = tokens;
  };
  assign(columnIndexMap['id'], atomId);
  assign(columnIndexMap['mol'], moleculeType);
  assign(columnIndexMap['type'], atomType);
  assign(xCol, x);
  assign(yCol, y);
  assign(zCol, z);
  if (hasImages) {
    assign(columnIndexMap['ix'], ix);
    assign(columnIndexMap['iy'], iy);
    assign(columnIndexMap['iz'], iz);
  }

  const { data, length } = tokenizer;
  let linesAlreadyRead = 0;
  await chunkedSubtask(
    state.runtimeCtx,
    100000,
    void 0,
    (chunkSize) => {
      const linesToRead = Math.min(count - linesAlreadyRead, chunkSize);
      // the tokenizer's skipWhitespace, eatValue and eatLine inlined with the position a local: the hot loop of a dump
      let p = tokenizer.position;
      let lineEnd = tokenizer.tokenEnd;
      for (let i = 0; i < linesToRead; ++i) {
        for (let j = 0; j < n; ++j) {
          while (p < length) {
            const c = data[p];
            if (c === 32 || c === 9 || c === 10 || c === 13) ++p;
            else break;
          }
          const start = p;
          while (p < length) {
            const c = data[p];
            if (c <= 32 && (c === 32 || c === 10 || c === 13 || c === 9)) break;
            ++p;
          }
          const tokens = targets[j];
          if (tokens !== undefined) TokenBuilder.addUnchecked(tokens, start, p);
        }
        // ignore any extra columns
        lineEnd = length;
        while (p < length) {
          const c = data[p];
          if (c === 10) {
            lineEnd = p++;
            break;
          }
          if (c === 13) {
            lineEnd = p++;
            if (data[p] === 10) ++p;
            break;
          }
          ++p;
        }
      }
      tokenizer.position = p;
      tokenizer.tokenStart = p;
      tokenizer.tokenEnd = lineEnd;
      linesAlreadyRead += linesToRead;
      return linesToRead;
    },
    (ctx) => ctx.update({ message: 'Parsing...', current: tokenizer.position, max: length }),
  );

  return {
    count,
    atomMode: atomMode,
    atomId: TokenColumn(atomId)(Column.Schema.int),
    moleculeId: TokenColumn(moleculeType)(Column.Schema.int),
    atomType: TokenColumn(atomType)(Column.Schema.int),
    x: TokenColumn(x)(Column.Schema.float),
    y: TokenColumn(y)(Column.Schema.float),
    z: TokenColumn(z)(Column.Schema.float),
    ...(hasImages
      ? {
          image: {
            ix: TokenColumn(ix)(Column.Schema.int),
            iy: TokenColumn(iy)(Column.Schema.int),
            iz: TokenColumn(iz)(Column.Schema.int),
          },
        }
      : {}),
  };
}

/**
 * Possible Attributes from Lammps Dump
 * see https://docs.lammps.org/dump.html fro more details
 * possible attributes = id, mol, proc, procp1, type, element, mass,
 *                     x, y, z, xs, ys, zs, xu, yu, zu,
 *                     xsu, ysu, zsu, ix, iy, iz,
 *                     vx, vy, vz, fx, fy, fz,
 *                     q, mux, muy, muz, mu,
 *                     radius, diameter, omegax, omegay, omegaz,
 *                     angmomx, angmomy, angmomz, tqx, tqy, tqz,
 *                     c_ID, c_ID[I], f_ID, f_ID[I], v_name,
 *                     i_name, d_name, i2_name[I], d2_name[I]
 * ITEM: BOX BOUNDS xx yy zz
 * xlo xhi
 * ylo yhi
 * zlo zhi
 * or, for a triclinic box:
 * ITEM: BOX BOUNDS xy xz yz xx yy zz
 * xlo_bound xhi_bound xy
 * ylo_bound yhi_bound xz
 * zlo_bound zhi_bound yz
 */
async function parseInternal(data: ReadonlyUint8Array, ctx: RuntimeContext): Promise<Result<LammpsTrajectoryFile>> {
  const tokenizer = TextByteTokenizer(data);
  const state = State(tokenizer, ctx);
  const f: LammpsTrajectoryFile = {
    frames: [],
    times: [],
    bounds: [],
    timeOffset: 0.0,
    deltaTime: 0.0,
  };
  const frames = f.frames;
  let numAtoms = 0;
  let timestep = 0;
  while (tokenizer.tokenEnd < tokenizer.length) {
    const line = readLine(state.tokenizer).trim();
    if (line.includes('ITEM: TIMESTEP')) {
      timestep = parseInt(readLine(state.tokenizer).trim());
      f.times.push(timestep);
    } else if (line.includes('ITEM: NUMBER OF ATOMS')) {
      numAtoms = parseInt(readLine(state.tokenizer).trim());
    } else if (line.includes('ITEM: ATOMS')) {
      // this line provide also the style of the output and will give the order of the columns
      const parts = line.split(' ').slice(2);
      const frame: LammpsFrame = await handleAtoms(state, numAtoms, parts);
      frames.push(frame);
    } else if (line.includes('ITEM: BOX BOUNDS')) {
      // `ITEM: BOX BOUNDS pp pp pp`, or for a triclinic box `ITEM: BOX BOUNDS xy xz yz pp pp pp`: then each line has
      // its tilt factor third and the bounds are of the box's bounding box (the LAMMPS manual, `dump`)
      const tokens = line.split('ITEM: BOX BOUNDS ')[1].trim().split(/\s+/);
      const triclinic = tokens[0] === 'xy';
      const flags = tokens.slice(-3);
      const bound = () => readLine(state.tokenizer).trim().split(/\s+/).map(parseFloat);
      const [xlb, xhb, xy = 0] = bound();
      const [ylb, yhb, xz = 0] = bound();
      const [zlo, zhi, yz = 0] = bound();
      const tilt: [number, number, number] = triclinic ? [xy, xz, yz] : [0, 0, 0];
      const xlo = xlb - Math.min(0, tilt[0], tilt[1], tilt[0] + tilt[1]);
      const xhi = xhb - Math.max(0, tilt[0], tilt[1], tilt[0] + tilt[1]);
      const ylo = ylb - Math.min(0, tilt[2]);
      const yhi = yhb - Math.max(0, tilt[2]);
      f.bounds.push({
        lower: [xlo, ylo, zlo],
        length: [xhi - xlo, yhi - ylo, zhi - zlo],
        tilt,
        periodicity: [flags[0], flags[1], flags[2]],
      });
    }
  }
  if (f.times.length >= 1) {
    f.timeOffset = f.times[0];
  }
  if (f.times.length >= 2) {
    f.deltaTime = f.times[1] - f.times[0];
  }
  return Result.success(f);
}

/** Parses a LAMMPS dump from its text or its UTF-8 bytes (scanned as bytes: `TextByteTokenizer`). */
export function parseLammpsTrajectory(data: StringLike | ReadonlyUint8Array) {
  return Task.create<Result<LammpsTrajectoryFile>>('Parse Lammp Trajectory', async (ctx) => {
    return await parseInternal(TextByteTokenizer.bytesOf(data), ctx);
  });
}
