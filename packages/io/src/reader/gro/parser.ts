/**
 * Copyright (c) 2017-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Column } from '@molstar/core/data/db';
import { TextByteTokenizer } from '../common/text/text-byte-tokenizer.js';
import { TextByteFixedColumnProvider as FixedColumn } from '../common/text/column/text-byte-fixed.js';
import type * as Schema from './schema.js';
import { ReaderResult as Result } from '../result.js';
import { Task, type RuntimeContext } from '@molstar/core/task';
import type { StringLike } from '@molstar/core/util/string-like';
import { parseIntBytes } from '@molstar/core/util/number-parser';
import type { ReadonlyUint8Array } from '@molstar/core/util/type-helpers';

interface State {
  tokenizer: TextByteTokenizer;
  header: Schema.GroHeader;
  numberOfAtoms: number;
  runtimeCtx: RuntimeContext;
}

/** A malformed file: a `ReaderResult` error, not a failed task */
class GroError extends Error {}

function createEmptyHeader(): Schema.GroHeader {
  return {
    title: '',
    timeInPs: 0,
    hasVelocities: false,
    precision: { position: 0, velocity: 0 },
    box: [0, 0, 0],
    boxVectors: [0, 0, 0, 0, 0, 0, 0, 0, 0],
  };
}

function State(tokenizer: TextByteTokenizer, runtimeCtx: RuntimeContext): State {
  return {
    tokenizer,
    header: createEmptyHeader(),
    numberOfAtoms: 0,
    runtimeCtx,
  };
}

/**
 * title string (free format string, optional time in ps after 't=')
 */
function handleTitleString(state: State) {
  const { tokenizer, header } = state;
  let line = TextByteTokenizer.readLine(tokenizer);

  // skip potential empty lines...
  if (line.trim().length === 0) {
    line = TextByteTokenizer.readLine(tokenizer);
  }

  const timeOffset = line.lastIndexOf('t=');
  if (timeOffset >= 0) {
    header.timeInPs = parseFloat(line.substring(timeOffset + 2));
    header.title = line.substring(0, timeOffset).trim();
    if (header.title && header.title[header.title.length - 1] === ',') {
      header.title = header.title.substring(0, header.title.length - 1);
    }
  } else {
    header.title = line;
  }
}

/**
 * number of atoms (free format integer)
 */
function handleNumberOfAtoms(state: State) {
  const { tokenizer } = state;
  TextByteTokenizer.markLine(tokenizer);
  // `parseInt` of the line, from its bytes: leading spaces skipped, NaN without a digit
  const { data, tokenEnd } = tokenizer;
  let s = tokenizer.tokenStart;
  while (s < tokenEnd && (data[s] === 32 || data[s] === 9)) s++;
  const d = data[s] === 45 || data[s] === 43 ? s + 1 : s;
  state.numberOfAtoms = d < tokenEnd && data[d] >= 48 && data[d] <= 57 ? parseIntBytes(data, s, tokenEnd) : NaN;
}

/**
 * This format is fixed, ie. all columns are in a fixed position.
 * Optionally (for now only yet with trjconv) you can write gro files
 * with any number of decimal places, the format will then be n+5
 * positions with n decimal places (n+1 for velocities) in stead
 * of 8 with 3 (with 4 for velocities). Upon reading, the precision
 * will be inferred from the distance between the decimal points
 * (which will be n+5). Columns contain the following information
 * (from left to right):
 *     residue number (5 positions, integer)
 *     residue name (5 characters)
 *     atom name (5 characters)
 *     atom number (5 positions, integer)
 *     position (in nm, x y z in 3 columns, each 8 positions with 3 decimal places)
 *     velocity (in nm/ps (or km/s), x y z in 3 columns, each 8 positions with 4 decimal places)
 */
async function handleAtoms(state: State): Promise<Schema.GroAtoms> {
  const { tokenizer, numberOfAtoms } = state;
  const lines = await TextByteTokenizer.readLinesAsync(tokenizer, numberOfAtoms, state.runtimeCtx, 100000);

  // the decimal places of the numbers of the first atom line after its position columns (byte 20): `/\.\d+/g`
  const precisions = decimalPlaces(tokenizer.data, Math.min(lines.indices[0] + 20, lines.indices[1]), lines.indices[1]);
  if (numberOfAtoms > 0 && precisions.length < 3) throw new GroError('No decimal coordinates in the first atom line.');
  const hasVelocities = precisions.length === 6;

  state.header.hasVelocities = hasVelocities;
  state.header.precision.position = precisions[0];
  state.header.precision.velocity = hasVelocities ? precisions[3] : 0;

  const pO = 20;
  const pW = state.header.precision.position + 5;
  const vO = pO + 3 * pW;
  const vW = state.header.precision.velocity + 4;

  const col = FixedColumn(lines);
  const undef = Column.Undefined(state.numberOfAtoms, Column.Schema.float);

  const ret = {
    count: state.numberOfAtoms,
    residueNumber: col(0, 5, Column.Schema.int),
    residueName: col(5, 5, Column.Schema.str),
    atomName: col(10, 5, Column.Schema.str),
    atomNumber: col(15, 5, Column.Schema.int),
    x: col(pO, pW, Column.Schema.float),
    y: col(pO + pW, pW, Column.Schema.float),
    z: col(pO + 2 * pW, pW, Column.Schema.float),
    vx: hasVelocities ? col(vO, vW, Column.Schema.float) : undef,
    vy: hasVelocities ? col(vO + vW, vW, Column.Schema.float) : undef,
    vz: hasVelocities ? col(vO + 2 * vW, vW, Column.Schema.float) : undef,
  };

  return ret;
}

/**
 * box vectors (free format, space separated reals), values:
 * v1(x) v2(y) v3(z) v1(y) v1(z) v2(x) v2(z) v3(x) v3(y),
 * the last 6 values may be omitted (they will be set to zero).
 * Gromacs only supports boxes with v1(y)=v1(z)=v2(z)=0.
 */
function handleBoxVectors(state: State) {
  const { tokenizer } = state;
  // one line a frame: decoded, the values read by `+value` as before
  const values = TextByteTokenizer.readLine(tokenizer).trim().split(/\s+/g);
  const v = (i: number) => (values.length > i ? +values[i] : 0);
  state.header.box = [v(0), v(1), v(2)];
  state.header.boxVectors = [v(0), v(3), v(4), v(5), v(1), v(6), v(7), v(8), v(2)];
}

/** The lengths of the runs of digits after each '.' in the bytes `start` to `end` (a run of at least one digit). */
function decimalPlaces(data: ReadonlyUint8Array, start: number, end: number) {
  const ret: number[] = [];
  for (let i = start; i < end; i++) {
    if (data[i] !== 46) continue;
    let j = i + 1;
    while (j < end && data[j] >= 48 && data[j] <= 57) j++;
    if (j > i + 1) ret.push(j - i - 1);
    i = j - 1;
  }
  return ret;
}

async function parseInternal(data: ReadonlyUint8Array, ctx: RuntimeContext): Promise<Result<Schema.GroFile>> {
  const tokenizer = TextByteTokenizer(data);

  await ctx.update({ message: 'Parsing...', current: 0, max: data.length });
  const structures: Schema.GroStructure[] = [];
  while (tokenizer.position < data.length) {
    const state = State(tokenizer, ctx);
    handleTitleString(state);
    handleNumberOfAtoms(state);
    const atoms = await handleAtoms(state);
    handleBoxVectors(state);
    structures.push({ header: state.header, atoms });
  }

  const result: Schema.GroFile = { structures };
  return Result.success(result);
}

/** Parses GRO from its text or its UTF-8 bytes (read, not modified). */
export function parseGRO(data: StringLike | ReadonlyUint8Array) {
  return Task.create<Result<Schema.GroFile>>('Parse GRO', async (ctx) => {
    try {
      return await parseInternal(TextByteTokenizer.bytesOf(data), ctx);
    } catch (e) {
      if (e instanceof GroError) return Result.error<Schema.GroFile>(e.message);
      throw e;
    }
  });
}
