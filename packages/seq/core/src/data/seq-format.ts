/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { parseJson, printJson } from '@molstar/core/repr/codecs/json';
import {
  ColumnTypes,
  columnType,
  type Issue,
  issueOf,
  isNode,
  isRecord,
  type Node,
  type TreePath,
  type Value,
  type ValueRecord,
} from '@molstar/core/repr/tree';
import { SourceIssues } from './issues.js';
import { DefaultTileBins } from './keys.js';
import type {
  Coverage,
  DynamicNode,
  FieldInfo,
  RegionAnswer,
  RegionRequest,
  Shape,
  ShapeInfo,
  SystemInfo,
  Tables,
} from './protocol.js';
import { FailureKinds, Shapes } from './protocol.js';
import { checkTableBounds, checkTables, DefaultMaxColumnLength } from './tables.js';

// The `seq` format (_spec/README.md "Data", The seq format): the protocol written down, a coverage and the
// request/answer pairs of a source. Written as a tree value (`(seq :coverage … :pairs […])`), so the JSON codec
// carries its typed columns and atoms losslessly; inline it is the same record as a spec value.

export interface SeqPair {
  readonly request: RegionRequest;
  readonly answer: RegionAnswer;
}

export interface SeqData {
  readonly coverage: Coverage;
  readonly pairs: readonly SeqPair[];
}

/** The node kind of a `seq` file. */
export const SeqKind = 'seq';

/** An answer as a value: a failure's `kind` is written `failure` (a record with `kind` reads as a node). */
function answerToValue(a: RegionAnswer): ValueRecord {
  if (a.ok) return a as unknown as ValueRecord;
  return { ok: false, failure: a.kind, issue: a.issue as unknown as ValueRecord };
}

/** The data as a tree value: the record a spec's inline `value` holds. */
export function seqToValue(data: SeqData): ValueRecord {
  return {
    coverage: data.coverage as unknown as ValueRecord,
    pairs: data.pairs.map((p) => ({ request: p.request as unknown as ValueRecord, answer: answerToValue(p.answer) })),
  };
}

/** The data as a `seq` node, the root of a file. */
export function seqToNode(data: SeqData): Node {
  return { kind: SeqKind, params: seqToValue(data) };
}

/** A `seq` file as JSON text. */
export function encodeSeq(data: SeqData, options: { indent?: number } = {}): string {
  return printJson(seqToNode(data), options);
}

/** A `seq` file read from JSON text: the data, or the issues (the codec's, or `seq.source.malformed`). */
export function decodeSeq(text: string): { data?: SeqData; issues: readonly Issue[] } {
  const decoded = parseJson(text);
  if (!decoded.node) return { issues: decoded.issues };
  return readSeq(decoded.node);
}

class Malformed extends Error {
  constructor(
    readonly path: TreePath,
    message: string,
  ) {
    super(message);
  }
}

function fail(path: TreePath, message: string): never {
  throw new Malformed(path, message);
}

function rec(v: unknown, path: TreePath): ValueRecord {
  if (!isRecord(v as Value)) fail(path, 'expected a record');
  return v as ValueRecord;
}

function list(v: unknown, path: TreePath): readonly Value[] {
  if (!Array.isArray(v)) fail(path, 'expected a list');
  return v;
}

function str(v: unknown, path: TreePath): string {
  if (typeof v !== 'string') fail(path, 'expected a string');
  return v;
}

function int(v: unknown, path: TreePath, min = Number.NEGATIVE_INFINITY): number {
  if (typeof v !== 'number' || !Number.isInteger(v) || v < min) fail(path, `expected an integer of at least ${min}`);
  return v;
}

function opt<T>(r: ValueRecord, key: string, path: TreePath, read: (v: unknown, p: TreePath) => T): T | undefined {
  return r[key] === undefined ? undefined : read(r[key], [...path, key]);
}

function oneOf<T extends string>(v: unknown, of: readonly T[], path: TreePath): T {
  if (!of.includes(v as T)) fail(path, `expected one of ${of.join(', ')}`);
  return v as T;
}

function readSystem(v: unknown, path: TreePath): SystemInfo {
  const r = rec(v, path);
  str(r.name, [...path, 'name']);
  opt(r, 'length', path, (x, p) => {
    const n = int(x, p, 0);
    if (n > DefaultMaxColumnLength) fail(p, `a length of more than ${DefaultMaxColumnLength}`);
    return n;
  });
  opt(r, 'alphabet', path, str);
  opt(r, 'cyclic', path, (x, p) => (typeof x === 'boolean' ? x : fail(p, 'expected a boolean')));
  opt(r, 'aliases', path, (x, p) => list(x, p).map((a, i) => str(a, [...p, i])));
  return r as unknown as SystemInfo;
}

function readField(v: unknown, path: TreePath): FieldInfo {
  const r = rec(v, path);
  str(r.name, [...path, 'name']);
  oneOf(r.type, ColumnTypes, [...path, 'type']);
  opt(r, 'doc', path, str);
  return r as unknown as FieldInfo;
}

function readDynamicNode(v: unknown, path: TreePath): DynamicNode {
  const r = rec(v, path);
  str(r.key, [...path, 'key']);
  if (!isNode(r.node as Value)) fail([...path, 'node'], 'expected a node');
  return r as unknown as DynamicNode;
}

function readNodes(r: ValueRecord, path: TreePath): readonly DynamicNode[] | undefined {
  return opt(r, 'nodes', path, (x, p) => list(x, p).map((n, i) => readDynamicNode(n, [...p, i])));
}

function readCoverage(v: unknown, path: TreePath): Coverage {
  const r = rec(v, path);
  const systems = list(r.systems, [...path, 'systems']).map((s, i) => readSystem(s, [...path, 'systems', i]));
  const names = new Set<string>();
  systems.forEach((s, i) => {
    if (names.has(s.name)) fail([...path, 'systems', i, 'name'], 'a system listed twice');
    names.add(s.name);
  });
  list(r.shapes, [...path, 'shapes']).forEach((s, i) => {
    const p = [...path, 'shapes', i];
    const sr = rec(s, p);
    if (!names.has(str(sr.system, [...p, 'system']))) fail([...p, 'system'], 'a system the coverage does not list');
    oneOf(sr.shape, Shapes, [...p, 'shape']);
    opt(sr, 'fields', p, (x, fp) => list(x, fp).map((f, j) => readField(f, [...fp, j])));
  });
  const levels = list(r.levels, [...path, 'levels']).map((l, i) => int(l, [...path, 'levels', i], 1));
  if (levels.length === 0) fail([...path, 'levels'], 'expected at least one level');
  for (let i = 1; i < levels.length; i++)
    if (levels[i] <= levels[i - 1]) fail([...path, 'levels', i], 'levels are not ascending');
  readNodes(r, path);
  return r as unknown as Coverage;
}

function readRequest(v: unknown, path: TreePath): RegionRequest {
  const r = rec(v, path);
  str(r.system, [...path, 'system']);
  const start = int(r.start, [...path, 'start'], 0);
  if (int(r.end, [...path, 'end'], 0) < start) fail([...path, 'end'], 'the end is before the start');
  int(r.level, [...path, 'level'], 1);
  list(r.shapes, [...path, 'shapes']).forEach((s, i) => oneOf(s, Shapes, [...path, 'shapes', i]));
  return r as unknown as RegionRequest;
}

const Severities = ['error', 'warning', 'info'] as const;

function readAnswer(v: unknown, path: TreePath): RegionAnswer {
  const r = rec(v, path);
  if (r.ok === true) {
    const tables = rec(r.tables, [...path, 'tables']) as unknown as Tables;
    const bad = checkTables(tables);
    if (bad) fail([...path, 'tables'], bad);
    for (const [shape, t] of Object.entries(tables)) {
      const big = t && checkTableBounds(t);
      if (big) fail([...path, 'tables', shape], big);
    }
    const nodes = readNodes(r, path);
    return nodes ? { ok: true, tables, nodes } : { ok: true, tables };
  }
  if (r.ok !== false) fail([...path, 'ok'], 'expected true or false');
  const kind = oneOf(r.failure, FailureKinds, [...path, 'failure']);
  const ir = rec(r.issue, [...path, 'issue']);
  str(ir.code, [...path, 'issue', 'code']);
  oneOf(ir.severity, Severities, [...path, 'issue', 'severity']);
  str(ir.message, [...path, 'issue', 'message']);
  return { ok: false, kind, issue: ir as unknown as Issue };
}

/**
 * `seq` data from a tree value: a `seq` node (a file's root) or its params' record (inline). Checked to the
 * protocol's shapes; the first misfit is a `seq.source.malformed` issue at its path.
 */
export function readSeq(value: Value): { data?: SeqData; issues: readonly Issue[] } {
  try {
    let r: ValueRecord;
    if (isNode(value)) {
      if (value.kind !== SeqKind) fail(['kind'], `expected a ${SeqKind} node`);
      r = value.params ?? {};
    } else r = rec(value, []);
    const coverage = readCoverage(r.coverage, ['coverage']);
    const pairs = list(r.pairs, ['pairs']).map((p, i) => {
      const pr = rec(p, ['pairs', i]);
      return {
        request: readRequest(pr.request, ['pairs', i, 'request']),
        answer: readAnswer(pr.answer, ['pairs', i, 'answer']),
      };
    });
    return { data: { coverage, pairs }, issues: [] };
  } catch (e) {
    if (!(e instanceof Malformed)) throw e;
    return { issues: [malformed(e)] };
  }
}

/**
 * A coverage checked as `readSeq` checks a file's (systems with unique names and lengths of at least 0, shapes of
 * listed systems, integer levels of at least 1, ascending, at least one): the first misfit as an issue.
 */
export function checkCoverage(coverage: Coverage): Issue | undefined {
  try {
    readCoverage(coverage, ['coverage']);
    return undefined;
  } catch (e) {
    if (!(e instanceof Malformed)) throw e;
    return malformed(e);
  }
}

function malformed(e: Malformed): Issue {
  return issueOf(SourceIssues, 'seq.source.malformed', `The coverage or data is malformed: ${e.message}`, {
    path: e.path,
  });
}

/** The shapes a table record holds, in the protocol's order. */
function shapesOf(tables: Tables): Shape[] {
  return Shapes.filter((s) => tables[s] !== undefined);
}

/** Levels between an inline list's level 1 and its whole level grow by this factor (the summaries' too). */
const WholeFactor = 16;

/**
 * The level whose one tile (`tileBins` bins) covers `length` positions: the least power of 16 that does, 1 for a
 * system within a tile at level 1.
 */
function wholeLevel(length: number, tileBins = DefaultTileBins): number {
  let level = 1;
  while (level * tileBins < length) level *= WholeFactor;
  return level;
}

/**
 * Inline data as `seq` data: one system, one region over all of it at level 1 (its length the system's, else the
 * tables' extent). The fields of a ranges table are listed in the coverage. Ranges alone as a list over a system
 * longer than a tile are served at the level whose one tile covers the system too (`wholeLevel`, the same rows): a
 * reader of the whole system reads that tile, never the system a tile at level 1 (81,307 keys at 83 Mb).
 */
export function inlineSeq(system: SystemInfo, tables: Tables, nodes?: readonly DynamicNode[]): SeqData {
  const shapes = shapesOf(tables);
  const length = system.length ?? extentOf(tables);
  const shapeInfo: ShapeInfo[] = shapes.map((shape) => {
    const r = tables.ranges;
    if (shape !== 'ranges' || r?.form !== 'list' || !r.fields) return { system: system.name, shape };
    const fields = r.fields.map((f) => ({ name: f.name, type: columnType(f.column) ?? 'float64' }));
    return { system: system.name, shape, fields };
  });
  const answer: RegionAnswer = nodes ? { ok: true, tables, nodes } : { ok: true, tables };
  const whole = shapes.length === 1 && tables.ranges?.form === 'list' ? wholeLevel(length) : 1;
  const levels = whole > 1 ? [1, whole] : [1];
  return {
    coverage: { systems: [{ ...system, length }], shapes: shapeInfo, levels },
    pairs: levels.map((level) => ({ request: { system: system.name, start: 0, end: length, level, shapes }, answer })),
  };
}

function extentOf(tables: Tables): number {
  let end = 0;
  const { letters: l, ranges: r, values: v } = tables;
  if (l) end = Math.max(end, l.start + l.letters.length);
  if (r?.form === 'list') for (let i = 0; i < r.end.length; i++) end = Math.max(end, r.end[i]);
  if (v?.form === 'positions') end = Math.max(end, v.start + v.value.length);
  return end;
}
