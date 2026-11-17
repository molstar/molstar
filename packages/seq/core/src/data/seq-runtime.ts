/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type Issue, issueOf } from '@molstar/core/repr/tree';
import { SourceIssues } from './issues.js';
import {
  abortError,
  type FormatRuntime,
  type OpenContext,
  type OpenedSource,
  type RegionAnswer,
  type RegionRequest,
  type SourceEntry,
  type Tables,
} from './protocol.js';
import { decodeSeq, readSeq, type SeqData, type SeqPair } from './seq-format.js';
import { sliceTable } from './tables.js';

// The `seq` format's runtime (_spec/README.md "Data", The seq format): a source of request/answer pairs answers a
// region from the pair of its system and level that holds it, cut to the region.

/** A source that could not be opened, with its issues. */
export class SourceOpenError extends Error {
  constructor(readonly issues: readonly Issue[]) {
    super(issues[0]?.message ?? 'The source could not be opened');
    this.name = 'SourceOpenError';
  }
}

function notServed(why: string): RegionAnswer {
  return {
    ok: false,
    kind: 'not-served',
    issue: issueOf(SourceIssues, 'seq.source.not_served', `The source does not serve this region: ${why}`),
  };
}

function covers(p: SeqPair, system: string, start: number, end: number, level: number): boolean {
  return p.request.system === system && p.request.level === level && p.request.start <= start && p.request.end >= end;
}

function overlaps(p: SeqPair, system: string, start: number, end: number, level: number): boolean {
  return p.request.system === system && p.request.level === level && p.request.start < end && p.request.end > start;
}

/**
 * An opened source over `seq` data. A request names a system by name or alias, at a level the coverage lists: else it
 * is not served. It is answered by the first pair that holds the region (else the first that overlaps it), cut to it
 * and to the shapes asked; a region of a listed system that no pair touches is empty.
 */
export function openSeq(data: SeqData): OpenedSource {
  const names = new Map<string, string>();
  for (const s of data.coverage.systems) {
    names.set(s.name, s.name);
    for (const a of s.aliases ?? []) if (!names.has(a)) names.set(a, s.name);
  }
  const answer = (request: RegionRequest): RegionAnswer => {
    const system = names.get(request.system);
    if (system === undefined) return notServed('a system the coverage does not list');
    if (!data.coverage.levels.includes(request.level)) return notServed('a level the coverage does not list');
    const { start, end, level } = request;
    const pair =
      data.pairs.find((p) => covers(p, system, start, end, level)) ??
      data.pairs.find((p) => overlaps(p, system, start, end, level));
    if (!pair) return { ok: true, tables: {} };
    if (!pair.answer.ok) return pair.answer;
    const tables: { -readonly [S in keyof Tables]: Tables[S] } = {};
    for (const shape of request.shapes) {
      const t = pair.answer.tables[shape];
      if (t) (tables as Record<string, unknown>)[shape] = sliceTable(t, start, end, level);
    }
    return pair.answer.nodes ? { ok: true, tables, nodes: pair.answer.nodes } : { ok: true, tables };
  };
  return {
    coverage: data.coverage,
    read(request, signal) {
      if (signal.aborted) return Promise.reject(abortError());
      return Promise.resolve(answer(request));
    },
    dispose() {},
  };
}

/** Reads a `seq` source: its inline `value`, else its file at `uri` (JSON) through the asset port. */
async function open(source: SourceEntry, ctx: OpenContext): Promise<OpenedSource> {
  let read: { data?: SeqData; issues: readonly Issue[] };
  if (source.value !== undefined) read = readSeq(source.value);
  else if (source.uri !== undefined) read = decodeSeq(await ctx.assets.readText(source.uri, ctx.signal));
  else read = { issues: [issueOf(SourceIssues, 'seq.source.no_data', 'A seq source needs a value or a uri')] };
  if (!read.data) throw new SourceOpenError(read.issues.map((i) => ({ ...i, address: i.address ?? source.address })));
  return openSeq(read.data);
}

/** The runtime of the `seq` format. */
export const SeqRuntime: FormatRuntime = { name: 'seq', formats: ['seq'], open };
