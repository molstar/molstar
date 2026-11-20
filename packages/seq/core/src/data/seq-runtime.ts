/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type Issue, issueOf } from '@molstar/core/repr/tree';
import { SourceIssues } from './issues.js';
import {
  abortError,
  type DynamicNode,
  type FormatRuntime,
  isAbortError,
  type OpenContext,
  type OpenedSource,
  type RegionAnswer,
  type RegionRequest,
  type SourceEntry,
  type Table,
  type Tables,
} from './protocol.js';
import { decodeSeq, readSeq, type SeqData, type SeqPair } from './seq-format.js';
import { joinTables, sliceTable } from './tables.js';
import { sourceNameOf } from './system-names.js';

// The `seq` format's runtime (_spec/README.md "Data", The seq format): a source of request/answer pairs answers a
// region by joining the pairs of its system and level that overlap it, each cut to its part of the region.

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

function malformed(why: string): RegionAnswer {
  return {
    ok: false,
    kind: 'malformed',
    issue: issueOf(SourceIssues, 'seq.source.malformed', `The source's pairs do not join: ${why}`),
  };
}

/**
 * An opened source over `seq` data. A request names a system by name or alias, at a level the coverage lists: else it
 * is not served. The pairs of that system and level overlapping the region are joined in position order, each
 * answering for its part of the region; where none does (a gap), the region is not served, naming the first gap. A
 * failure among them is the answer. A list without ids gets `row` from its pair's place, stable across tiles.
 */
export function openSeq(data: SeqData): OpenedSource {
  const names = new Map<string, string>();
  for (const s of data.coverage.systems) {
    names.set(s.name, s.name);
    for (const a of s.aliases ?? []) if (!names.has(a)) names.set(a, s.name);
  }
  // A stable row identity: each pair's list rows numbered after the rows of the pairs before it.
  const rowBase = new Map<SeqPair, number>();
  let rows = 0;
  for (const p of data.pairs) {
    rowBase.set(p, rows);
    const r = p.answer.ok ? p.answer.tables.ranges : undefined;
    if (r?.form === 'list') rows += r.start.length;
  }

  const answer = (request: RegionRequest): RegionAnswer => {
    const system = sourceNameOf(names, request.system);
    if (system === undefined) return notServed('a system the coverage does not list');
    if (!data.coverage.levels.includes(request.level)) return notServed('a level the coverage does not list');
    const { start, end, level } = request;
    const pairs = data.pairs
      .filter((p) => p.request.system === system && p.request.level === level)
      .filter((p) => p.request.start < end && p.request.end > start)
      .sort((a, b) => a.request.start - b.request.start);
    const parts: { pair: SeqPair; from: number; to: number }[] = [];
    let at = start;
    for (const pair of pairs) {
      if (at >= end) break;
      if (pair.request.end <= at) continue;
      if (pair.request.start > at) return notServed(`no pair holds positions ${at} to ${pair.request.start}`);
      const to = Math.min(end, pair.request.end);
      parts.push({ pair, from: at, to });
      at = to;
    }
    if (at < end) return notServed(`no pair holds positions ${at} to ${end}`);
    for (const p of parts) if (!p.pair.answer.ok) return p.pair.answer;

    const tables: Record<string, Table> = {};
    for (const shape of request.shapes) {
      const cuts: Table[] = [];
      for (const { pair, from, to } of parts) {
        const t = pair.answer.ok ? pair.answer.tables[shape] : undefined;
        if (t) cuts.push(sliceTable(t, level, { start, end, from, to, rowBase: rowBase.get(pair) }));
      }
      if (!cuts.length) continue;
      const joined = joinTables(cuts, level);
      if (typeof joined === 'string') return malformed(joined);
      tables[shape] = joined;
    }
    const nodes = new Map<string, DynamicNode>();
    for (const { pair } of parts) if (pair.answer.ok) for (const n of pair.answer.nodes ?? []) nodes.set(n.key, n);
    const out = tables as Tables;
    return nodes.size ? { ok: true, tables: out, nodes: [...nodes.values()] } : { ok: true, tables: out };
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

/** A MessagePack map's first byte (fixmap, map 16, map 32): the only start of a `seq` file in MessagePack. */
const isMapMarker = (b: number | undefined) => b !== undefined && ((b & 0xf0) === 0x80 || b === 0xde || b === 0xdf);

/**
 * A `seq` file's bytes: MessagePack when its first byte is a map's marker, else JSON (an HTML page or truncated JSON
 * reports JSON's error), each through the tree's codec. The MessagePack codec is loaded only when a file needs it, so
 * a reader of inline specs does not pay for it.
 */
async function decodeSeqBytes(bytes: Uint8Array): Promise<{ data?: SeqData; issues: readonly Issue[] }> {
  if (!isMapMarker(bytes[0])) return decodeSeq(new TextDecoder('utf-8').decode(bytes));
  const { MsgPackCodec } = await import('@molstar/core/repr/codecs/msgpack');
  const decoded = MsgPackCodec.decode(bytes);
  return decoded.node ? readSeq(decoded.node) : { issues: decoded.issues };
}

/** Reads a `seq` source: its inline `value`, else its file at `uri` (JSON or MessagePack) through the asset port. */
async function open(source: SourceEntry, ctx: OpenContext): Promise<OpenedSource> {
  let read: { data?: SeqData; issues: readonly Issue[] };
  if (source.value !== undefined) read = readSeq(source.value);
  else if (source.uri !== undefined) {
    let bytes: Uint8Array;
    try {
      bytes = await ctx.assets.readBytes(source.uri, ctx.signal);
    } catch (e) {
      if (isAbortError(e) || ctx.signal.aborted) throw e;
      const issue = issueOf(SourceIssues, 'seq.source.unreadable', 'The source’s file could not be read', {
        address: source.address,
      });
      throw new SourceOpenError([issue]);
    }
    // bounded before it is parsed (_spec/README.md "Data", Protocol: limits)
    if (bytes.byteLength > ctx.limits.maxBytes)
      throw new SourceOpenError([
        issueOf(
          SourceIssues,
          'seq.source.too_large',
          `The file has ${bytes.byteLength} bytes, more than ${ctx.limits.maxBytes}`,
          { address: source.address, hint: 'Split the file, or serve it by region' },
        ),
      ]);
    read = await decodeSeqBytes(bytes);
  } else read = { issues: [issueOf(SourceIssues, 'seq.source.no_data', 'A seq source needs a value or a uri')] };
  if (!read.data) throw new SourceOpenError(read.issues.map((i) => ({ ...i, address: i.address ?? source.address })));
  return openSeq(read.data);
}

/** The runtime of the `seq` format. */
export const SeqRuntime: FormatRuntime = { name: 'seq', formats: ['seq'], open };
