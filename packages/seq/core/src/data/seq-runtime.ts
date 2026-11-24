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
  SeqNodeKeyField,
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
// region by joining the pairs of its system and level that overlap it, each cut to its part of the region; its
// dynamic nodes are found by key when its coverage says so.

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
 * A `seq` source's dynamic nodes found by key (_spec/README.md "Data", The seq format), when its coverage lists `key` in
 * `names`: the index a lazy chunk (`seq-node-keys.ts`), loaded on the first find.
 */
function nodeKeys(data: SeqData): OpenedSource['findNames'] {
  if (!data.coverage.names?.includes(SeqNodeKeyField)) return undefined;
  return (text, options, signal) =>
    signal.aborted
      ? Promise.reject(abortError())
      : import('./seq-node-keys.js').then((m) => m.findNodeKeys(data, text, options, signal));
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
    // a sided system by its side and name (`a:chr1`), an unsided one by its name
    const at = (n: string) => (s.side ? `${s.side}:${n}` : n);
    names.set(at(s.name), s.name);
    for (const a of s.aliases ?? []) if (!names.has(at(a))) names.set(at(a), s.name);
  }
  // A stable row identity: each pair's list rows numbered after the rows of the pairs before it.
  const rowBase = new Map<SeqPair, number>();
  let rows = 0;
  for (const p of data.pairs) {
    rowBase.set(p, rows);
    const t = p.answer.ok ? p.answer.tables : undefined;
    if (t?.ranges?.form === 'list') rows += t.ranges.start.length;
    if (t?.pairs?.form === 'list') rows += t.pairs.aStart.length;
    if (t?.matrix?.form === 'items') rows += t.matrix.start.length;
  }

  const answer = (request: RegionRequest): RegionAnswer => {
    const side = request.side;
    const system = sourceNameOf(names, side ? `${side}:${request.system}` : request.system);
    if (system === undefined) return notServed('a system the coverage does not list');
    if (!data.coverage.levels.includes(request.level)) return notServed('a level the coverage does not list');
    const { start, end } = request;
    // an items tile is the finest level's items over the region: the level is only its width
    const level = request.items ? data.coverage.levels[0] : request.level;
    const ry = request.y;
    const ySystem = ry && sourceNameOf(names, ry.system);
    if (ry && ySystem === undefined) return notServed('a second system the coverage does not list');
    // a grid's pair holds the whole of the asked second axis (grid parts join along x only)
    const yHolds = (p: SeqPair) => {
      const py = p.request.y;
      return ry ? !!py && py.system === ySystem && py.start <= ry.start && py.end >= ry.end : !py;
    };
    const pairs = data.pairs
      .filter((p) => p.request.system === system && p.request.level === level && p.request.side === side)
      .filter((p) => p.request.start < end && p.request.end > start && yHolds(p))
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
        if (t)
          cuts.push(
            sliceTable(t, level, {
              start,
              end,
              from,
              to,
              rowBase: rowBase.get(pair),
              system,
              ...(side && { side }),
              ...(ry && { y: { start: ry.start, end: ry.end } }),
            }),
          );
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
  const findNames = nodeKeys(data);
  return {
    coverage: data.coverage,
    read(request, signal) {
      if (signal.aborted) return Promise.reject(abortError());
      return Promise.resolve(answer(request));
    },
    ...(findNames ? { findNames } : {}),
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
