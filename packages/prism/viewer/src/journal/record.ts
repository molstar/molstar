/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import {
  createFiberHost,
  createFiberLimits,
  createMemoryRows,
  FIBER_LIMIT_HARD_BOUNDS,
  FIBER_PAYLOAD_VERSION,
  type FiberEntry,
  type FiberHost,
  type FiberId,
  type FiberRows,
  fiberCommandId,
  fiberEntryId,
  fiberId,
  fiberOwnerId,
} from '@molstar/ai-fibers';
import type { Json } from '@molstar/ai-agent';
import type { HistoryStep } from '@molstar/prism-ui/core';
import type { Observable } from 'rxjs';
import type { RunOrigin, ViewerApi } from '../api/area.js';
import type { RunRecord } from '../api/run.js';

// The session's record (journal.md): what the person did and what the agent did, in order, in a fiber
// (`@molstar/ai-fibers`). One record per viewer session, in memory by default; kept on purpose (IndexedDB rows) with a
// small index of kept records (journal.md §1; fibers have no listing by owner). One writer: every change the API runs
// (a `run` entry; an undo is an `undo` entry naming its step; a session opened, cleared or saved a `session` entry), and
// the agent's turns (`turn`, begun at the prompt, settled at the end). Append-only: nothing is edited in place.

export const RecordCodec = 'prism-viewer.record.v1';
const Owner = fiberOwnerId('prism-viewer');

export interface RunEntry {
  readonly kind: 'run';
  readonly operation: string;
  readonly title: string;
  readonly values: Json;
  readonly origin: { readonly kind: RunOrigin['kind']; readonly turn?: string };
  readonly result:
    | { readonly ok: true; readonly did: string }
    | { readonly ok: false; readonly code: string; readonly message: string };
  readonly at: number;
  readonly ms: number;
  readonly step?: number;
}

export interface TurnEntry {
  readonly kind: 'turn';
  readonly turn: string;
  readonly prompt: string;
  readonly focus?: string;
  readonly reply?: string;
  readonly outcome?: 'done' | 'stopped' | 'failed';
  readonly model?: string;
  readonly usage?: Json;
  readonly at: number;
}

export interface UndoEntry {
  readonly kind: 'undo';
  /** The step undone (a run's `step`), and its words. */
  readonly step: number;
  readonly what: string;
  /** The agent's turn it was part of. */
  readonly turn?: string;
  readonly at: number;
}

export interface SessionEntry {
  readonly kind: 'session';
  readonly event: 'opened' | 'cleared' | 'saved';
  readonly detail?: string;
  readonly at: number;
}

export type RecordEntry = RunEntry | TurnEntry | UndoEntry | SessionEntry;

/** An entry as History reads it: its value, and whether it is still being written (a turn under way). */
export interface RecordItem {
  readonly id: string;
  readonly position: number;
  readonly entry: RecordEntry;
  readonly open: boolean;
}

export interface TurnWriter {
  settle(end: { reply: string; outcome: 'done' | 'stopped' | 'failed'; model?: string; usage?: Json }): Promise<void>;
}

export interface SessionRecord {
  readonly id: string;
  /** Kept in IndexedDB (or the rows given to `keep`): it outlives the page. */
  readonly kept: boolean;
  /** Entries newest first: `limit` of them before `before` (a position), and whether there are more. */
  read(limit: number, before?: number): Promise<{ items: RecordItem[]; more: boolean }>;
  /** Called after each write. */
  onChange(listener: () => void): () => void;
  /** An agent's turn, begun at its prompt; settled at its end. */
  beginTurn(turn: string, prompt: string, focus?: string): TurnWriter;
  /** Copies the record into these rows and goes on there; it is listed in the index (journal.md §1). */
  keep(rows: FiberRows, index: RecordIndex, title: string): Promise<void>;
  /** Writes in flight are done. */
  flush(): Promise<void>;
  dispose(): void;
}

/** Operations whose runs are not changes worth a line: reads, hover, panels. */
const isQuiet = (r: RunRecord, read: boolean) =>
  read || r.operation.startsWith('highlight.') || r.operation.startsWith('view.') || r.operation === 'export.preview';

/** A run's values as strict JSON, small: what JSON cannot carry (a loci, a file) said in words; long text cut. */
export function summarizeValues(values: Readonly<Record<string, unknown>>, max = 2000): Json {
  const seen = new WeakSet<object>();
  const walk = (v: unknown, depth: number): Json => {
    if (v === null || typeof v === 'boolean') return v;
    if (typeof v === 'number') return Number.isFinite(v) ? v : String(v);
    if (typeof v === 'string') return v.length > 200 ? `${v.slice(0, 199)}…` : v;
    if (typeof v !== 'object') return String(v);
    if (seen.has(v)) return '(cycle)';
    seen.add(v);
    if (typeof File !== 'undefined' && v instanceof File) return `(file ${v.name})`;
    if (typeof Blob !== 'undefined' && v instanceof Blob) return `(${v.size} bytes)`;
    if ('loci' in (v as object)) return '(atoms picked in the view)';
    if (depth > 4) return '(…)';
    if (Array.isArray(v)) return v.slice(0, 20).map((x) => walk(x, depth + 1));
    const proto = Object.getPrototypeOf(v);
    if (proto !== Object.prototype && proto !== null) return `(${proto?.constructor?.name ?? 'object'})`;
    return Object.fromEntries(Object.entries(v as object).map(([k, x]) => [k, walk(x, depth + 1)]));
  };
  const out = walk(values, 0);
  return JSON.stringify(out).length > max ? '(values too large to keep)' : out;
}

const envelope = (value: Json) => ({ codec: RecordCodec, version: FIBER_PAYLOAD_VERSION as 1, value });
const limits = () => createFiberLimits({ ...FIBER_LIMIT_HARD_BOUNDS, maxActiveStreams: 8 });

export interface SessionRecordOptions {
  readonly api: ViewerApi;
  /** The steps the history undoes (an `undo` entry each, naming it). */
  readonly undone?: Observable<HistoryStep>;
  /** Where it is written (memory by default). */
  readonly rows?: FiberRows;
  /** Its id (one per viewer session); new by default. */
  readonly id?: string;
  readonly clock?: () => number;
}

export async function createSessionRecord(o: SessionRecordOptions): Promise<SessionRecord> {
  const now = o.clock ?? Date.now;
  const id = o.id ?? `session-${now().toString(36)}-${Math.random().toString(36).slice(2, 8)}`;
  const fid = fiberId(id);
  let host: FiberHost = createFiberHost({ rows: o.rows ?? createMemoryRows() });
  let kept = false;
  let seq = 0;
  let keptIndex: { index: RecordIndex; title: string } | undefined;
  const listeners = new Set<() => void>();
  const changed = () => {
    for (const l of [...listeners]) l();
  };
  const ok = <T extends { status: string }>(r: T, what: string) => {
    if (r.status !== 'accepted') throw new Error(`The record could not ${what}: ${JSON.stringify(r)}`);
    return r;
  };
  ok(
    await host.create({ commandId: fiberCommandId(`${id}:create`), fiberId: fid, ownerId: Owner, limits: limits() }),
    'start',
  );

  // One chain of writes: the record's order is the order things happened.
  let chain: Promise<void> = Promise.resolve();
  const write = (work: () => Promise<void>) => {
    chain = chain.then(work).then(changed, (e) => console.error('[prism-viewer] record:', e));
    return chain;
  };
  const append = (entry: RecordEntry) => {
    const n = ++seq;
    return write(async () => {
      ok(
        await host.append({
          commandId: fiberCommandId(`${id}:${n}`),
          fiberId: fid,
          entryId: fiberEntryId(`${entry.kind}-${n}`),
          descriptor: envelope({ kind: entry.kind }),
          payload: envelope(entry as unknown as Json),
        }),
        'append',
      );
      if (keptIndex) await keptIndex.index.touch(id, keptIndex.title, now());
    });
  };

  const sub = o.api.events.finished.subscribe((r) => {
    if (r.nested || !r.result) return;
    const def = o.api.describe(r.operation)?.definition;
    if (isQuiet(r, !!def?.read)) return;
    const at = r.started;
    const origin = r.origin.kind === 'agent' ? { kind: r.origin.kind, turn: r.origin.turn } : { kind: r.origin.kind };
    // Undos are written from the history's own event, naming the step (below).
    if (r.operation === 'history.undo' || r.operation === 'history.undo-turn') return;
    if (
      r.result.ok &&
      (r.operation === 'session.open' || r.operation === 'session.clear' || r.operation === 'session.save')
    ) {
      const event = r.operation === 'session.open' ? 'opened' : r.operation === 'session.clear' ? 'cleared' : 'saved';
      void append({ kind: 'session', event, detail: r.result.did, at });
      return;
    }
    void append({
      kind: 'run',
      operation: r.operation,
      title: r.title,
      values: summarizeValues(r.values),
      origin,
      result: r.result.ok
        ? { ok: true, did: r.result.did }
        : { ok: false, code: r.result.code, message: r.result.message },
      at,
      ms: Math.round(r.ms ?? 0),
      ...(r.step !== undefined ? { step: r.step } : {}),
    });
  });

  const undoneSub = o.undone?.subscribe(
    (step) =>
      void append({
        kind: 'undo',
        step: step.id,
        what: step.label,
        ...(step.origin.turn ? { turn: step.origin.turn } : {}),
        at: now(),
      }),
  );

  const record: SessionRecord = {
    id,
    get kept() {
      return kept;
    },
    async read(limit, before) {
      await chain;
      const page = await host.window({
        fiberId: fid,
        limit,
        direction: 'backward',
        ...(before !== undefined ? { start: before } : {}),
      });
      if ('status' in page) return { items: [], more: false };
      const items = page.entries
        .map(itemOf)
        .filter((x): x is RecordItem => !!x)
        .reverse();
      return { items, more: page.hasMore };
    },
    onChange(listener) {
      listeners.add(listener);
      return () => listeners.delete(listener);
    },
    beginTurn(turn, prompt, focus) {
      const n = ++seq;
      const entryId = fiberEntryId(`turn-${n}`);
      const begun: TurnEntry = { kind: 'turn', turn, prompt, ...(focus ? { focus } : {}), at: now() };
      void write(async () => {
        ok(
          await host.begin({
            commandId: fiberCommandId(`${id}:${n}`),
            fiberId: fid,
            entryId,
            descriptor: envelope(begun as unknown as Json),
          }),
          'begin a turn',
        );
      });
      return {
        settle: (end) =>
          write(async () => {
            const done: TurnEntry = {
              ...begun,
              reply: end.reply,
              outcome: end.outcome,
              ...(end.model ? { model: end.model } : {}),
              ...(end.usage ? { usage: end.usage } : {}),
            };
            ok(
              await host.settle({
                commandId: fiberCommandId(`${id}:${n}:settle`),
                fiberId: fid,
                entryId,
                payload: envelope(done as unknown as Json),
              }),
              'settle a turn',
            );
          }),
      };
    },
    async keep(rows, index, title) {
      await chain;
      if (kept) return;
      // Copied in order, the same entry ids; a turn under way is copied as begun and settles there.
      const target = createFiberHost({ rows });
      ok(
        await target.create({
          commandId: fiberCommandId(`${id}:create`),
          fiberId: fid,
          ownerId: Owner,
          limits: limits(),
        }),
        'be kept',
      );
      let start = 0;
      for (;;) {
        const page = await host.window({ fiberId: fid, start, limit: 200, direction: 'forward' });
        if ('status' in page) break;
        for (const e of page.entries) await copyEntry(target, fid, e, id);
        if (!page.hasMore) break;
        start = page.endExclusive;
      }
      host = target;
      kept = true;
      keptIndex = { index, title };
      await index.touch(id, title, now());
      changed();
    },
    flush: () => chain,
    dispose: () => {
      sub.unsubscribe();
      undoneSub?.unsubscribe();
    },
  };
  return record;
}

/** An entry as History reads it (settled or under way); an aborted one is left out. */
function itemOf(e: FiberEntry): RecordItem | undefined {
  if (e.state.kind === 'aborted') return undefined;
  const value = (e.state.kind === 'settled' ? e.state.payload : e.descriptor).value as unknown as RecordEntry;
  if (!value || typeof value !== 'object' || !('kind' in value) || !('at' in value)) return undefined;
  return { id: e.entryId, position: e.position, entry: value, open: e.state.kind !== 'settled' };
}

async function copyEntry(target: FiberHost, fid: FiberId, e: FiberEntry, id: string): Promise<void> {
  const commandId = fiberCommandId(`${id}:copy:${e.entryId}`);
  if (e.state.kind === 'settled') {
    await target.append({
      commandId,
      fiberId: fid,
      entryId: e.entryId,
      descriptor: e.descriptor,
      payload: e.state.payload,
    });
  } else if (e.state.kind === 'streaming') {
    await target.begin({ commandId, fiberId: fid, entryId: e.entryId, descriptor: e.descriptor });
  }
}

// ---- The index of kept records (journal.md §1) ----

/** A kept record as Home lists it. */
export interface KeptRecord {
  readonly id: string;
  readonly title: string;
  readonly first: number;
  readonly last: number;
}

export interface RecordIndex {
  list(): Promise<KeptRecord[]>;
  /** A record kept or appended to. */
  touch(id: string, title: string, at: number): Promise<void>;
  /** A record forgotten. */
  forget(id: string): Promise<void>;
}

/** The index as a fiber of its own in the same rows: append-only changes, folded when listed. */
export async function createRecordIndex(rows: FiberRows): Promise<RecordIndex> {
  const host = createFiberHost({ rows });
  const fid = fiberId('prism-viewer-records');
  const created = await host.create({
    commandId: fiberCommandId('records:create'),
    fiberId: fid,
    ownerId: Owner,
    limits: limits(),
  });
  if (created.status !== 'accepted' && created.status !== 'conflict')
    throw new Error('The records index could not start.');
  let n = 0;
  const put = (value: Json) =>
    host.append({
      commandId: fiberCommandId(`records:${Date.now().toString(36)}:${++n}:${Math.random().toString(36).slice(2, 8)}`),
      fiberId: fid,
      entryId: fiberEntryId(`r-${Date.now().toString(36)}-${n}-${Math.random().toString(36).slice(2, 6)}`),
      descriptor: envelope({ kind: 'index' }),
      payload: envelope(value),
    });
  return {
    async list() {
      const out = new Map<string, KeptRecord>();
      let start = 0;
      for (;;) {
        const page = await host.window({ fiberId: fid, start, limit: 200, direction: 'forward' });
        if ('status' in page) break;
        for (const e of page.entries) {
          if (e.state.kind !== 'settled') continue;
          const v = e.state.payload.value as { id: string; title?: string; at?: number; forget?: true };
          if (v.forget) out.delete(v.id);
          else {
            const had = out.get(v.id);
            out.set(v.id, {
              id: v.id,
              title: v.title ?? had?.title ?? v.id,
              first: had?.first ?? v.at ?? 0,
              last: v.at ?? 0,
            });
          }
        }
        if (!page.hasMore) break;
        start = page.endExclusive;
      }
      return [...out.values()].sort((a, b) => b.last - a.last);
    },
    async touch(id, title, at) {
      await put({ id, title, at });
    },
    async forget(id) {
      await put({ id, forget: true });
    },
  };
}
