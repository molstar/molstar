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
  type FiberEntryId,
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
/** The IndexedDB database kept records go to when the spec gives no rows (journal.md §1; History's Keep). */
export const RecordsDatabase = 'prism-viewer-records';
const Owner = fiberOwnerId('prism-viewer');

/** The fiber's largest window (`maxWindowEntries`): longer reads go in pages of it. */
export const RecordWindowMax = FIBER_LIMIT_HARD_BOUNDS.maxWindowEntries;
/** A turn under way says it is alive this often, well inside the fiber's inactivity timeout (journal.md §2, F52). */
const HeartbeatMs = 5 * 60_000;
/** A kept record's index entry is touched at most this often while it is appended to (and at a turn's end). */
const TouchMs = 60_000;

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

/** A read: the items newest first, the positions it covered (`start` inclusive, `end` exclusive), and whether there are earlier ones. */
export interface RecordPage {
  readonly items: RecordItem[];
  readonly more: boolean;
  readonly start: number;
  readonly end: number;
}

export interface TurnWriter {
  /** Ends the turn's entry (once; later calls do nothing). Never rejects: a failure is logged. */
  settle(end: { reply: string; outcome: 'done' | 'stopped' | 'failed'; model?: string; usage?: Json }): Promise<void>;
}

export interface SessionRecord {
  readonly id: string;
  /** Kept in IndexedDB (or the rows given to `keep`): it outlives the page. */
  readonly kept: boolean;
  /** Where "Keep" puts it: the spec's `persist.rows` (journal.md §1), when it gives them. */
  readonly keepIn?: FiberRows;
  /** Entries newest first: `limit` of them before `before` (a position), read in windows of at most `RecordWindowMax`. */
  read(limit: number, before?: number): Promise<RecordPage>;
  /**
   * Entries at `from` (a position) and after, newest first: the oldest `limit` of them, read in windows of at most
   * `RecordWindowMax`. `more`: later entries were left out. What History reads after a change, instead of a window.
   */
  readAfter(from: number, limit: number): Promise<RecordPage>;
  /** Called after each write; a listener that throws is logged and does not stop the others. */
  onChange(listener: () => void): () => void;
  /** An agent's turn, begun at its prompt; settled at its end. Alive (a heartbeat) until then. */
  beginTurn(turn: string, prompt: string, focus?: string): TurnWriter;
  /**
   * Copies the record into these rows and goes on there; it is listed in the index (journal.md §1). In the write order:
   * what is appended meanwhile waits for the copy and goes to the kept record. One keep at a time (a second call gets the
   * first's promise).
   */
  keep(rows: FiberRows, index: RecordIndex, title: string): Promise<void>;
  /**
   * Forgets the kept copy (its rows and its index entry: the Storage view's Forget, journal.md §7) and goes on in memory
   * with what it holds, so History still shows it and it can be kept again. In the write order; nothing when not kept.
   */
  unkeep(): Promise<void>;
  /** Writes in flight are done. */
  flush(): Promise<void>;
  /** Stops recording: a turn still open is settled as stopped, the index touched, the writes flushed. */
  dispose(): Promise<void>;
}

/** Operations whose runs are not changes worth a line: reads, hover, panels. */
const isQuiet = (r: RunRecord, read: boolean) =>
  read || r.operation.startsWith('highlight.') || r.operation.startsWith('view.') || r.operation === 'export.preview';

/** Looks like a URL: a scheme and `//`. */
const UrlLike = /^[a-z][a-z0-9+.-]*:\/\//i;

/**
 * A run's values as strict JSON, small: what JSON cannot carry (a loci, a file, a function) said in words; long text
 * cut; a URL without its query and fragment (signed URLs and tokens stay out of kept records and `.molx`).
 */
export function summarizeValues(values: Readonly<Record<string, unknown>>, max = 2000): Json {
  const seen = new WeakSet<object>();
  const walk = (v: unknown, depth: number): Json => {
    if (v === null || typeof v === 'boolean') return v;
    if (typeof v === 'number') return Number.isFinite(v) ? v : String(v);
    if (typeof v === 'string') {
      const s = UrlLike.test(v) ? v.replace(/[?#][\s\S]*$/, (q) => (q.length > 1 ? `${q[0]}…` : '')) : v;
      return s.length > 200 ? `${s.slice(0, 199)}…` : s;
    }
    if (typeof v === 'function') return '(function)';
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
  /** Where "Keep" puts it (the spec's `persist.rows`); the browser's IndexedDB when absent. */
  readonly keepIn?: FiberRows;
  /** Its id (one per viewer session); new by default. */
  readonly id?: string;
  /** The clock of the entries and of the fiber (its inactivity timeout); `Date.now` by default. */
  readonly clock?: () => number;
  /** Runs `fn` every `ms` until the returned function is called (`setInterval` by default): the turns' heartbeat. */
  readonly every?: (fn: () => void, ms: number) => () => void;
}

const everyInterval = (fn: () => void, ms: number) => {
  const handle = setInterval(fn, ms);
  return () => clearInterval(handle);
};

export async function createSessionRecord(o: SessionRecordOptions): Promise<SessionRecord> {
  const now = o.clock ?? Date.now;
  const clock = { now };
  const every = o.every ?? everyInterval;
  const id = o.id ?? `session-${now().toString(36)}-${Math.random().toString(36).slice(2, 8)}`;
  const fid = fiberId(id);
  let host: FiberHost = createFiberHost({ rows: o.rows ?? createMemoryRows(), clock });
  let kept = false;
  let disposed = false;
  let seq = 0;
  let keptIndex: { index: RecordIndex; title: string } | undefined;
  /** The index's last touch, and whether entries were appended since. */
  let touched = 0;
  let dirty = false;
  const listeners = new Set<() => void>();
  const changed = () => {
    for (const l of [...listeners]) {
      try {
        l();
      } catch (e) {
        console.error('[prism-viewer] record listener:', e);
      }
    }
  };
  const ok = <T extends { status: string }>(r: T, what: string) => {
    if (r.status !== 'accepted') throw new Error(`The record could not ${what}: ${JSON.stringify(r)}`);
    return r;
  };
  ok(
    await host.create({ commandId: fiberCommandId(`${id}:create`), fiberId: fid, ownerId: Owner, limits: limits() }),
    'start',
  );

  // One chain of writes: the record's order is the order things happened. `write` gives the work's own outcome (keep
  // needs its failure); the chain itself never rejects, so one failed write does not skip the next.
  let chain: Promise<void> = Promise.resolve();
  const write = <T>(work: () => Promise<T>, notify = true): Promise<T> => {
    const run = chain.then(work);
    chain = run.then(
      () => {
        if (notify) changed();
      },
      (e) => console.error('[prism-viewer] record:', e),
    );
    return run;
  };
  const quiet = (p: Promise<unknown>): Promise<void> =>
    p.then(
      () => undefined,
      () => undefined,
    );

  /** The kept record's index entry: on keep, then at most once a `TouchMs` while appended to, and at a turn's end. */
  const touch = async (force: boolean) => {
    if (!keptIndex || !dirty || (!force && now() - touched < TouchMs)) return;
    dirty = false;
    touched = now();
    await keptIndex.index.touch(id, keptIndex.title, touched);
  };

  const append = (entry: RecordEntry) => {
    if (disposed) return;
    const n = ++seq;
    void quiet(
      write(async () => {
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
        dirty = true;
        await touch(false);
      }),
    );
  };

  // The turns under way, each alive by a heartbeat while open (a turn parked for an approval waits as long as the
  // person does; the fiber would abort it after its inactivity timeout otherwise).
  const open = new Map<number, { readonly entryId: FiberEntryId; readonly stop: () => Promise<void> }>();
  let stopBeat: (() => void) | undefined;
  const beat = () =>
    void quiet(
      write(async () => {
        for (const t of open.values()) {
          const r = await host.heartbeat({ fiberId: fid, entryId: t.entryId });
          if ('status' in r) console.warn('[prism-viewer] record: a turn could not be kept alive', r.status);
        }
      }, false),
    );
  const beating = () => {
    if (open.size && !stopBeat) stopBeat = every(beat, HeartbeatMs);
    else if (!open.size && stopBeat) {
      stopBeat();
      stopBeat = undefined;
    }
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
      append({ kind: 'session', event, detail: r.result.did, at });
      return;
    }
    append({
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

  const undoneSub = o.undone?.subscribe((step) =>
    append({
      kind: 'undo',
      step: step.id,
      what: step.label,
      ...(step.origin.turn ? { turn: step.origin.turn } : {}),
      at: now(),
    }),
  );

  let keeping: Promise<void> | undefined;

  const record: SessionRecord = {
    id,
    get kept() {
      return kept;
    },
    ...(o.keepIn ? { keepIn: o.keepIn } : {}),
    async read(limit, before) {
      await chain;
      const items: RecordItem[] = [];
      let start = before;
      let end: number | undefined;
      let more = false;
      // Backward windows of at most the fiber's limit, until `limit` entries are read or the record begins.
      for (let left = Math.max(1, limit); left > 0; ) {
        const page = await host.window({
          fiberId: fid,
          limit: Math.min(left, RecordWindowMax),
          direction: 'backward',
          ...(start !== undefined ? { start } : {}),
        });
        if ('status' in page) {
          console.error('[prism-viewer] record: could not read', page);
          break;
        }
        end ??= page.endExclusive;
        for (let i = page.entries.length - 1; i >= 0; i--) {
          const item = itemOf(page.entries[i]);
          if (item) items.push(item);
        }
        more = page.hasMore;
        start = page.start;
        left -= page.entries.length;
        if (!page.hasMore || !page.entries.length) break;
      }
      return { items, more, start: start ?? 0, end: end ?? start ?? 0 };
    },
    async readAfter(from, limit) {
      await chain;
      const entries: FiberEntry[] = [];
      let start = from;
      let more = false;
      // Forward windows of at most the fiber's limit, until `limit` entries are read or the record ends.
      for (let left = Math.max(1, limit); left > 0; ) {
        const page = await host.window({
          fiberId: fid,
          start,
          limit: Math.min(left, RecordWindowMax),
          direction: 'forward',
        });
        if ('status' in page) {
          console.error('[prism-viewer] record: could not read', page);
          break;
        }
        for (const e of page.entries) entries.push(e);
        more = page.hasMore;
        if (page.entries.length) start = page.endExclusive;
        left -= page.entries.length;
        if (!page.hasMore || !page.entries.length) break;
      }
      const items: RecordItem[] = [];
      for (let i = entries.length - 1; i >= 0; i--) {
        const item = itemOf(entries[i]);
        if (item) items.push(item);
      }
      return { items, more, start: from, end: start };
    },
    onChange(listener) {
      listeners.add(listener);
      return () => listeners.delete(listener);
    },
    beginTurn(turn, prompt, focus) {
      if (disposed) return { settle: async () => {} };
      const n = ++seq;
      const entryId = fiberEntryId(`turn-${n}`);
      const begun: TurnEntry = { kind: 'turn', turn, prompt, ...(focus ? { focus } : {}), at: now() };
      void quiet(
        write(async () => {
          ok(
            await host.begin({
              commandId: fiberCommandId(`${id}:${n}`),
              fiberId: fid,
              entryId,
              descriptor: envelope(begun as unknown as Json),
            }),
            'begin a turn',
          );
          dirty = true;
        }),
      );
      const settle: TurnWriter['settle'] = (end) => {
        if (!open.delete(n)) return chain;
        beating();
        return quiet(
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
            dirty = true;
            await touch(true);
          }),
        );
      };
      open.set(n, { entryId, stop: () => settle({ reply: '', outcome: 'stopped' }) });
      beating();
      return { settle };
    },
    keep(rows, index, title) {
      if (kept) return Promise.resolve();
      if (disposed) return Promise.reject(new Error('The record has ended.'));
      keeping ??= write(async () => {
        if (kept) return;
        // Copied in order, the same entry ids and positions; a turn under way is copied as begun and settles there.
        const target = createFiberHost({ rows, clock });
        ok(
          await target.create({
            commandId: fiberCommandId(`${id}:create`),
            fiberId: fid,
            ownerId: Owner,
            limits: limits(),
          }),
          'be kept',
        );
        try {
          await copyFiber(host, target, fid, id, ok);
        } catch (e) {
          // Nothing half-kept: the copy goes, and the record goes on in memory.
          await target.forget(fid).catch(() => undefined);
          throw e;
        }
        host = target;
        kept = true;
        keptIndex = { index, title };
        dirty = true;
        try {
          await touch(true);
        } catch (e) {
          console.error('[prism-viewer] record: the index could not list it', e);
        }
      }).finally(() => {
        keeping = undefined;
      });
      return keeping;
    },
    unkeep() {
      if (!kept) return Promise.resolve();
      return write(async () => {
        if (!kept || !keptIndex) return;
        // Back into memory first, in order (as keep copies out); then the kept copy and its index entry go.
        const target = createFiberHost({ rows: createMemoryRows(), clock });
        ok(
          await target.create({
            commandId: fiberCommandId(`${id}:create`),
            fiberId: fid,
            ownerId: Owner,
            limits: limits(),
          }),
          'go on in memory',
        );
        await copyFiber(host, target, fid, id, ok);
        const from = host;
        const { index } = keptIndex;
        host = target;
        kept = false;
        keptIndex = undefined;
        dirty = false;
        await from.forget(fid);
        await index.forget(id);
      });
    },
    flush: () => chain,
    async dispose() {
      if (!disposed) {
        disposed = true;
        sub.unsubscribe();
        undoneSub?.unsubscribe();
        // A turn still open ends as stopped (journal.md §4), so it is neither "under way" forever nor left to abort.
        for (const t of [...open.values()]) void t.stop();
        stopBeat?.();
        stopBeat = undefined;
        void quiet(write(() => touch(true), false));
      }
      await chain;
      listeners.clear();
    },
  };
  return record;
}

/**
 * An entry as History reads it (settled or under way). An aborted turn (its turn never settled: the page went away while
 * it ran, the fiber's inactivity timeout) reads as stopped, from what it was begun with; other aborted entries are left
 * out (only turns stream).
 */
function itemOf(e: FiberEntry): RecordItem | undefined {
  if (e.state.kind === 'aborted') {
    const begun = e.descriptor.value as TurnEntry | undefined;
    if (!begun || typeof begun !== 'object' || begun.kind !== 'turn') return undefined;
    return { id: e.entryId, position: e.position, entry: { ...begun, outcome: 'stopped' }, open: false };
  }
  const value = (e.state.kind === 'settled' ? e.state.payload : e.descriptor).value as unknown as RecordEntry;
  if (!value || typeof value !== 'object' || !('kind' in value) || !('at' in value)) return undefined;
  return { id: e.entryId, position: e.position, entry: value, open: e.state.kind !== 'settled' };
}

/** Every entry of a fiber into another, in order, at the same positions (keep, and unkeep back into memory). */
async function copyFiber(
  from: FiberHost,
  to: FiberHost,
  fid: FiberId,
  id: string,
  ok: <T extends { status: string }>(r: T, what: string) => T,
): Promise<void> {
  let start = 0;
  for (;;) {
    const page = await from.window({ fiberId: fid, start, limit: RecordWindowMax, direction: 'forward' });
    if ('status' in page) throw new Error(`The record could not be read (${page.reason}).`);
    for (const e of page.entries) await copyEntry(to, fid, e, id, ok);
    if (!page.hasMore) break;
    start = page.endExclusive;
  }
}

/** One entry into the kept record, at the same position: settled as settled, under way as begun, aborted as aborted. */
async function copyEntry(
  target: FiberHost,
  fid: FiberId,
  e: FiberEntry,
  id: string,
  ok: <T extends { status: string }>(r: T, what: string) => T,
): Promise<void> {
  const commandId = fiberCommandId(`${id}:copy:${e.entryId}`);
  if (e.state.kind === 'settled') {
    ok(
      await target.append({
        commandId,
        fiberId: fid,
        entryId: e.entryId,
        descriptor: e.descriptor,
        payload: e.state.payload,
      }),
      'be kept',
    );
    return;
  }
  ok(await target.begin({ commandId, fiberId: fid, entryId: e.entryId, descriptor: e.descriptor }), 'be kept');
  if (e.state.kind === 'aborted') {
    ok(
      await target.abort({
        commandId: fiberCommandId(`${id}:copy:${e.entryId}:abort`),
        fiberId: fid,
        entryId: e.entryId,
        reason: e.state.reason,
      }),
      'be kept',
    );
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
  const put = async (value: Json) => {
    const r = await host.append({
      commandId: fiberCommandId(`records:${Date.now().toString(36)}:${++n}:${Math.random().toString(36).slice(2, 8)}`),
      fiberId: fid,
      entryId: fiberEntryId(`r-${Date.now().toString(36)}-${n}-${Math.random().toString(36).slice(2, 6)}`),
      descriptor: envelope({ kind: 'index' }),
      payload: envelope(value),
    });
    if (r.status !== 'accepted') throw new Error(`The records index could not be written (${r.status}).`);
  };
  // Every touch is an entry, folded here; folding the fiber into one entry per record when it grows is for later.
  return {
    async list() {
      const out = new Map<string, KeptRecord>();
      let start = 0;
      for (;;) {
        const page = await host.window({ fiberId: fid, start, limit: RecordWindowMax, direction: 'forward' });
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

/**
 * Forgets a kept record (the Storage view's Forget, journal.md §7): its rows, then its index entry. Not the current
 * session's record while it is kept: that one is `unkeep`'s.
 */
export async function forgetKept(rows: FiberRows, index: RecordIndex, id: string): Promise<void> {
  await createFiberHost({ rows }).forget(fiberId(id));
  await index.forget(id);
}

/**
 * A kept record's size, about: the characters of its rows as JSON (the fibers' default scope, `fiber:<id>`).
 * Undefined when its rows are not there.
 */
export async function keptSize(rows: FiberRows, id: string): Promise<number | undefined> {
  const loaded = await rows.load(`fiber:${id}`);
  if (!loaded) return undefined;
  let size = 0;
  for (const r of loaded.rows) size += r.key.length + JSON.stringify(r.row).length;
  return size;
}
