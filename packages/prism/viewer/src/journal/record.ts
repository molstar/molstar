/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import {
  createThreadHost,
  createThreadLimits,
  createMemoryRows,
  THREAD_LIMIT_HARD_BOUNDS,
  THREAD_PAYLOAD_VERSION,
  threadCommandId,
  type ThreadEntry,
  type ThreadEntryId,
  threadEntryId,
  type ThreadHost,
  type ThreadId,
  threadId,
  threadOwnerId,
  type ThreadRows,
  type ThreadTable,
  type Json,
} from '@molstar/ai';
import type { HistoryStep } from '@molstar/prism-ui/core';
import type { Observable } from 'rxjs';
import type { Changed, RunOrigin, ViewerApi, OperationDefinition } from '../api/area.js';
import type { RunRecord } from '../api/run.js';

// The session's record (journal.md): what the person did and what the agent did, in order, in a thread
// (`@molstar/ai`). One record per viewer session, in memory by default; kept on purpose (IndexedDB rows) with a
// small index of kept records (journal.md §1; threads have no listing by owner). One writer: every change the API runs
// (a `run` entry; an undo is an `undo` entry naming its step; a session opened, cleared or saved a `session` entry), and
// the agent's turns (`turn`, begun at the prompt, settled at the end). Append-only: nothing is edited in place.

export const RecordCodec = 'prism-viewer.record.v1';
/** The viewer's default persist key: what the records database is named after when the spec has none. */
const DefaultPersistKey = 'prism-viewer';
/** The IndexedDB database of a viewer without a `persist.key` (and of the app's, whose key is the default). */
export const RecordsDatabase = `${DefaultPersistKey}-records`;
/**
 * The IndexedDB database kept records go to when the spec gives no rows (journal.md §1, §6; History's Keep): named after
 * the viewer's `persist.key`, so viewers on a page keep their records apart (README decision 29). Viewers with no key
 * share the default one.
 */
export const recordsDatabase = (persistKey?: string): string => `${persistKey || DefaultPersistKey}-records`;
const Owner = threadOwnerId('prism-viewer');

/** The thread's largest window (`maxWindowEntries`): longer reads go in pages of it. */
export const RecordWindowMax = THREAD_LIMIT_HARD_BOUNDS.maxWindowEntries;
/** A turn under way says it is alive this often, well inside the thread's inactivity timeout (journal.md §2). */
const HeartbeatMs = 5 * 60_000;
/** A kept record's index entry is touched at most this often while it is appended to (and at a turn's end). */
const TouchMs = 60_000;
/** A reply streams into its turn at most this often (journal.md §2). */
export const ReplyMs = 300;
/** A reply longer than this is not previewed (the final text carries it), nor are more deltas than this written for one turn. */
const ReplyPreviewMax = 32_000;
const ReplyDeltasMax = 1_000;
/** A run's `changed` keeps this many node refs; the rest are counted. */
export const ChangedNodesMax = 8;
/** A run's target text is cut at this length. */
const TargetMax = 80;
/** The values that say what a run was pointed at, in the order History reads them. */
const TargetKeys = ['structure', 'component', 'node', 'target'] as const;

export interface RunEntry {
  readonly kind: 'run';
  readonly operation: string;
  readonly title: string;
  readonly values: Json;
  readonly origin: RecordActor;
  readonly result:
    | { readonly ok: true; readonly did: string }
    | { readonly ok: false; readonly code: string; readonly message: string };
  /** What it was pointed at, in words: the `structure`, `component`, `node` and `target` values it ran with (`targetsOf`). */
  readonly targets?: readonly string[];
  /** What it changed, as the operation said (`changedOf`): the nodes (a few refs, the rest counted), the selection, the view. */
  readonly changed?: RunChanged;
  readonly at: number;
  readonly ms: number;
  readonly step?: number;
}

/** A run's `changed`, small (journal.md §2): at most `ChangedNodesMax` node refs, `more` counts the rest. */
export interface RunChanged {
  readonly nodes?: readonly string[];
  readonly more?: number;
  readonly selection?: true;
  readonly view?: true;
  readonly sequence?: true;
}

/** Who did a thing: a run's origin, and the actor of an undo or a redo. */
export interface RecordActor {
  readonly kind: RunOrigin['kind'];
  readonly turn?: string;
}

export interface TurnEntry {
  readonly kind: 'turn';
  readonly turn: string;
  readonly prompt: string;
  readonly focus?: string;
  readonly reply?: string;
  readonly outcome?: 'done' | 'stopped' | 'failed';
  /** The model that ran it, `provider:model` (known at the prompt, so a turn under way has it). */
  readonly model?: string;
  readonly usage?: Json;
  readonly at: number;
}

export interface UndoEntry {
  readonly kind: 'undo';
  /** The step undone (a run's `step`), and its words. An Undo turn is one entry: the newest step here, all in `steps`. */
  readonly step: number;
  /** Every step an Undo turn took back (newest first), when it was more than one. */
  readonly steps?: readonly number[];
  readonly what: string;
  /** The agent's turn it was part of. */
  readonly turn?: string;
  /** Who undid it (the origin of the call); absent when the history was driven outside the API. */
  readonly by?: RecordActor;
  readonly at: number;
}

/** A step Redo put back: after its `undo`, History shows it done again. */
export interface RedoEntry {
  readonly kind: 'redo';
  readonly step: number;
  readonly what: string;
  readonly turn?: string;
  /** Who redid it. */
  readonly by?: RecordActor;
  readonly at: number;
}

export interface SessionEntry {
  readonly kind: 'session';
  readonly event: 'opened' | 'cleared' | 'saved';
  readonly detail?: string;
  readonly at: number;
}

export type RecordEntry = RunEntry | TurnEntry | UndoEntry | RedoEntry | SessionEntry;

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
  /**
   * The reply so far (the model's text of the call under way): written to the turn's stream as the text appended
   * since, coalesced (at most once per `ReplyMs`), so History shows it as it arrives. The final text of `settle`
   * replaces it. A reply too long to preview, or a stream that cannot be written, is simply not shown until the end.
   */
  reply(text: string): void;
  /** Ends the turn's entry (once; later calls do nothing). Never rejects: a failure is logged. */
  settle(end: { reply: string; outcome: 'done' | 'stopped' | 'failed'; model?: string; usage?: Json }): Promise<void>;
}

export interface SessionRecord {
  readonly id: string;
  /** Kept in IndexedDB (or the rows given to `keep`): it outlives the page. */
  readonly kept: boolean;
  /** Where "Keep" puts it: the spec's `persist.rows` (journal.md §1), when it gives them. */
  readonly keepIn?: ThreadRows;
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
  beginTurn(turn: string, prompt: string, focus?: string, model?: string): TurnWriter;
  /**
   * Copies the record into these rows and goes on there; it is listed in the index (journal.md §1). In the write order:
   * what is appended meanwhile waits for the copy and goes to the kept record. One keep at a time (a second call gets the
   * first's promise). `release`: the rows are the record's to let go (a connection to close): called when the keep
   * fails or is not needed, and when the kept copy is forgotten or the record ends.
   */
  keep(rows: ThreadRows, index: RecordIndex, title: string, release?: () => void): Promise<void>;
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
const isQuiet = (r: RunRecord, def: OperationDefinition<any, any> | undefined) =>
  !!def?.read ||
  def?.record === false ||
  r.operation.startsWith('highlight.') ||
  r.operation.startsWith('view.') ||
  r.operation === 'export.preview';

/** A URL anywhere in a text: a scheme, `//` and what follows up to a space or a quote. */
const UrlIn = /(?<![a-z0-9+.-])[a-z][a-z0-9+.-]*:\/\/[^\s"'<>]*/gi;

/** A URL without its credentials (`user:token@`), query and fragment. */
const bareUrl = (url: string) =>
  url
    .replace(/^([a-z][a-z0-9+.-]*:\/\/)[^/?#]*@/i, '$1')
    .replace(/[?#][\s\S]*$/, (q) => (q.length > 1 ? `${q[0]}…` : ''));

/**
 * A run's values as strict JSON, small: what JSON cannot carry (a loci, a file, a function) said in words; long text
 * cut; a URL (anywhere in a text) without its credentials, query and fragment (signed URLs and tokens stay out of kept records and `.molx`).
 */
export function summarizeValues(values: Readonly<Record<string, unknown>>, max = 2000): Json {
  const seen = new WeakSet<object>();
  const walk = (v: unknown, depth: number): Json => {
    if (v === null || typeof v === 'boolean') return v;
    if (typeof v === 'number') return Number.isFinite(v) ? v : String(v);
    if (typeof v === 'string') {
      const s = v.includes('://') ? v.replace(UrlIn, bareUrl) : v;
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

/** A text with its URLs bare (no credentials, query, fragment), cut at `max`. */
function plainText(v: string, max: number): string {
  const s = v.includes('://') ? v.replace(UrlIn, bareUrl) : v;
  return s.length > max ? `${s.slice(0, max - 1)}…` : s;
}

/**
 * What a run was pointed at, in words (journal.md §2): the text of its `structure`, `component`, `node` and `target`
 * values, in that order, each once and cut; loci picked in the page and query objects said in a word. The values are
 * the runner's own (checked, merged over the defaults), so nothing is resolved again.
 */
export function targetsOf(values: Readonly<Record<string, unknown>>): string[] {
  const out: string[] = [];
  for (const key of TargetKeys) {
    const v = values[key];
    let text: string | undefined;
    if (typeof v === 'string') text = v.trim() ? plainText(v.trim(), TargetMax) : undefined;
    else if (v && typeof v === 'object') text = 'loci' in v ? 'atoms picked in the view' : '(a selector)';
    if (text && !out.includes(text)) out.push(text);
  }
  return out;
}

/** A result's `changed`, small (journal.md §2); undefined when it changed nothing a reader re-reads. */
export function changedOf(changed: Changed | undefined): RunChanged | undefined {
  if (!changed) return undefined;
  const nodes = changed.nodes?.length ? changed.nodes : undefined;
  const out: {
    nodes?: readonly string[];
    more?: number;
    selection?: true;
    view?: true;
    sequence?: true;
  } = {};
  if (nodes) {
    out.nodes = nodes.slice(0, ChangedNodesMax);
    if (nodes.length > ChangedNodesMax) out.more = nodes.length - ChangedNodesMax;
  }
  if (changed.selection) out.selection = true;
  if (changed.view) out.view = true;
  if (changed.sequence) out.sequence = true;
  return Object.keys(out).length ? out : undefined;
}

const envelope = (value: Json) => ({ codec: RecordCodec, version: THREAD_PAYLOAD_VERSION as 1, value });
const limits = () => createThreadLimits({ ...THREAD_LIMIT_HARD_BOUNDS, maxActiveStreams: 8 });

export interface SessionRecordOptions {
  readonly api: ViewerApi;
  /** The steps the history undoes (an `undo` entry each, naming it). */
  readonly undone?: Observable<HistoryStep>;
  /** The steps Redo puts back (a `redo` entry each). */
  readonly redone?: Observable<HistoryStep>;
  /** Where it is written (memory by default). */
  readonly rows?: ThreadRows;
  /** Where "Keep" puts it (the spec's `persist.rows`); the browser's IndexedDB when absent. */
  readonly keepIn?: ThreadRows;
  /** Its id (one per viewer session); new by default. */
  readonly id?: string;
  /** The clock of the entries and of the thread (its inactivity timeout); `Date.now` by default. */
  readonly clock?: () => number;
  /** Runs `fn` every `ms` until the returned function is called (`setInterval` by default): the turns' heartbeat. */
  readonly every?: (fn: () => void, ms: number) => () => void;
  /** Runs `fn` once after `ms`, unless the returned function is called first (`setTimeout` by default): reply deltas' coalescing. */
  readonly after?: (fn: () => void, ms: number) => () => void;
}

const everyInterval = (fn: () => void, ms: number) => {
  const handle = setInterval(fn, ms);
  return () => clearInterval(handle);
};

const afterTimeout = (fn: () => void, ms: number) => {
  const handle = setTimeout(fn, ms);
  return () => clearTimeout(handle);
};

/** The operations whose runs undo or redo steps: the actor of an undo is the origin of its run. */
const UndoOperations = new Set(['history.undo', 'history.undo-turn', 'history.redo']);

export async function createSessionRecord(o: SessionRecordOptions): Promise<SessionRecord> {
  const now = o.clock ?? Date.now;
  const clock = { now };
  const every = o.every ?? everyInterval;
  const after = o.after ?? afterTimeout;
  const id = o.id ?? `session-${now().toString(36)}-${Math.random().toString(36).slice(2, 8)}`;
  const fid = threadId(id);
  let host: ThreadHost = createThreadHost({ rows: o.rows ?? createMemoryRows(), clock });
  let kept = false;
  let disposed = false;
  let seq = 0;
  let keptIndex: { index: RecordIndex; title: string } | undefined;
  /** Lets go of the kept rows (closes their connection) once, when the kept copy is forgotten or the record ends. */
  let keptRelease: (() => void) | undefined;
  const letGo = () => {
    const release = keptRelease;
    keptRelease = undefined;
    release?.();
  };
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
    await host.create({ commandId: threadCommandId(`${id}:create`), threadId: fid, ownerId: Owner, limits: limits() }),
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
            commandId: threadCommandId(`${id}:${n}`),
            threadId: fid,
            entryId: threadEntryId(`${entry.kind}-${n}`),
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
  // person does; the thread would abort it after its inactivity timeout otherwise).
  const open = new Map<number, { readonly entryId: ThreadEntryId; readonly stop: () => Promise<void> }>();
  let stopBeat: (() => void) | undefined;
  const beat = () =>
    void quiet(
      write(async () => {
        for (const t of open.values()) {
          const r = await host.heartbeat({ threadId: fid, entryId: t.entryId });
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

  // The runs of Undo, Undo turn and Redo under way: their origin is the actor of the steps the history reports meanwhile
  // (it emits them inside the run), and an Undo turn's steps are one entry, written when the run ends.
  const undoing = new Map<
    number,
    { readonly operation: string; readonly by: RecordActor; readonly steps: HistoryStep[] }
  >();
  const actorOf = (origin: RunOrigin): RecordActor =>
    origin.kind === 'agent' ? { kind: origin.kind, turn: origin.turn } : { kind: origin.kind };
  const startedSub = o.api.events.started?.subscribe((r) => {
    if (UndoOperations.has(r.operation))
      undoing.set(r.id, { operation: r.operation, by: actorOf(r.origin), steps: [] });
  });
  const latestUndoing = () => {
    let last: ReturnType<typeof undoing.get>;
    for (const u of undoing.values()) last = u;
    return last;
  };

  const sub = o.api.events.finished.subscribe((r) => {
    const turnUndo = undoing.get(r.id);
    if (turnUndo) {
      undoing.delete(r.id);
      // An Undo turn: every step it took back, one entry.
      if (turnUndo.operation === 'history.undo-turn' && turnUndo.steps.length) {
        const steps = turnUndo.steps;
        const oldest = steps[steps.length - 1];
        const turn = oldest.origin.turn ?? (typeof r.values.turn === 'string' ? r.values.turn : undefined);
        append({
          kind: 'undo',
          step: steps[0].id,
          ...(steps.length > 1 ? { steps: steps.map((st) => st.id) } : {}),
          what: oldest.label,
          ...(turn ? { turn } : {}),
          by: turnUndo.by,
          at: now(),
        });
      }
    }
    if (r.nested || !r.result) return;
    // The definition only (not `describe`, which resolves params and availability): this runs for every hover.
    const def = o.api.definitions.value.get(r.operation);
    if (isQuiet(r, def)) return;
    const at = r.started;
    const origin = actorOf(r.origin);
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
    const targets = targetsOf(r.values);
    const changed = r.result.ok ? changedOf(r.result.changed) : undefined;
    append({
      kind: 'run',
      operation: r.operation,
      title: r.title,
      values: summarizeValues(r.values),
      origin,
      result: r.result.ok
        ? { ok: true, did: r.result.did }
        : { ok: false, code: r.result.code, message: r.result.message },
      ...(targets.length ? { targets } : {}),
      ...(changed ? { changed } : {}),
      at,
      ms: Math.round(r.ms ?? 0),
      ...(r.step !== undefined ? { step: r.step } : {}),
    });
  });

  const undoneSub = o.undone?.subscribe((step) => {
    const u = latestUndoing();
    // Inside an Undo turn: collected, one entry at the run's end.
    if (u?.operation === 'history.undo-turn') {
      u.steps.push(step);
      return;
    }
    append({
      kind: 'undo',
      step: step.id,
      what: step.label,
      ...(step.origin.turn ? { turn: step.origin.turn } : {}),
      ...(u ? { by: u.by } : {}),
      at: now(),
    });
  });

  const redoneSub = o.redone?.subscribe((step) => {
    const u = latestUndoing();
    append({
      kind: 'redo',
      step: step.id,
      what: step.label,
      ...(step.origin.turn ? { turn: step.origin.turn } : {}),
      ...(u ? { by: u.by } : {}),
      at: now(),
    });
  });

  let keeping: Promise<void> | undefined;

  const record: SessionRecord = {
    id,
    get kept() {
      return kept;
    },
    ...(o.keepIn ? { keepIn: o.keepIn } : {}),
    async read(limit, before) {
      await chain;
      return readBackward(host, fid, limit, before);
    },
    async readAfter(from, limit) {
      await chain;
      return readForward(host, fid, from, limit);
    },
    onChange(listener) {
      listeners.add(listener);
      return () => listeners.delete(listener);
    },
    beginTurn(turn, prompt, focus, model) {
      if (disposed) return { reply: () => {}, settle: async () => {} };
      const n = ++seq;
      const entryId = threadEntryId(`turn-${n}`);
      const begun: TurnEntry = {
        kind: 'turn',
        turn,
        prompt,
        ...(focus ? { focus } : {}),
        ...(model ? { model } : {}),
        at: now(),
      };
      /** The stream's next offset (from `begin`), what was streamed, and the newest text waiting to be. */
      let offset: number | undefined;
      let sent = '';
      let latest = '';
      let deltas = 0;
      let cancel: (() => void) | undefined;
      let streaming = true;
      void quiet(
        write(async () => {
          const began = ok(
            await host.begin({
              commandId: threadCommandId(`${id}:${n}`),
              threadId: fid,
              entryId,
              descriptor: envelope(begun as unknown as Json),
            }),
            'begin a turn',
          );
          const state = began.status === 'accepted' ? began.value.entry.state : undefined;
          offset = state?.kind === 'streaming' ? state.nextOffset : undefined;
          dirty = true;
        }),
      );
      /** The text since the last delta, or the whole of it when it is not an extension of what was sent. */
      const flush = () => {
        cancel = undefined;
        if (!streaming || !open.has(n)) return;
        void quiet(
          write(async () => {
            // Decided here, in the write order: the newest text against what the last delta left.
            const text = latest;
            if (!streaming || !open.has(n) || offset === undefined || text === sent) return;
            if (text.length > ReplyPreviewMax || deltas >= ReplyDeltasMax) {
              streaming = false;
              return;
            }
            const r = await host.appendDelta({
              commandId: threadCommandId(`${id}:${n}:d${++deltas}`),
              threadId: fid,
              entryId,
              start: offset,
              payload: envelope(text.startsWith(sent) ? { text: text.slice(sent.length) } : { replace: text }),
            });
            // A stream that cannot be written (kept meanwhile, or at its limits) is not previewed any more; the
            // final text still arrives with `settle`.
            if (r.status !== 'accepted' || r.value.kind !== 'appended') {
              streaming = false;
              return;
            }
            offset = r.value.nextOffset;
            sent = text;
          }),
        );
      };
      const reply: TurnWriter['reply'] = (text) => {
        if (!streaming || !open.has(n)) return;
        latest = text;
        cancel ??= after(flush, ReplyMs);
      };
      const settle: TurnWriter['settle'] = (end) => {
        if (!open.delete(n)) return chain;
        streaming = false;
        cancel?.();
        cancel = undefined;
        beating();
        // Settled with nothing to say (stopped while parked, or at the end of the page): what had streamed stays.
        const said = end.reply || latest;
        return quiet(
          write(async () => {
            const done: TurnEntry = {
              ...begun,
              reply: said,
              outcome: end.outcome,
              ...(end.model ? { model: end.model } : {}),
              ...(end.usage ? { usage: end.usage } : {}),
            };
            ok(
              await host.settle({
                commandId: threadCommandId(`${id}:${n}:settle`),
                threadId: fid,
                entryId,
                payload: envelope(done as unknown as Json),
              }),
              'settle a turn',
            );
            forgetPreview(host, entryId);
            dirty = true;
            await touch(true);
          }),
        );
      };
      open.set(n, { entryId, stop: () => settle({ reply: '', outcome: 'stopped' }) });
      beating();
      return { reply, settle };
    },
    keep(rows, index, title, release) {
      if (kept || disposed) release?.();
      if (kept) return Promise.resolve();
      if (disposed) return Promise.reject(new Error('The record has ended.'));
      keeping ??= write(async () => {
        if (kept) return release?.();
        // Copied in order, the same entry ids and positions; a turn under way is copied as begun and settles there.
        const target = createThreadHost({ rows, clock });
        try {
          ok(
            await target.create({
              commandId: threadCommandId(`${id}:create`),
              threadId: fid,
              ownerId: Owner,
              limits: limits(),
            }),
            'be kept',
          );
          await copyThread(host, target, fid, id, ok);
        } catch (e) {
          // Nothing half-kept: the copy goes, and the record goes on in memory.
          await target.forget(fid).catch(() => undefined);
          release?.();
          throw e;
        }
        host = target;
        kept = true;
        keptIndex = { index, title };
        keptRelease = release;
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
        const target = createThreadHost({ rows: createMemoryRows(), clock });
        ok(
          await target.create({
            commandId: threadCommandId(`${id}:create`),
            threadId: fid,
            ownerId: Owner,
            limits: limits(),
          }),
          'go on in memory',
        );
        await copyThread(host, target, fid, id, ok);
        const from = host;
        const { index } = keptIndex;
        host = target;
        kept = false;
        keptIndex = undefined;
        dirty = false;
        await from.forget(fid);
        await index.forget(id);
        letGo();
      });
    },
    flush: () => chain,
    async dispose() {
      if (!disposed) {
        disposed = true;
        sub.unsubscribe();
        startedSub?.unsubscribe();
        undoneSub?.unsubscribe();
        redoneSub?.unsubscribe();
        // A turn still open ends as stopped (journal.md §4), so it is neither "under way" forever nor left to abort.
        for (const t of [...open.values()]) void t.stop();
        stopBeat?.();
        stopBeat = undefined;
        void quiet(write(() => touch(true), false));
      }
      await chain;
      letGo();
      listeners.clear();
    },
  };
  return record;
}

/**
 * An entry as History reads it (settled or under way). An aborted turn (its turn never settled: the page went away while
 * it ran, the thread's inactivity timeout) reads as stopped, from what it was begun with; other aborted entries are left
 * out (only turns stream).
 */
/** Entries newest first: `limit` of them before `before`, in backward windows of at most the thread's limit. */
async function readBackward(host: ThreadHost, fid: ThreadId, limit: number, before?: number): Promise<RecordPage> {
  const items: RecordItem[] = [];
  let start = before;
  let end: number | undefined;
  let more = false;
  for (let left = Math.max(1, limit); left > 0; ) {
    const page = await host.window({
      threadId: fid,
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
      const item = await itemOf(host, fid, page.entries[i]);
      if (item) items.push(item);
    }
    more = page.hasMore;
    start = page.start;
    left -= page.entries.length;
    if (!page.hasMore || !page.entries.length) break;
  }
  return { items, more, start: start ?? 0, end: end ?? start ?? 0 };
}

/** Entries at `from` and after, newest first: the oldest `limit` of them, in forward windows. */
async function readForward(host: ThreadHost, fid: ThreadId, from: number, limit: number): Promise<RecordPage> {
  const entries: ThreadEntry[] = [];
  let start = from;
  let more = false;
  for (let left = Math.max(1, limit); left > 0; ) {
    const page = await host.window({
      threadId: fid,
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
    const item = await itemOf(host, fid, entries[i]);
    if (item) items.push(item);
  }
  return { items, more, start: from, end: start };
}

/**
 * What History reads of a record (journal.md §3): the session's own, or a kept one opened after a reload (Home's list),
 * read-only.
 */
export interface RecordReader {
  /** It outlives the page. */
  readonly kept: boolean;
  read(limit: number, before?: number): Promise<RecordPage>;
  readAfter(from: number, limit: number): Promise<RecordPage>;
  onChange(listener: () => void): () => void;
}

/** A kept record, read-only: what Home opens (journal.md §3, "After a reload"); nothing in it is undone or replayed. */
export function openKeptRecord(rows: ThreadRows, id: string): RecordReader {
  const host = createThreadHost({ rows });
  const fid = threadId(id);
  return {
    kept: true,
    read: (limit, before) => readBackward(host, fid, limit, before),
    readAfter: (from, limit) => readForward(host, fid, from, limit),
    onChange: () => () => {},
  };
}

/**
 * The record as a session file keeps it (journal.md §5, `record.jsonl`): an entry a line, oldest first, each in its
 * envelope (`{ codec, version, value }`).
 */
export async function recordLines(reader: RecordReader): Promise<string> {
  const { items } = await reader.readAfter(0, Number.MAX_SAFE_INTEGER);
  const lines: string[] = [];
  for (let i = items.length - 1; i >= 0; i--) lines.push(JSON.stringify(envelope(items[i].entry as unknown as Json)));
  return lines.join('\n');
}

/** A record read from a session file: read-only, as a kept one (journal.md §5), its span and size. */
export interface FileRecord {
  readonly reader: RecordReader;
  readonly count: number;
  readonly first: number;
  readonly last: number;
}

/** A session file's record from its lines; a line that is not an entry of ours is skipped. None: no entries. */
export function recordFromLines(text: string): FileRecord | undefined {
  const entries: RecordEntry[] = [];
  for (const line of text.split(/\r?\n/)) {
    if (!line.trim()) continue;
    try {
      const e = JSON.parse(line) as { codec?: string; value?: RecordEntry };
      const v = e?.codec === RecordCodec ? e.value : undefined;
      if (v && typeof v === 'object' && 'kind' in v && 'at' in v) entries.push(v);
    } catch {
      // not a line of ours
    }
  }
  if (!entries.length) return undefined;
  const item = (i: number): RecordItem => ({ id: `file-${i}`, position: i, entry: entries[i], open: false });
  const n = entries.length;
  const reader: RecordReader = {
    kept: true,
    async read(limit, before) {
      const end = Math.min(before ?? n, n);
      const start = Math.max(0, end - Math.max(1, limit));
      const items: RecordItem[] = [];
      for (let i = end - 1; i >= start; i--) items.push(item(i));
      return { items, more: start > 0, start, end };
    },
    async readAfter(from, limit) {
      const end = Math.min(n, from + Math.max(1, limit));
      const items: RecordItem[] = [];
      for (let i = end - 1; i >= from; i--) items.push(item(i));
      return { items, more: end < n, start: from, end };
    },
    onChange: () => () => {},
  };
  return { reader, count: n, first: entries[0].at, last: entries[n - 1].at };
}

async function itemOf(host: ThreadHost, fid: ThreadId, e: ThreadEntry): Promise<RecordItem | undefined> {
  if (e.state.kind !== 'streaming') forgetPreview(host, e.entryId);
  if (e.state.kind === 'aborted') {
    const begun = e.descriptor.value as TurnEntry | undefined;
    if (!begun || typeof begun !== 'object' || begun.kind !== 'turn') return undefined;
    return { id: e.entryId, position: e.position, entry: { ...begun, outcome: 'stopped' }, open: false };
  }
  const value = (e.state.kind === 'settled' ? e.state.payload : e.descriptor).value as unknown as RecordEntry;
  if (!value || typeof value !== 'object' || !('kind' in value) || !('at' in value)) return undefined;
  if (e.state.kind === 'streaming' && value.kind === 'turn') {
    // A turn under way: the reply as it streams (the final text replaces it when the turn settles).
    const said = await previewOf(host, fid, e.entryId);
    return { id: e.entryId, position: e.position, entry: said ? { ...value, reply: said } : value, open: true };
  }
  return { id: e.entryId, position: e.position, entry: value, open: e.state.kind !== 'settled' };
}

/**
 * The reply streamed so far per stream, read incrementally: a read of a turn under way costs what was added since.
 * `tail` is the stream's reads in order (journal.md §2): two readers at once (History and the assistant, an export
 * racing a refresh) would both fold the deltas after the same offset, and the reply would double.
 */
type Preview = { offset: number; text: string; tail: Promise<string> };
const previews = new WeakMap<ThreadHost, Map<string, Preview>>();

/** @internal For tests: how many streams' previews are held for a host. */
export const previewsHeld = (host: ThreadHost): number => previews.get(host)?.size ?? 0;

/** A stream's preview is let go once its turn is not under way (settled, aborted, or gone with the host). */
const forgetPreview = (host: ThreadHost, entryId: string) => void previews.get(host)?.delete(entryId);

/** A turn's streamed reply: the deltas folded (`{ text }` appends, `{ replace }` starts over). Empty when there is none. */
function previewOf(host: ThreadHost, fid: ThreadId, entryId: ThreadEntryId): Promise<string> {
  let cache = previews.get(host);
  if (!cache) previews.set(host, (cache = new Map()));
  let seen = cache.get(entryId);
  if (!seen) cache.set(entryId, (seen = { offset: 0, text: '', tail: Promise.resolve('') }));
  const s = seen;
  const fold = async () => {
    try {
      for (;;) {
        const page = await host.deltas({ threadId: fid, entryId, afterOffset: s.offset, limit: 256 });
        if ('status' in page) break;
        for (const b of page.batches) {
          const v = b.payload.value as { text?: unknown; replace?: unknown } | null;
          if (typeof v?.replace === 'string') s.text = v.replace;
          else if (typeof v?.text === 'string') s.text += v.text;
        }
        s.offset = page.nextOffset;
        if (!page.hasMore || !page.batches.length) break;
      }
    } catch {
      // a stream that cannot be read shows nothing more
    }
    return s.text;
  };
  s.tail = s.tail.then(fold, fold);
  return s.tail;
}

/** Every entry of a thread into another, in order, at the same positions (keep, and unkeep back into memory). */
async function copyThread(
  from: ThreadHost,
  to: ThreadHost,
  fid: ThreadId,
  id: string,
  ok: <T extends { status: string }>(r: T, what: string) => T,
): Promise<void> {
  let start = 0;
  for (;;) {
    const page = await from.window({ threadId: fid, start, limit: RecordWindowMax, direction: 'forward' });
    if ('status' in page) throw new Error(`The record could not be read (${page.reason}).`);
    for (const e of page.entries) await copyEntry(to, fid, e, id, ok);
    if (!page.hasMore) break;
    start = page.endExclusive;
  }
}

/** One entry into the kept record, at the same position: settled as settled, under way as begun, aborted as aborted. */
async function copyEntry(
  target: ThreadHost,
  fid: ThreadId,
  e: ThreadEntry,
  id: string,
  ok: <T extends { status: string }>(r: T, what: string) => T,
): Promise<void> {
  const commandId = threadCommandId(`${id}:copy:${e.entryId}`);
  if (e.state.kind === 'settled') {
    ok(
      await target.append({
        commandId,
        threadId: fid,
        entryId: e.entryId,
        descriptor: e.descriptor,
        payload: e.state.payload,
      }),
      'be kept',
    );
    return;
  }
  ok(await target.begin({ commandId, threadId: fid, entryId: e.entryId, descriptor: e.descriptor }), 'be kept');
  if (e.state.kind === 'aborted') {
    ok(
      await target.abort({
        commandId: threadCommandId(`${id}:copy:${e.entryId}:abort`),
        threadId: fid,
        entryId: e.entryId,
        reason: e.state.reason,
      }),
      'be kept',
    );
  }
}

/** The browser's kept records' rows, only when the database exists (listing never makes an empty one). */
export async function openRecordsRows(
  database: string = RecordsDatabase,
): Promise<(ThreadRows & { close?(): void }) | undefined> {
  let idb: IDBFactory | undefined;
  try {
    idb = globalThis.indexedDB;
  } catch {
    return undefined; // blocked site data
  }
  if (!idb) return undefined;
  if (typeof idb.databases === 'function') {
    const names = (await idb.databases()).map((d) => d.name);
    if (!names.includes(database)) return undefined;
  }
  // Loaded when asked, as History's Keep does: a viewer that never lists them never pays for it.
  const { createIndexedDbRows } = await import('./indexeddb-rows.js');
  return createIndexedDbRows<ThreadTable>(database);
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
  /** A fold of the index under way is done (journal.md §1): the index folds itself in the background. */
  idle(): Promise<void>;
}

/**
 * The index thread is folded when it holds more than this many entries (journal.md §1): every keep, touch and forget is
 * an entry, so it would grow with the time records are appended to. Folded, it holds one entry per kept record.
 */
export const IndexFoldAt = 64;

/** One index entry as the index folds them: a record's newest title and last time, and when it was first listed. */
interface IndexValue {
  readonly id: string;
  readonly title?: string;
  readonly at?: number;
  readonly first?: number;
  readonly forget?: true;
}

/**
 * The index as a thread of its own in the same rows: append-only changes, folded when listed. Past `IndexFoldAt` entries
 * it is rewritten as one entry per kept record (its newest title and time, and the time it was first listed), which
 * drops the superseded touches and the forgotten records' tombstones. A record is never dropped by folding (a kept record
 * idle for 90 days is only marked in Home's list, and the person forgets it, journal.md §7).
 */
export async function createRecordIndex(rows: ThreadRows): Promise<RecordIndex> {
  const host = createThreadHost({ rows });
  const fid = threadId('prism-viewer-records');
  const create = async () => {
    const created = await host.create({
      commandId: threadCommandId('records:create'),
      threadId: fid,
      ownerId: Owner,
      limits: limits(),
    });
    if (created.status !== 'accepted' && created.status !== 'conflict')
      throw new Error('The records index could not start.');
  };
  await create();
  let n = 0;
  const unique = () => `${Date.now().toString(36)}-${++n}-${Math.random().toString(36).slice(2, 8)}`;
  const append = async (value: Json) =>
    host.append({
      commandId: threadCommandId(`records:${unique()}`),
      threadId: fid,
      entryId: threadEntryId(`r-${unique()}`),
      descriptor: envelope({ kind: 'index' }),
      payload: envelope(value),
    });

  /** The entries folded: one per record, oldest listed first; and how many entries there were. */
  const read = async () => {
    const out = new Map<string, KeptRecord>();
    let count = 0;
    let start = 0;
    for (;;) {
      const page = await host.window({ threadId: fid, start, limit: RecordWindowMax, direction: 'forward' });
      if ('status' in page) break;
      for (const e of page.entries) {
        if (e.state.kind !== 'settled') continue;
        count++;
        const v = e.state.payload.value as IndexValue;
        if (v.forget) out.delete(v.id);
        else {
          const had = out.get(v.id);
          out.set(v.id, {
            id: v.id,
            title: v.title ?? had?.title ?? v.id,
            first: had?.first ?? v.first ?? v.at ?? 0,
            last: v.at ?? had?.last ?? 0,
          });
        }
      }
      if (!page.hasMore) break;
      start = page.endExclusive;
    }
    return { records: [...out.values()], count };
  };

  // One thing at a time in this index (a fold rewrites the thread: nothing is appended or read meanwhile). Other pages
  // writing to the same rows are not fenced; a touch they lose is the next one's to write.
  let tail: Promise<unknown> = Promise.resolve();
  const exclusive = <T>(work: () => Promise<T>): Promise<T> => {
    const run = tail.then(work, work);
    tail = run.catch(() => undefined);
    return run;
  };

  /** How many records the index held when last read: folding is for when the entries are far more than that. */
  let known = 0;
  const foldAt = () => Math.max(IndexFoldAt, 2 * known);
  let folding: Promise<void> | undefined;

  /** Appends a change; on a thread folded away by another page, once more on the new one. */
  const put = async (value: IndexValue): Promise<number> => {
    let r = await append(value as unknown as Json);
    if (r.status === 'not-found') {
      await create();
      r = await append(value as unknown as Json);
    }
    if (r.status !== 'accepted') throw new Error(`The records index could not be written (${r.status}).`);
    return r.value.thread.nextPosition;
  };

  /** Rewrites the thread as one entry per record, queued behind what is under way; once at a time. */
  const fold = () =>
    (folding ??= exclusive(async () => {
      const { records, count } = await read();
      known = records.length;
      if (count <= foldAt()) return;
      await host.forget(fid);
      await create();
      for (const r of records) await put({ id: r.id, title: r.title, at: r.last, first: r.first });
    })
      .catch((e) => console.error('[prism-viewer] the records index could not be folded:', e))
      .finally(() => {
        folding = undefined;
      }));

  const change = (value: IndexValue) =>
    exclusive(async () => {
      const size = await put(value);
      if (size > foldAt()) void fold();
    });

  return {
    list: () =>
      exclusive(async () => {
        const { records, count } = await read();
        known = records.length;
        if (count > foldAt()) void fold();
        return records.sort((a, b) => b.last - a.last);
      }),
    touch: (id, title, at) => change({ id, title, at }),
    forget: (id) => change({ id, forget: true }),
    idle: async () => {
      await folding;
    },
  };
}

/**
 * Forgets a kept record (the Storage view's Forget, journal.md §7): its rows, then its index entry. Not the current
 * session's record while it is kept: that one is `unkeep`'s.
 */
export async function forgetKept(rows: ThreadRows, index: RecordIndex, id: string): Promise<void> {
  await createThreadHost({ rows }).forget(threadId(id));
  await index.forget(id);
}

/**
 * A kept record's size, about: the characters of its rows as JSON (the threads' default scope, `thread:<id>`).
 * Undefined when its rows are not there.
 */
export async function keptSize(rows: ThreadRows, id: string): Promise<number | undefined> {
  const loaded = await rows.load(`thread:${id}`);
  if (!loaded) return undefined;
  let size = 0;
  for (const r of loaded.rows) size += r.key.length + JSON.stringify(r.row).length;
  return size;
}
