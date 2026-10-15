/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { observableState, type ReadonlyState } from '@molstar/prism-ui/core';
import type { RunOriginKind } from '../api/area.js';
import type { RecordEntry, RecordItem, RunEntry, SessionRecord, TurnEntry } from './record.js';

// What History shows of the session's record (journal.md §3), DOM-free. The record is read in windows, newest first (the
// newest `Window` entries, and "Show earlier" reads the window before the oldest loaded). It is read again after a
// change, once per frame, and only while the tab is active: a hidden tab holds no subscription and reads nothing.
// The rows: a person's (script's, link's) run is one row; an agent's turn is one row holding its runs; an `undo` entry is
// no row, it marks the step it undid; a `session` entry is a small separator.

/** Entries read at first, and by each "Show earlier". */
export const HistoryWindow = 200;

export interface HistoryRunRow {
  readonly kind: 'run';
  readonly key: string;
  readonly origin: RunOriginKind;
  /** What it did, or why it failed. */
  readonly text: string;
  readonly failed: boolean;
  readonly undone: boolean;
  readonly at: number;
}

export interface HistoryTurnRow {
  readonly kind: 'turn';
  readonly key: string;
  readonly turn: string;
  readonly prompt: string;
  /** The reply, whole; the view shows its first line and expands it. */
  readonly reply?: string;
  readonly outcome?: TurnEntry['outcome'];
  /** Under way: begun, not settled. */
  readonly open: boolean;
  /** Every run of it that changed something was undone. */
  readonly undone: boolean;
  readonly runs: readonly HistoryRunRow[];
  readonly at: number;
}

export interface HistorySessionRow {
  readonly kind: 'session';
  readonly key: string;
  readonly text: string;
  readonly at: number;
}

export type HistoryRow = HistoryRunRow | HistoryTurnRow | HistorySessionRow;

export interface HistoryModel {
  /** Newest first. */
  readonly rows: ReadonlyState<readonly HistoryRow[]>;
  /** Earlier entries are in the record than the rows hold. */
  readonly more: ReadonlyState<boolean>;
  /** The record outlives the page. */
  readonly kept: ReadonlyState<boolean>;
  /** Active: reads the record now, and again after each change. Inactive: nothing is read or subscribed. */
  setActive(active: boolean): void;
  /** Reads the window before the oldest loaded one. */
  showEarlier(): Promise<void>;
  /** Resolves when no read is waiting or running. */
  idle(): Promise<void>;
  dispose(): void;
}

export interface HistoryModelOptions {
  /** Entries per window (200). */
  readonly window?: number;
  /** Runs `fn` once, soon: a frame (the default where there is one) or a microtask. */
  readonly schedule?: (fn: () => void) => void;
}

const nextFrame = (fn: () => void) => {
  if (typeof requestAnimationFrame === 'function') requestAnimationFrame(() => fn());
  else queueMicrotask(fn);
};

export function createHistoryModel(record: SessionRecord, o: HistoryModelOptions = {}): HistoryModel {
  const size = o.window ?? HistoryWindow;
  const schedule = o.schedule ?? nextFrame;
  const rows = observableState<readonly HistoryRow[]>([]);
  const more = observableState(false);
  const kept = observableState(record.kept);

  /** The loaded entries, by position. */
  let loaded: RecordItem[] = [];
  let hasMore = false;
  let active = false;
  let disposed = false;
  let unlisten: (() => void) | undefined;
  /** A read is wanted / running; `chain` ends when none is. */
  let wanted = false;
  let scheduled = false;
  let chain: Promise<void> = Promise.resolve();
  let signatures = new Map<string, string>();
  let previous = new Map<string, HistoryRow>();

  const publish = () => {
    // Rows whose content did not change keep their object, so a keyed list leaves them alone.
    const next = rowsOf(loaded);
    const sigs = new Map<string, string>();
    const prev = new Map<string, HistoryRow>();
    const out = next.map((row) => {
      const sig = JSON.stringify(row);
      sigs.set(row.key, sig);
      const same = signatures.get(row.key) === sig ? previous.get(row.key) : undefined;
      const use = same ?? row;
      prev.set(row.key, use);
      return use;
    });
    signatures = sigs;
    previous = prev;
    rows.set(out);
    more.set(hasMore);
    kept.set(record.kept);
  };

  /** Reads the newest window and merges it into what is loaded. */
  const refresh = async () => {
    const fresh = await record.read(size);
    if (disposed || !active) return;
    if (!fresh.more || !loaded.length) {
      // The window holds the whole record, or nothing is loaded yet.
      loaded = fresh.items;
      hasMore = fresh.more;
    } else {
      const oldest = fresh.items.length ? fresh.items[fresh.items.length - 1].position : Number.POSITIVE_INFINITY;
      if (oldest > loaded[0].position + 1) {
        // More new entries than a window: start over from the newest; earlier ones come by "Show earlier".
        loaded = fresh.items;
        hasMore = true;
      } else if (loaded.some((i) => i.open && i.position < oldest)) {
        // A turn under way, older than this window, may have settled since: read all that is loaded again.
        const all = await record.read(Math.max(size, loaded.length));
        if (disposed || !active) return;
        loaded = all.items;
        hasMore = all.more;
      } else {
        const older = loaded.filter((i) => i.position < oldest);
        loaded = [...fresh.items, ...older];
        if (!older.length) hasMore = fresh.more;
      }
    }
    publish();
  };

  const run = () => {
    chain = chain
      .then(async () => {
        while (wanted && active && !disposed) {
          wanted = false;
          await refresh();
        }
      })
      .catch((e) => console.error('[prism-viewer] history:', e));
  };

  const want = () => {
    wanted = true;
    if (scheduled || !active || disposed) return;
    scheduled = true;
    // One read per frame, however many writes came.
    schedule(() => {
      scheduled = false;
      if (wanted && active && !disposed) run();
    });
  };

  return {
    rows: rows.asReadonly,
    more: more.asReadonly,
    kept: kept.asReadonly,
    setActive(next) {
      if (next === active || disposed) return;
      active = next;
      if (active) {
        unlisten = record.onChange(want);
        wanted = true;
        run();
      } else {
        unlisten?.();
        unlisten = undefined;
        wanted = false;
      }
    },
    async showEarlier() {
      await chain;
      if (!active || disposed || !hasMore || !loaded.length) return;
      const before = loaded[loaded.length - 1].position;
      const page = await record.read(size, before);
      if (!active || disposed) return;
      loaded = [...loaded, ...page.items];
      hasMore = page.more && page.items.length > 0;
      publish();
    },
    async idle() {
      // A scheduled read has not started; wait for it, then for the chain.
      for (let i = 0; i < 100 && (scheduled || wanted); i++) {
        await chain;
        await new Promise<void>((r) => setTimeout(r, 0));
      }
      await chain;
    },
    dispose() {
      disposed = true;
      unlisten?.();
      unlisten = undefined;
    },
  };
}

const runRow = (item: RecordItem, e: RunEntry, undone: ReadonlySet<number>): HistoryRunRow => ({
  kind: 'run',
  key: item.id,
  origin: e.origin.kind,
  text: e.result.ok ? e.result.did : e.result.message,
  failed: !e.result.ok,
  undone: e.step !== undefined && undone.has(e.step),
  at: e.at,
});

/** The rows of the loaded entries (newest first in, newest first out). */
export function rowsOf(items: readonly RecordItem[]): HistoryRow[] {
  const undone = new Set<number>();
  for (const i of items) if (i.entry.kind === 'undo') undone.add(i.entry.step);

  interface Turn {
    readonly item: RecordItem;
    readonly entry: TurnEntry;
    readonly slot: number;
    readonly runs: RunEntry[];
    readonly rows: HistoryRunRow[];
  }
  // Oldest first: a run of an agent's turn joins the latest turn of that id begun before it (without one, it is a row).
  const out: (HistoryRow | undefined)[] = [];
  const turns: Turn[] = [];
  const latest = new Map<string, Turn>();
  for (let n = items.length - 1; n >= 0; n--) {
    const item = items[n];
    const e: RecordEntry = item.entry;
    if (e.kind === 'turn') {
      const t: Turn = { item, entry: e, slot: out.length, runs: [], rows: [] };
      turns.push(t);
      latest.set(e.turn, t);
      out.push(undefined);
    } else if (e.kind === 'run') {
      const t = e.origin.kind === 'agent' && e.origin.turn !== undefined ? latest.get(e.origin.turn) : undefined;
      if (t) {
        t.runs.push(e);
        t.rows.push(runRow(item, e, undone));
      } else out.push(runRow(item, e, undone));
    } else if (e.kind === 'session') {
      out.push({ kind: 'session', key: item.id, text: sessionText(e.event, e.detail), at: e.at });
    }
  }
  for (const t of turns) {
    const changed = t.runs.filter((r) => r.step !== undefined);
    out[t.slot] = {
      kind: 'turn',
      key: t.item.id,
      turn: t.entry.turn,
      prompt: t.entry.prompt,
      ...(t.entry.reply !== undefined ? { reply: t.entry.reply } : {}),
      ...(t.entry.outcome ? { outcome: t.entry.outcome } : {}),
      open: t.item.open,
      undone: changed.length > 0 && changed.every((r) => undone.has(r.step as number)),
      runs: t.rows,
      at: t.entry.at,
    };
  }
  return (out as HistoryRow[]).reverse();
}

function sessionText(event: 'opened' | 'cleared' | 'saved', detail: string | undefined): string {
  if (detail) return detail;
  return event === 'opened' ? 'Opened a session' : event === 'cleared' ? 'Cleared the session' : 'Saved the session';
}
