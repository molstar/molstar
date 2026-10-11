/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { DatabaseSync } from 'node:sqlite';
import type { Json } from '@molstar/ai-agent';
import type { HarnessEvent, HarnessStorage, StepEvent } from '@molstar/ai-harness';
import { createDurableStorage, type DurableBackend } from './durable.js';

const SCHEMA = `
  create table if not exists meta (h text primary key, version integer not null, tombstone integer not null default 0);
  create table if not exists records (h text, kind text, id text, value text not null, primary key (h, kind, id));
  create table if not exists events (h text, seq integer, step text, turn text, key text, ended integer,
    value text not null, primary key (h, seq));
  create index if not exists events_step on events (h, step);
  create index if not exists events_turn on events (h, turn, seq);
  create index if not exists events_key on events (h, key);
  create index if not exists events_ended on events (h, ended, seq);
  create table if not exists compactions (h text, pointer text, seq integer, primary key (h, pointer, seq));
  create table if not exists blobs (h text, ref text, body text not null, primary key (h, ref));
`;

type Row = Record<string, unknown>;
const events = (rows: readonly Row[]) => rows.map((row) => JSON.parse(row.value as string) as HarnessEvent);

/**
 * The server's backend: one SQLite file (harness.md §14), through Node's built-in `node:sqlite`. Writes
 * apply in one `BEGIN IMMEDIATE` transaction, so several processes may share the file.
 */
export function createSqliteBackend(path: string): DurableBackend & { close(): void } {
  const db = new DatabaseSync(path);
  db.exec('pragma journal_mode = wal; pragma synchronous = normal;');
  db.exec(SCHEMA);
  const all = (sql: string, ...params: (string | number | null)[]) => db.prepare(sql).all(...params) as Row[];
  const one = (sql: string, ...params: (string | number | null)[]) => db.prepare(sql).get(...params) as Row | undefined;
  const run = (sql: string, ...params: (string | number | null)[]) => db.prepare(sql).run(...params);

  return {
    async get(h, kind, ids) {
      return ids.map((id) => {
        const row = one('select value from records where h = ? and kind = ? and id = ?', h, kind, id);
        return row === undefined ? undefined : JSON.parse(row.value as string);
      });
    },
    after: async (h, seq, limit) =>
      events(all('select value from events where h = ? and seq > ? order by seq limit ?', h, seq, limit)),
    byStep: async (h, ids) =>
      events(
        ids.length === 0
          ? []
          : all(
              `select value from events where h = ? and step in (${ids.map(() => '?').join(',')}) order by seq`,
              h,
              ...ids,
            ),
      ) as StepEvent[],
    byTurn: async (h, turn, after, limit) =>
      events(
        all(
          'select value from events where h = ? and turn = ? and seq > ? order by seq limit ?',
          h,
          turn,
          after ?? -1,
          limit,
        ),
      ),
    endedTurns: async (h, before, limit) =>
      events(
        all(
          'select value from events where h = ? and ended = 1 and seq < ? order by seq desc limit ?',
          h,
          before ?? Number.MAX_SAFE_INTEGER,
          limit,
        ),
      ),
    byKey: async (h, key) => events(all('select value from events where h = ? and key = ? order by seq', h, key)),
    async lastCompaction(h, pointer) {
      const row = one(
        'select e.value from compactions c join events e on e.h = c.h and e.seq = c.seq where c.h = ? and c.pointer = ? order by c.seq desc limit 1',
        h,
        pointer,
      );
      return row === undefined ? undefined : (JSON.parse(row.value as string) as StepEvent);
    },
    tombstone: async (h) => one('select tombstone from meta where h = ?', h)?.tombstone === 1,
    version: async (h) => (one('select version from meta where h = ?', h)?.version as number | undefined) ?? 0,
    async apply(h, expected, writes) {
      db.exec('begin immediate');
      try {
        const version = (one('select version from meta where h = ?', h)?.version as number | undefined) ?? 0;
        if (version !== expected) {
          db.exec('rollback');
          return false;
        }
        run(
          'insert into meta (h, version, tombstone) values (?, ?, ?) on conflict (h) do update set version = excluded.version, tombstone = max(tombstone, excluded.tombstone)',
          h,
          version + 1,
          writes.tombstone ? 1 : 0,
        );
        for (const { kind, id, value } of writes.records)
          run(
            'insert or replace into records (h, kind, id, value) values (?, ?, ?, ?)',
            h,
            kind,
            id,
            JSON.stringify(value),
          );
        for (const event of writes.events) {
          const step = event.kind === 'step' ? event : undefined;
          const key = event.kind === 'step' ? event.commitKey : (event.key ?? null);
          const ended = event.kind === 'step' ? event.end !== undefined : event.summary !== undefined;
          run(
            'insert into events (h, seq, step, turn, key, ended, value) values (?, ?, ?, ?, ?, ?, ?)',
            h,
            event.seq,
            step?.step ?? null,
            event.turnId ?? null,
            key,
            ended ? 1 : 0,
            JSON.stringify(event),
          );
          for (const pointer of new Set((step?.compaction ?? []).map((c) => c.pointer)))
            run('insert into compactions (h, pointer, seq) values (?, ?, ?)', h, pointer, event.seq);
        }
        db.exec('commit');
        return true;
      } catch (error) {
        db.exec('rollback');
        throw error;
      }
    },
    async putBlobs(h, blobs) {
      for (const { ref, body } of blobs)
        run('insert or ignore into blobs (h, ref, body) values (?, ?, ?)', h, ref, JSON.stringify(body));
    },
    getBlobs: async (h, refs) =>
      refs.map((ref) => {
        const row = one('select body from blobs where h = ? and ref = ?', h, ref);
        return row === undefined ? undefined : (JSON.parse(row.body as string) as Json);
      }),
    async drop(h) {
      for (const table of ['records', 'events', 'compactions', 'blobs']) run(`delete from ${table} where h = ?`, h);
    },
    close: () => db.close(),
  };
}

/** A durable storage in one SQLite file. */
export function createSqliteStorage(path: string): HarnessStorage & { close(): void } {
  const backend = createSqliteBackend(path);
  return Object.assign(createDurableStorage(backend), { close: () => backend.close() });
}
