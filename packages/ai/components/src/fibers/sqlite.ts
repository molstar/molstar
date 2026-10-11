/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { DatabaseSync, type SQLInputValue } from 'node:sqlite';
import { type Row, RowConflict, type RowStore } from '@molstar/ai-fibers';

const SCHEMA = `
  create table if not exists scopes (scope text primary key, generation integer not null);
  create table if not exists rows (scope text, tbl text, key text, row text not null, primary key (scope, tbl, key));
`;

/**
 * The Node backend (analysis/viewer-mode.md §6): the fibers' rows in one SQLite file, through Node's built-in
 * `node:sqlite`, as the harness's (`harness/storage/sqlite`). A write is one `BEGIN IMMEDIATE` transaction that
 * compares the scope's generation, writes every row and moves the generation, so it is atomic and several
 * processes may share the file (WAL); one that wrote first makes this write a `RowConflict`, with nothing written.
 * A row is stored as JSON. Another process's write cannot be heard here, so there is no `watch`: a host sees it at
 * its next command, by the generation.
 */
export function createSqliteRows<T extends string = string>(path: string): RowStore<T> & { close(): void } {
  const db = new DatabaseSync(path);
  db.exec('pragma journal_mode = wal; pragma synchronous = normal; pragma busy_timeout = 5000;');
  db.exec(SCHEMA);
  const all = (sql: string, ...params: SQLInputValue[]) => db.prepare(sql).all(...params) as Record<string, unknown>[];
  const one = (sql: string, ...params: SQLInputValue[]) =>
    db.prepare(sql).get(...params) as Record<string, unknown> | undefined;
  const run = (sql: string, ...params: SQLInputValue[]) => db.prepare(sql).run(...params);
  const generationOf = (scope: string) =>
    (one('select generation from scopes where scope = ?', scope)?.generation as number | undefined) ?? 0;

  return {
    async load(scope) {
      // One read transaction, so the generation and the rows are of the same moment.
      db.exec('begin');
      try {
        const generation = one('select generation from scopes where scope = ?', scope)?.generation as
          | number
          | undefined;
        if (generation === undefined) return undefined;
        const rows = all('select tbl, key, row from rows where scope = ?', scope).map((r) => ({
          table: r.tbl as T,
          key: r.key as string,
          row: JSON.parse(r.row as string) as Row,
        }));
        return { generation, rows };
      } finally {
        db.exec('commit');
      }
    },
    async apply(scope, expected, writes) {
      db.exec('begin immediate');
      try {
        if (generationOf(scope) !== expected) {
          db.exec('rollback');
          throw new RowConflict(scope);
        }
        for (const { table, key, row } of writes) {
          if (row)
            run(
              'insert or replace into rows (scope, tbl, key, row) values (?, ?, ?, ?)',
              scope,
              table,
              key,
              JSON.stringify(row),
            );
          else run('delete from rows where scope = ? and tbl = ? and key = ?', scope, table, key);
        }
        run(
          'insert into scopes (scope, generation) values (?, ?) on conflict (scope) do update set generation = excluded.generation',
          scope,
          expected + 1,
        );
        db.exec('commit');
      } catch (error) {
        if (db.isTransaction) db.exec('rollback');
        throw error;
      }
    },
    async generation(scope) {
      return generationOf(scope);
    },
    async drop(scope) {
      db.exec('begin immediate');
      try {
        run('delete from rows where scope = ?', scope);
        run('delete from scopes where scope = ?', scope);
        db.exec('commit');
      } catch (error) {
        if (db.isTransaction) db.exec('rollback');
        throw error;
      }
    },
    close: () => db.close(),
  };
}
