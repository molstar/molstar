/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { digest, hashMessage, type Json, type Message } from '@molstar/ai-agent';
import {
  type HarnessClient,
  HarnessError,
  type HarnessHost,
  type HarnessReader,
  type RunnerPort,
} from '../contracts/ports.js';
import {
  type HarnessEvent,
  type HarnessId,
  type Head,
  MAIN,
  SCHEMA_VERSION,
  type TurnSummary,
} from '../contracts/records.js';
import type { Clock } from '../contracts/runner.js';
import type { HarnessStorage, StorageTx } from '../contracts/storage.js';
import { emptyHead } from '../head/index.js';
import { clientOperations } from './client.js';
import { runnerOperations } from './runner.js';
import { type HostSettings, Raise, Writer } from './writer.js';

export interface HostOptions {
  readonly clock?: Clock;
  readonly leaseMs?: number;
  readonly launchMs?: number;
  readonly schemaVersion?: number;
}

const systemClock: Clock = {
  now: () => Date.now(),
  sleep: (ms, signal) =>
    new Promise((resolve, reject) => {
      const timer = setTimeout(resolve, ms);
      signal?.addEventListener('abort', () => {
        clearTimeout(timer);
        reject(signal.reason);
      });
    }),
};

export interface Operation<R> {
  readonly key?: string;
  readonly input?: Json;
  /** The first call's result, rebuilt from the events written under its key. */
  replay?(events: readonly HarnessEvent[], tx: StorageTx): Promise<R>;
  /** A `Raise` commits the writes, then throws its error. */
  run(writer: Writer): Promise<R | Raise>;
}

export interface Host {
  readonly storage: HarnessStorage;
  readonly clock: Clock;
  readonly settings: HostSettings;
  mutate<R>(harnessId: HarnessId, operation: Operation<R>): Promise<R>;
  read<R>(harnessId: HarnessId, run: (tx: StorageTx) => Promise<R>): Promise<R>;
  /** Stores the messages whose hashes are new; hashing happens before the transaction. */
  hashed(messages: readonly Message[]): Promise<{ id: string; message: Message }[]>;
}

/** The host protocol over a storage backend (harness.md §12). */
export function createHarnessHost(storage: HarnessStorage, options: HostOptions = {}): HarnessHost {
  const clock = options.clock ?? systemClock;
  const settings: HostSettings = {
    leaseMs: options.leaseMs ?? 30_000,
    launchMs: options.launchMs ?? 60_000,
    schemaVersion: options.schemaVersion ?? SCHEMA_VERSION,
  };

  const open = async (tx: StorageTx, harnessId: HarnessId) => {
    const [execution] = await tx.execution.get([harnessId]);
    if (execution === undefined || (await tx.store.tombstone(harnessId)))
      throw new HarnessError('invalid', `no harness ${harnessId}`);
    return execution;
  };

  const host: Host = {
    storage,
    clock,
    settings,
    async mutate<R>(harnessId: HarnessId, operation: Operation<R>): Promise<R> {
      const inputDigest = operation.key !== undefined ? await digest(operation.input ?? null) : undefined;
      const outcome = await storage.transact(
        harnessId,
        async (tx): Promise<{ result: R } | { error: HarnessError }> => {
          const execution = await open(tx, harnessId);
          if (operation.key !== undefined) {
            const prior = await tx.events.byKey(operation.key);
            if (prior.length > 0) {
              if (prior[0]?.inputDigest !== inputDigest)
                throw new HarnessError('invalid', `key ${operation.key} was used with other input`);
              if (operation.replay === undefined) throw new HarnessError('invalid', `key ${operation.key} was used`);
              return { result: await operation.replay(prior, tx) };
            }
          }
          const writer = new Writer(tx, execution, clock.now(), settings, operation.key, inputDigest);
          const value = await operation.run(writer);
          if (writer.changed) tx.execution.put(harnessId, writer.exec);
          return value instanceof Raise ? { error: value.error } : { result: value as R };
        },
      );
      if ('error' in outcome) throw outcome.error;
      return outcome.result;
    },
    // Reads check nothing: a removed harness's records are dropped, so they read as empty.
    read: (harnessId, run) => storage.transact(harnessId, run),
    hashed: (messages) => Promise.all(messages.map(async (message) => ({ id: await hashMessage(message), message }))),
  };

  const reader = (harnessId: HarnessId): HarnessReader => ({
    id: harnessId,
    execution: () => host.read(harnessId, (tx) => open(tx, harnessId)),
    heads: (ids) =>
      host.read(harnessId, async (tx) => {
        const records = await tx.heads.get(ids);
        const out: Head[] = [];
        for (const [i, record] of records.entries()) {
          if (record === undefined) throw new HarnessError('invalid', `no head ${ids[i]}`);
          const segments = await tx.segments.get(record.segments);
          out.push({ id: record.id, entries: segments.flatMap((segment) => segment?.entries ?? []) });
        }
        return out;
      }),
    pointers: (scopes) =>
      host.read(harnessId, async (tx) => {
        const [execution] = await tx.execution.get([harnessId]);
        const stored = await tx.pointers.get(scopes.filter((scope) => scope !== MAIN));
        return scopes.map((scope) =>
          scope === MAIN ? { scope, head: execution?.main ?? '' } : stored.find((p) => p?.scope === scope),
        );
      }),
    entries: (ids) =>
      host.read(harnessId, async (tx) => {
        const messages = await tx.messages.get(ids);
        return messages.map((message, i) => {
          if (message === undefined) throw new HarnessError('integrity', `no message ${ids[i]}`);
          return message;
        });
      }),
    events: ({ after = 0, limit = 100 }) => host.read(harnessId, (tx) => tx.events.after(after, limit)),
    steps: (ids) => host.read(harnessId, (tx) => tx.events.byStep(ids)),
    turns: ({ before, limit = 20 }) =>
      host.read(harnessId, async (tx) => {
        const out: TurnSummary[] = [];
        const execution = await open(tx, harnessId);
        if (before === undefined && execution?.turn !== undefined) {
          const turn = execution.turn;
          const [begin] = await tx.events.byTurn(turn.id, undefined, 1);
          out.push({
            id: turn.id,
            status: 'active',
            startedBy: begin?.kind === 'step' ? (begin.begun?.items ?? []) : [],
            firstStep: begin?.kind === 'step' ? begin.step : turn.pending.id,
            usage: turn.usage,
          });
        }
        for (const event of await tx.events.endedTurns(before, limit - out.length))
          if (event.summary !== undefined) out.push(event.summary);
        return out;
      }),
    blobs: async (refs) => {
      const found = await storage.getBlobs(harnessId, refs);
      return found.map((body, i) => {
        if (body === undefined) throw new HarnessError('integrity', `no blob ${refs[i]}`);
        return body;
      });
    },
    watch: () => storage.watch(harnessId),
  });

  return {
    async create(harnessId, input) {
      const empty = await emptyHead();
      await storage.transact(harnessId, async (tx) => {
        if (await tx.store.tombstone(harnessId)) throw new HarnessError('invalid', `harness ${harnessId} was removed`);
        const [existing] = await tx.execution.get([harnessId]);
        if (existing !== undefined) return;
        tx.heads.put(empty.id, empty);
        tx.execution.put(harnessId, {
          harnessId,
          seq: 0,
          schemaVersion: settings.schemaVersion,
          main: empty.id,
          inbox: [],
          runner: { phase: 'none' },
        });
        void input;
      });
    },
    async remove(harnessId) {
      await storage.transact(harnessId, async (tx) => {
        tx.store.setTombstone(harnessId);
      });
      await storage.drop(harnessId);
    },
    client: (harnessId): HarnessClient => ({ ...reader(harnessId), ...clientOperations(host, harnessId) }),
    runnerPort: (harnessId): RunnerPort => ({ ...reader(harnessId), ...runnerOperations(host, harnessId) }),
  };
}
