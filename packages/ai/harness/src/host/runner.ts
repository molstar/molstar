/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { canonicalJson, type Json, type Message } from '@molstar/ai-agent';
import {
  type CommitBatch,
  HarnessError,
  type HarnessReader,
  type Lease,
  type Move,
  type RunnerPort,
} from '../contracts/ports.js';
import { type HarnessEvent, type HarnessId, type ItemKey, MAIN, type Progress } from '../contracts/records.js';
import { launchResult } from './client.js';
import type { Host } from './index.js';
import { addUsage, Raise, usageOf, type Writer, ZERO } from './writer.js';

type RunnerOperations = Omit<RunnerPort, keyof HarnessReader>;

/** The live lease, or `stale`; then extended, with the runner's progress (harness.md §4). */
function hold(writer: Writer, lease: Lease, progress?: Progress): Lease {
  const runner = writer.exec.runner;
  if (runner.phase !== 'live' || runner.runnerId !== lease.runnerId || runner.leaseUntil <= writer.now)
    throw new HarnessError('stale', 'the lease is no longer current');
  const until = writer.now + writer.settings.leaseMs;
  writer.set('runner', {
    phase: 'live',
    runnerId: runner.runnerId,
    leaseUntil: until,
    signalAt: writer.now,
    ...(progress !== undefined ? { progress } : runner.progress !== undefined ? { progress: runner.progress } : {}),
  });
  return { harnessId: writer.harnessId, runnerId: runner.runnerId, until };
}

/** The current lease for a replayed operation, if the runner still holds it. */
function currentLease(
  writer: { exec: { harnessId: HarnessId; runner: Writer['exec']['runner'] } },
  lease: Lease,
): Lease {
  const runner = writer.exec.runner;
  if (runner.phase === 'live' && runner.runnerId === lease.runnerId)
    return { harnessId: writer.exec.harnessId, runnerId: runner.runnerId, until: runner.leaseUntil };
  return lease;
}

/** Stores the new messages, head records, segments and pointer records of a commit (harness.md §12). */
async function applyMoves(
  writer: Writer,
  messages: readonly { id: string; message: Message }[],
  moves: readonly Move[],
) {
  const seen = new Set<string>();
  for (const move of moves) {
    if (seen.has(move.pointer)) throw new HarnessError('invalid', `two moves of ${move.pointer} in one commit`);
    seen.add(move.pointer);
    const current = move.pointer === MAIN ? writer.exec.main : (await writer.tx.pointers.get([move.pointer]))[0]?.head;
    if (move.from !== current) throw new HarnessError('head-changed', `${move.pointer} moved since ${move.from}`);
  }
  const ids = messages.map((m) => m.id);
  const existing = await writer.tx.messages.get(ids);
  for (const [i, entry] of messages.entries())
    if (existing[i] === undefined) writer.tx.messages.put(entry.id, entry.message);
  for (const move of moves) {
    for (const segment of move.segments) writer.tx.segments.put(segment.id, segment);
    writer.tx.heads.put(move.head.id, move.head);
    if (move.pointer === MAIN) writer.set('main', move.head.id);
    else writer.tx.pointers.put(move.pointer, { scope: move.pointer, head: move.head.id });
  }
}

function consume(writer: Writer, keys: readonly ItemKey[]): void {
  for (const key of keys)
    if (!writer.exec.inbox.some((item) => item.key === key))
      throw new HarnessError('inbox-changed', `item ${key} no longer waits`);
  writer.set(
    'inbox',
    writer.exec.inbox.filter((item) => !keys.includes(item.key)),
  );
}

/** A commit's contents must match its ending (harness.md §6). */
function checkShape(batch: CommitBatch): void {
  const bad = (why: string) => {
    throw new HarnessError('invalid', `a ${batch.ending} commit ${why}`);
  };
  if (batch.ending === 'settled' && (batch.next === undefined) === (batch.end === undefined))
    bad('names exactly one of next and end');
  if (batch.ending === 'parked') {
    if (batch.moves.length > 0 || batch.consumed.length > 0) bad('moves nothing and consumes nothing');
    if (batch.park === undefined || batch.next !== undefined || batch.end !== undefined) bad('carries a park only');
  }
  if (batch.ending === 'interrupted') {
    if (batch.moves.some((move) => move.pointer !== MAIN)) bad('moves only main');
    if (batch.next?.kind !== 'model' || batch.end !== undefined) bad('names a model step');
  }
  if (batch.ending === 'abandoned') {
    if (batch.moves.length > 0 || batch.consumed.length > 0) bad('moves nothing and consumes nothing');
    if (batch.end === undefined || batch.next !== undefined) bad('ends the turn');
  }
}

const committed = (events: readonly HarnessEvent[]) => events.at(-1)?.seq ?? 0;

export function runnerOperations(host: Host, harnessId: HarnessId): RunnerOperations {
  const records = <T>(
    read: (ids: readonly string[]) => Promise<(T | undefined)[]>,
    ids: readonly string[],
    what: string,
  ) =>
    read(ids).then((found) =>
      found.map((record, i) => {
        if (record === undefined) throw new HarnessError('integrity', `no ${what} ${ids[i]}`);
        return record;
      }),
    );
  return {
    headRecords: (ids) => host.read(harnessId, (tx) => records((x) => tx.heads.get(x), ids, 'head')),
    segments: (ids) => host.read(harnessId, (tx) => records((x) => tx.segments.get(x), ids, 'segment')),
    lastCompaction: (pointer) =>
      host.read(harnessId, async (tx) => {
        const event = await tx.events.lastCompaction(pointer);
        const entry = event?.compaction?.findLast((c) => c.pointer === pointer);
        return event !== undefined && entry !== undefined ? { seq: event.seq, sizeAfter: entry.sizeAfter } : undefined;
      }),
    claim: (input) =>
      host.mutate(harnessId, {
        key: `claim:${input.launch}:${input.runnerId}`,
        input: input as unknown as Json,
        replay: async (events, tx) => {
          const [execution] = await tx.execution.get([harnessId]);
          const claimed = events.find((event) => event.kind === 'control' && event.control === 'claimed');
          const until = claimed?.kind === 'control' ? (claimed.detail as { until?: number }).until : undefined;
          if (execution === undefined || until === undefined)
            throw new HarnessError('stale', 'the launch was not claimed');
          return { lease: { harnessId, runnerId: input.runnerId, until }, execution };
        },
        async run(writer) {
          const runner = writer.exec.runner;
          if (runner.phase !== 'launching' || runner.key !== input.launch)
            throw new HarnessError('stale', `launch ${input.launch} is not current`);
          const turn = writer.exec.turn;
          if (turn?.stopRequested !== undefined) await writer.endByControl('stopped', turn.stopRequested.reason);
          const active = writer.exec.turn;
          if (active !== undefined && active.blocked === undefined && active.schemaVersion !== input.schemaVersion) {
            writer.set('turn', { ...active, blocked: 'stale-schema' });
            writer.set('runner', { phase: 'none' });
            writer.control('blocked', { turnId: active.id, detail: { schemaVersion: input.schemaVersion } });
            return new Raise(new HarnessError('stale-schema', `turn ${active.id} has schema ${active.schemaVersion}`));
          }
          const until = writer.now + writer.settings.leaseMs;
          writer.set('runner', { phase: 'live', runnerId: input.runnerId, leaseUntil: until, signalAt: writer.now });
          writer.control('claimed', { detail: { runnerId: input.runnerId, launch: input.launch, until } });
          return { lease: { harnessId, runnerId: input.runnerId, until }, execution: writer.exec };
        },
      }),

    async putBlobs(lease, blobs) {
      // The lease is checked, not extended: blobs are content-addressed, so a late write is harmless.
      await host.read(harnessId, async (tx) => {
        const [execution] = await tx.execution.get([harnessId]);
        const runner = execution?.runner;
        if (runner?.phase !== 'live' || runner.runnerId !== lease.runnerId || runner.leaseUntil <= host.clock.now())
          throw new HarnessError('stale', 'the lease is no longer current');
      });
      await host.storage.putBlobs(harnessId, blobs);
    },

    async begin(lease, input) {
      const messages = await host.hashed(input.messages);
      return host.mutate(harnessId, {
        key: input.key,
        input: input as unknown as Json,
        replay: async (events, tx) => {
          const [execution] = await tx.execution.get([harnessId]);
          return { seq: committed(events), lease: currentLease({ exec: execution as never }, lease) };
        },
        async run(writer) {
          const extended = hold(writer, lease, input.progress);
          if (writer.exec.turn !== undefined) throw new HarnessError('invalid', 'a turn is already active');
          if (input.move.pointer !== MAIN) throw new HarnessError('invalid', 'begin moves main');
          consume(writer, input.items);
          await applyMoves(writer, messages, [input.move]);
          const { all } = usageOf(input.calls);
          writer.set('turn', {
            id: input.turnId,
            schemaVersion: writer.settings.schemaVersion,
            configuration: input.configuration,
            budget: { limits: input.configuration.limits, consumed: ZERO },
            pending: input.next,
            usage: all,
          });
          writer.step({
            turnId: input.turnId,
            step: input.step,
            stepKind: 'begin',
            ending: 'settled',
            calls: input.calls,
            moves: [
              {
                pointer: MAIN,
                ...(input.move.from !== undefined ? { from: input.move.from } : {}),
                to: input.move.head.id,
              },
            ],
            ...(input.compaction !== undefined ? { compaction: input.compaction } : {}),
            consumed: input.items,
            next: input.next,
            begun: { items: input.items, configuration: input.configuration },
          });
          return { seq: writer.exec.seq, lease: extended };
        },
      });
    },

    async commit(lease, batch) {
      const messages = await host.hashed(batch.messages);
      return host.mutate(harnessId, {
        key: batch.key,
        input: batch as unknown as Json,
        replay: async (events, tx) => {
          const [execution] = await tx.execution.get([harnessId]);
          return { seq: committed(events), lease: currentLease({ exec: execution as never }, lease) };
        },
        async run(writer) {
          const extended = hold(writer, lease, batch.progress);
          const turn = writer.turn(batch.turnId);
          if (turn.pending.id !== batch.step) throw new HarnessError('stale', `step ${batch.step} is not pending`);
          if (turn.stopRequested !== undefined && (batch.ending !== 'abandoned' || batch.end?.status !== 'stopped'))
            throw new HarnessError('stale', 'a Stop was requested: abandon the step');
          checkShape(batch);
          const steered = batch.consumed.flatMap(
            (key) => writer.exec.inbox.find((item) => item.key === key)?.messages ?? [],
          );
          consume(writer, batch.consumed);
          if (
            batch.ending === 'interrupted' &&
            canonicalJson(steered as unknown as Json) !== canonicalJson(batch.messages as unknown as Json)
          )
            throw new HarnessError('invalid', 'an interrupted commit adds only the steered items');
          await applyMoves(writer, messages, batch.moves);
          const usage = usageOf(batch.calls);
          const limits = batch.limits ?? turn.budget.limits;
          const updated = {
            ...turn,
            usage: addUsage(turn.usage, usage.all),
            budget: { limits, consumed: addUsage(turn.budget.consumed, usage.charged) },
          };
          writer.set('turn', updated);
          const summary =
            batch.end !== undefined ? await writer.summary(batch.end.status, updated.usage, batch.step) : undefined;
          writer.step({
            turnId: turn.id,
            step: batch.step,
            stepKind: turn.pending.kind,
            ending: batch.ending,
            calls: batch.calls,
            moves: batch.moves.map((move) => ({
              pointer: move.pointer,
              ...(move.from !== undefined ? { from: move.from } : {}),
              to: move.head.id,
            })),
            ...(batch.compaction !== undefined ? { compaction: batch.compaction } : {}),
            consumed: batch.consumed,
            ...(batch.next !== undefined ? { next: batch.next } : {}),
            ...(batch.end !== undefined ? { end: batch.end } : {}),
            ...(summary !== undefined ? { summary } : {}),
            ...(batch.limits !== undefined ? { limits: batch.limits } : {}),
          });
          if (batch.end !== undefined) writer.set('turn', undefined);
          else if (batch.ending === 'parked' && batch.park !== undefined)
            writer.set('turn', { ...updated, park: batch.park });
          else if (batch.next !== undefined) writer.set('turn', { ...updated, pending: batch.next });
          return { seq: writer.exec.seq, lease: extended };
        },
      });
    },

    renew: (lease, input) => host.mutate(harnessId, { run: async (writer) => hold(writer, lease, input?.progress) }),

    release: (lease, input) =>
      host.mutate(harnessId, {
        key: `release:${lease.runnerId}:${input.seenSeq}`,
        input: input as unknown as Json,
        replay: async (events) => launchResult(events),
        async run(writer) {
          hold(writer, lease);
          const since = await writer.tx.events.after(input.seenSeq, 10_000);
          if (since.some((event) => event.kind === 'control' && event.control === 'submitted'))
            throw new HarnessError('inbox-changed', 'an item arrived after the runner last looked');
          if (writer.exec.turn?.stopRequested !== undefined)
            await writer.endByControl('stopped', writer.exec.turn.stopRequested.reason);
          writer.set('runner', { phase: 'none' });
          writer.control('released', { detail: { runnerId: lease.runnerId } });
          const launch = writer.launchIfNeeded();
          return { seq: writer.exec.seq, ...(launch !== undefined ? { launch } : {}) };
        },
      }),
  };
}
