/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { expect, it } from 'vitest';
import type { CallRecord, HarnessEvent, HarnessReader } from '../../harness/index.js';
import type { Json, Message } from '../../messages/index.js';
import { assistant, declareTools, system, toolDefinition, user } from '../../testing/index.js';
import {
  checkLinear,
  quadraticControl,
  scaling,
  scalingAsync,
} from '../../../../sequences/core/_test/model/scaling.js';
import { createTraceCache, readTrace } from '../index.js';

// The Inspector's reads (`test:scaling`, `sequences/core/_test/model/scaling.ts`; trace.md, Incremental reads): one turn of N model
// calls, each with a tool call, each seeing the conversation so far (a window that grows by two messages a call), read
// again as the turn goes on. A read that maps every call's window again is calls × messages a read; one that keeps the
// views of the calls it has read maps the new call's window only (one window: about the factor, 4×, not its square).

const TURN = 'turn-1';
const attempt = (at: number) => ({
  id: `a-${at}`,
  startedAt: at,
  endedAt: at + 5,
  outcome: 'succeeded' as const,
  usage: { modelCalls: 1, toolCalls: 0, uncachedInput: 100, cacheRead: 0, cacheWrite: 0, output: 10, reasoning: 0 },
});

/** A harness reader over a synthetic log of one turn: calls are appended with `grow`, as a running turn writes them. */
function syntheticReader() {
  const messages: Message[] = [
    system('Be brief.'),
    declareTools([toolDefinition('echo')]),
    {
      ...user('Request: go'),
      label: { kind: 'session:request', data: { prompt: 'go' } },
    } as Message,
  ];
  const ids = messages.map((_, i) => `m${i}`);
  const heads = new Map<string, readonly string[]>();
  const blobs = new Map<string, Json>();
  const events: HarnessEvent[] = [];
  let calls = 0;
  const entries = new Map(messages.map((m, i) => [`m${i}`, m]));
  const head = (n: number) => {
    const id = `h${n}`;
    if (!heads.has(id)) heads.set(id, ids.slice(0, n));
    return id;
  };
  const step = (kind: 'begin' | 'model' | 'tools', callRecords: CallRecord[], to?: string): HarnessEvent =>
    ({
      kind: 'step',
      seq: events.length + 1,
      at: 1000 + events.length * 10,
      turnId: TURN,
      step: `s${events.length}`,
      stepKind: kind,
      ending: 'settled',
      calls: callRecords,
      moves: to ? [{ pointer: 'main', to }] : [],
    }) as HarnessEvent;
  events.push(step('begin', [], head(3)));
  const push = (message: Message) => {
    const id = `m${ids.length}`;
    ids.push(id);
    entries.set(id, message);
  };
  const reader: HarnessReader & { grow(count: number): void; calls(): number } = {
    id: 'synthetic',
    execution: async () => ({ harnessId: 'synthetic', seq: events.length, main: head(ids.length) }) as never,
    heads: async (hs) => hs.map((h) => ({ id: h, entries: heads.get(h) ?? [] })) as never,
    pointers: async () => [],
    entries: async (es) => es.map((e) => entries.get(e) as Message),
    // Seqs are positions: a reader that finds the events after one without scanning the log.
    events: async ({ after = 0, limit = 500 }) => events.slice(after, after + limit),
    steps: async () => [],
    turns: async () => [],
    blobs: async (refs) => refs.map((r) => blobs.get(r) as Json),
    watch: (() => {
      throw new Error('not watched');
    }) as never,
    calls: () => calls,
    grow(count) {
      for (let n = 0; n < count; n += 1) {
        const i = calls++;
        const h = head(ids.length);
        const at = 1000 + events.length * 10;
        const request = `req-${i}`;
        const prefix = [`tools`, `system`, ...ids.slice(2).map((m) => `${m}`)];
        blobs.set(request, { head: h, configuration: { target: { model: 'scripted-1' } }, prefix });
        blobs.set(`res-${i}`, assistant(`call ${i}`) as unknown as Json);
        blobs.set(`treq-${i}`, { toolId: 'echo', arguments: { n: i } });
        blobs.set(`tres-${i}`, { status: 'succeeded', parts: [{ kind: 'text', text: 'echoed' }] });
        const model: CallRecord = {
          id: `c-${i}`,
          key: `k-${i}`,
          kind: 'model',
          request,
          result: `res-${i}`,
          attempts: [attempt(at)],
        } as never;
        push(assistant(`call ${i}`));
        events.push(step('model', [model], head(ids.length)));
        const tool: CallRecord = {
          id: `t-${i}`,
          key: `tk-${i}`,
          kind: 'tool',
          request: `treq-${i}`,
          result: `tres-${i}`,
          attempts: [attempt(at + 6)],
        } as never;
        push(user(`result ${i}`));
        events.push(step('tools', [tool], head(ids.length)));
      }
    },
  };
  return reader;
}

const read = (reader: HarnessReader, cache = createTraceCache()) => readTrace(reader, { cache, window: 200_000 });

it('a read after one more call of a turn of N calls', async () => {
  const s = await scalingAsync(
    {
      setup: async (n) => {
        const reader = syntheticReader();
        reader.grow(n);
        const cache = createTraceCache();
        const first = await read(reader, cache);
        expect(first.turns[0]?.steps).toHaveLength(2 * n);
        reader.grow(1);
        return { reader, cache, n };
      },
      run: async ({ reader, cache, n }) => {
        const trace = await read(reader, cache);
        expect(trace.turns[0]?.steps).toHaveLength(2 * n + 2);
      },
    },
    250,
    { reps: 7, alloc: true },
  );
  checkLinear('trace: a read after one more call of N', s);
});

it('a read when nothing happened, in a turn of N calls', async () => {
  const s = await scalingAsync(
    {
      setup: async (n) => {
        const reader = syntheticReader();
        reader.grow(n);
        const cache = createTraceCache();
        await read(reader, cache);
        return { reader, cache };
      },
      run: async ({ reader, cache }) => {
        for (let i = 0; i < 20; i += 1) await read(reader, cache);
      },
    },
    250,
    { reps: 7, alloc: true },
  );
  checkLinear('trace: 20 reads, nothing new, N calls', s);
});

it('a read from nothing of a turn of N calls is the windows read once each: caught as growing with the square', async () => {
  // The positive control: the same read as the first, with no cache kept (every window mapped again).
  const s = await scalingAsync(
    {
      setup: (n) => {
        const reader = syntheticReader();
        reader.grow(n);
        return { reader };
      },
      run: async ({ reader }) => {
        await read(reader);
      },
    },
    100,
    { reps: 3, alloc: true },
  );
  checkLinear('trace: control, a read from nothing of N calls', s, false);
});

it('the quadratic control is caught', async () => {
  checkLinear('control', await scaling(quadraticControl, 2000), false);
});
