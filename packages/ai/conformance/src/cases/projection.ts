/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { digest, type Json, type Message, type ProjectInput } from '@molstar/ai-agent';
import { type ConformanceFixture, createSession, type Session } from '../session.js';
import { settle } from '../testing/clock.js';
import { declareTools, system, textOf, toolDefinition, user } from '../testing/messages.js';
import { lastStep, stepEvents, with_ } from './helpers.js';

const WINDOWS = { leaseMs: 30_000, launchMs: 10_000 };
/** A small window, so a long prompt does not fit (harness.md §11, Admission). */
const SMALL_WINDOW = { contextWindowTokens: 10_000, maxOutputTokens: 1_000, defaultMaxOutputTokens: 1_000 };
const LONG = user('x'.repeat(60_000));

/** A crash and the relaunch after the lease runs out (harness.md §4); returns the new launch key. */
async function crash(s: Session, launch: string): Promise<string> {
  const runner = s.runner();
  const running = runner.run(launch);
  await settle(40);
  runner.abort('crash');
  await s.clock.runUntil(running);
  await s.clock.advance(WINDOWS.leaseMs + 1_000);
  const checked = await s.client.check();
  expect(checked.launch).toBeDefined();
  return checked.launch as string;
}

async function sha256Hex(text: string): Promise<string> {
  const bytes = new Uint8Array(await crypto.subtle.digest('SHA-256', new TextEncoder().encode(text)));
  return Array.from(bytes, (byte) => byte.toString(16).padStart(2, '0')).join('');
}

const toolIds = (tools: readonly { readonly tool: { readonly id: string } }[]) => tools.map((tool) => tool.tool.id);

/** The serialized messages up to, not including, the closing bracket: what a longer list must start with. */
const bytesOf = (messages: readonly Message[]) => JSON.stringify(messages).slice(0, -1);

/** C1–C2 (Projection). */
export function projectionCases(fixture: ConformanceFixture): void {
  describe('Projection', () => {
    // C1 (harness.md §11)
    it('C1 same-input-same-digest: one head projected twice (a step re-run after a crash) gives one digest', async () => {
      const s = await createSession(fixture, {
        host: WINDOWS,
        model: [{ text: 'lost', waitMs: 60_000 }, { text: 'again' }],
      });
      const { launch } = await s.send('go', { run: false });
      await s.run(await crash(s, launch as string));
      expect(s.model.invocations).toHaveLength(2);
      const [first, second] = s.model.invocations;
      expect(second?.digest).toBe(first?.digest);
      expect(second?.payload).toEqual(first?.payload);
      expect(s.model.projections[1]?.messages).toEqual(s.model.projections[0]?.messages);
    });

    it("C1 append-prefix: a later request's messages start with an earlier one's", async () => {
      const s = await createSession(fixture, { model: [{ text: 'one' }, { text: 'two' }] });
      await s.send('first');
      await s.send('second');
      const [first, second] = s.model.projections;
      expect(first).toBeDefined();
      expect(second?.messages.slice(0, first?.messages.length)).toEqual(first?.messages);
    });

    it('C1 append-prefix: after a model-and-tool round in one turn, the earlier payload is a byte prefix of the new one', async () => {
      const s = await createSession(fixture, {
        tools: { look: () => ({ text: 'seen' }) },
        model: [{ calls: [{ tool: 'look', arguments: {} }] }, { text: 'done' }],
      });
      await s.send('go');
      const [first, second] = s.model.projections;
      expect(second?.messages.length).toBeGreaterThan(first?.messages.length ?? 0);
      expect(JSON.stringify(second?.messages).startsWith(bytesOf(first?.messages ?? []))).toBe(true);
      expect(toolIds(second?.tools ?? []).slice(0, toolIds(first?.tools ?? []).length)).toEqual(
        toolIds(first?.tools ?? []),
      );
    });

    it('C1 deferred-load-appended: a deferred tool loaded later is appended to the head and the tools, not hoisted', async () => {
      const later = toolDefinition('later');
      const s = await createSession(fixture, {
        tools: { loader: () => ({ text: 'loaded', contributions: [declareTools([later])] }), later: () => ({}) },
        declare: ['loader'],
        deferred: ['later'],
        model: [{ calls: [{ tool: 'loader', arguments: {} }] }, { text: 'done' }],
      });
      await s.send('go');
      const [before, after] = s.model.projections;
      expect(toolIds(before?.tools ?? [])).toEqual(['loader']);
      expect(toolIds(after?.tools ?? [])).toEqual(['loader', 'later']);
      expect(JSON.stringify(after?.messages).startsWith(bytesOf(before?.messages ?? []))).toBe(true);
      const loadedAt = (after?.messages ?? []).findIndex((m) =>
        m.parts.some((p) => p.kind === 'tool-definition' && p.definition.tool.id === 'later' && !p.deferred),
      );
      expect(loadedAt).toBeGreaterThanOrEqual(before?.messages.length ?? Number.POSITIVE_INFINITY);
    });

    it("C1 digest-jcs: digest matches RFC 8785's vectors, integer-like keys included (defect 32)", async () => {
      // RFC 8785 section 3.2.3: keys sort by UTF-16 code units; the expected text lists them in that order.
      const sorting = {
        '€': 'Euro Sign',
        '\r': 'Carriage Return',
        דּ: 'Hebrew Letter Dalet With Dagesh',
        '1': 'One',
        '😀': 'Emoji: Grinning Face',
        '\u0080': 'Control',
        ö: 'Latin Small Letter O With Diaeresis',
      };
      const order = ['\r', '1', '\u0080', 'ö', '€', '😀', 'דּ'] as const;
      const sorted = `{${order.map((key) => `${JSON.stringify(key)}:${JSON.stringify(sorting[key])}`).join(',')}}`;
      expect(await digest(sorting)).toBe(await sha256Hex(sorted));
      // Integer-like keys sort as strings ("10" before "2"), not as numbers, and not in insertion order.
      expect(await digest({ b: 1, 10: 2, 9: 3, 2: 4 })).toBe(await sha256Hex('{"10":2,"2":4,"9":3,"b":1}'));
      // RFC 8785 section 3.2.2: ECMAScript number and string serialization.
      const values = JSON.parse(
        '{"numbers":[333333333.33333329,1E30,4.50,2e-3,0.000000000000000000000000001],"string":"\\u20ac$\\u000F\\u000aA\'B\\u0022\\u005c\\\\\\"\\/","literals":[null,true,false]}',
      ) as Json;
      const expected = String.raw`{"literals":[null,true,false],"numbers":[333333333.3333333,1e+30,4.5,0.002,1e-27],"string":"€$\u000f\nA'B\"\\\\\"/"}`;
      expect(await digest(values)).toBe(await sha256Hex(expected));
    });

    it('C1 scope-ignored-when-unsupported: without the capability, the payload is the same with or without a scope', async () => {
      const payload = (input: ProjectInput): Json => ({
        messages: input.messages as unknown as Json,
        // An adapter maps the scope only when its capabilities say so (harness.md §11).
        ...(input.capabilities.cacheScope && input.configuration.cacheScopeId !== undefined
          ? { scope: input.configuration.cacheScopeId }
          : {}),
      });
      const target = { provider: 'scripted', model: 'scripted-1' };
      const digestOf = async (cacheScope: boolean, cacheScopeId?: string) => {
        const s = await createSession(fixture, {
          model: [{ text: 'ok' }],
          modelOptions: { capabilities: { cacheScope }, payload },
          configuration: {
            model: { target, parameters: {}, ...(cacheScopeId !== undefined ? { cacheScopeId } : {}) },
          },
        });
        await s.send('go');
        return s.model.invocations[0]?.digest;
      };
      const without = await digestOf(false);
      expect(without).toBeDefined();
      expect(await digestOf(false, 'scope-a')).toBe(without);
      // The scope does reach the adapter: with the capability it changes the request.
      expect(await digestOf(true, 'scope-a')).not.toBe(await digestOf(true));
    });

    // C2 (harness.md §11, Admission; §7.1)
    it('C2 over-window: failure with capacity opens before any attempt; compact applies the edit and the step projects again', async () => {
      const codes: string[] = [];
      const s = await createSession(fixture, {
        modelOptions: { capabilities: SMALL_WINDOW },
        model: [{ text: 'ok' }],
        drive: with_({
          failure: ({ payload }) => {
            codes.push(payload.failure.code);
            return { compact: { entries: [user('short')] } };
          },
        }),
      });
      await s.send([LONG]);
      expect(codes).toEqual(['capacity']);
      expect(s.model.projections).toHaveLength(2);
      expect(s.model.invocations).toHaveLength(1);
      expect(s.model.projections[1]?.messages).toEqual([user('short')]);
      expect(stepEvents(await s.events()).some((e) => e.compaction?.some((c) => c.pointer === 'main'))).toBe(true);
      const main = await s.main();
      expect(textOf(main[0] as Message)).toBe('short');
      expect((await s.turns())[0]?.status).toBe('completed');
    });

    it('C2 invalid-edit: a compact whose edit leaves an unanswered call is an invalid answer', async () => {
      const dangling: Message = {
        role: 'assistant',
        parts: [{ kind: 'tool-call', requestId: 'r1', tool: { id: 'ghost' }, arguments: {} }],
      };
      const s = await createSession(fixture, {
        modelOptions: { capabilities: SMALL_WINDOW },
        model: [{ text: 'ok' }],
        drive: with_({ failure: () => ({ compact: { entries: [dangling] } }) }),
      });
      await s.send([LONG]);
      expect((await s.turns())[0]?.status).toBe('failed');
      expect(JSON.stringify(lastStep(await s.events()).end?.reason)).toContain('invalid-answer');
      expect(s.model.invocations).toHaveLength(0);
      expect(await s.main()).toEqual([LONG]);
    });

    it('C2 unrepresentable: a role the adapter cannot place yields unrepresentable, never a plain error (defect 15c)', async () => {
      const codes: string[] = [];
      const s = await createSession(fixture, {
        modelOptions: {
          unrepresentable: (input) => (input.messages.some((m) => m.role === 'system') ? 'no system role' : undefined),
        },
        model: [{ text: 'ok' }],
        drive: with_({
          failure: ({ payload }) => {
            codes.push(payload.failure.code);
            return 'pass';
          },
        }),
      });
      await s.send([system('be terse'), user('go')]);
      expect(codes).toEqual(['unrepresentable']);
      expect(s.model.invocations).toHaveLength(0);
      expect((await s.turns())[0]?.status).toBe('failed');
    });

    it('C2 unknown-tool-call: a call to a tool not in the request keeps its name and opens repair with unknown-tool (defect 15a)', async () => {
      const seen: { reason: string; toolId: string }[] = [];
      const s = await createSession(fixture, {
        model: [{ calls: [{ tool: 'ghost', arguments: {} }] }, { text: 'ok' }],
        drive: with_({
          repair: ({ payload }) => {
            seen.push({ reason: payload.reason, toolId: payload.call.toolId });
            return { report: { parts: [{ kind: 'text', text: 'no such tool' }] } };
          },
        }),
      });
      await s.send('go');
      expect(seen).toEqual([{ reason: 'unknown-tool', toolId: 'ghost' }]);
      const main = await s.main();
      const parts = main.flatMap((m): Message['parts'][number][] => [...m.parts]);
      const call = parts.find((p) => p.kind === 'tool-call');
      expect(call?.kind === 'tool-call' ? call.tool.id : undefined).toBe('ghost');
      const result = parts.find((p) => p.kind === 'tool-result');
      expect(result?.kind === 'tool-result' ? result.status : undefined).toBe('failed');
      expect((await s.turns())[0]?.status).toBe('completed');
    });
  });
}
