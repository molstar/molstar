/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { HarnessReader } from '@molstar/ai-harness';
import { declareTools, toolDefinition, user } from '@molstar/ai-conformance/testing';
import { formatUsage, type ModelStepView, readTrace, type ToolStepView, type TurnView } from '../index.js';
import {
  adapterState,
  configured,
  createSession,
  first,
  fixture,
  request,
  runnerOf,
  turn,
  until,
  with_,
} from './helpers.js';

const models = (t: TurnView) => t.steps.filter((s): s is ModelStepView => s.kind === 'model');
const tools = (t: TurnView) => t.steps.filter((s): s is ToolStepView => s.kind === 'tool');
const shape = (t: TurnView) => t.steps.map((s) => (s.kind === 'tool' ? s.name : 'model'));

/** echo succeeds, boom fails: a model that calls both, then answers. */
async function basic() {
  const s = await createSession(fixture, {
    declare: [],
    tools: {
      echo: () => ({ text: 'echoed' }),
      boom: () => ({ status: 'failed', text: 'it broke', effect: 'none' }),
    },
    model: [
      {
        calls: [
          { tool: 'echo', arguments: { a: 1 } },
          { tool: 'boom', arguments: {} },
        ],
        usage: { inputTokens: 1000, outputTokens: 20, cacheReadInputTokens: 800, reasoningOutputTokens: 5 },
      },
      {
        text: 'Done.',
        usage: { inputTokens: 1200, outputTokens: 4, cacheReadInputTokens: 1000, cacheWriteInputTokens: 100 },
        waitMs: 2000,
      },
    ],
  });
  await turn(s, first(s, request('Make it red', { state: 'Counter at 0.', focus: 'the display' })));
  return s;
}

describe('readTrace: a turn over the harness (harness.md §3, §11, §13)', () => {
  it('reads the turn, its status, request and focus, and its calls in the order they began', async () => {
    const s = await basic();
    const trace = await readTrace(s.client);
    expect(trace.harnessId).toBe(s.harnessId);
    expect(trace.turns).toHaveLength(1);
    const [t] = trace.turns as [TurnView];
    expect(t.status).toBe('done');
    expect(t.prompt).toBe('Make it red');
    expect(t.focus).toBe('the display');
    expect(shape(t)).toEqual(['model', 'echo', 'boom', 'model']);
    expect(t.failure).toBeUndefined();
    expect(trace.execution).toMatchObject({ runner: 'none' });
    expect(trace.execution.pendingStep).toBeUndefined();
  });

  it('shows a tool call with its arguments and result, and a failed one with its failure', async () => {
    const s = await basic();
    const [t] = (await readTrace(s.client)).turns as [TurnView];
    const [echo, boom] = tools(t) as [ToolStepView, ToolStepView];
    expect(echo).toMatchObject({ toolId: 'echo', args: { a: 1 }, status: 'succeeded', result: 'echoed' });
    expect(echo.attempts).toHaveLength(1);
    expect(echo.attempts[0]?.status).toBe('succeeded');
    expect(boom).toMatchObject({ toolId: 'boom', status: 'failed', result: 'it broke' });
    expect(boom.attempts[0]).toMatchObject({ ordinal: 1, status: 'failed', failure: 'it broke' });
  });

  it('shows what each model call asked for, how it ended (finish), and its usage', async () => {
    const s = await basic();
    const [t] = (await readTrace(s.client)).turns as [TurnView];
    const [a, b] = models(t) as [ModelStepView, ModelStepView];
    expect(a.model).toBe('scripted-1');
    expect(a.toolCalls).toEqual(['echo(a=1)', 'boom()']);
    expect(a.reply).toBeUndefined();
    expect(a.finish).toBe('tool-calls');
    expect(a.tokens).toEqual({ input: 1000, output: 20, cached: 800, reasoning: 5 });
    expect(b.reply).toBe('Done.');
    expect(b.toolCalls).toBeUndefined();
    expect(b.finish).toBe('complete');
    expect(b.tokens).toEqual({ input: 1200, output: 4, cached: 1000, written: 100 });
  });

  it('sums the turn’s tokens, with the cached share, and formats them in one line', async () => {
    const s = await basic();
    const [t] = (await readTrace(s.client)).turns as [TurnView];
    // Uncached 200 + 100, cache reads 1800, cache writes 100.
    expect(t.tokens).toMatchObject({ input: 2200, output: 24, cached: 1800, written: 100, reasoning: 5 });
    expect(formatUsage(t.tokens)).toBe('2.2k in (82% cached, 100 written) · 24 out');
  });

  it('measures times in ms from the first step event, and a wait shows in the attempt', async () => {
    const s = await basic();
    const [t] = (await readTrace(s.client)).turns as [TurnView];
    const last = models(t)[1] as ModelStepView;
    const attempt = last.attempts[0];
    expect((attempt?.end ?? 0) - (attempt?.start ?? 0)).toBe(2000);
    expect(t.start).toBeGreaterThanOrEqual(0);
    expect(t.end ?? 0).toBeGreaterThanOrEqual(attempt?.end ?? 0);
    expect(t.milestones.length).toBeGreaterThan(1);
    expect(t.milestones[0]?.kind).toBe('begin settled');
  });

  it('marks a count the provider did not report as estimated', async () => {
    const s = await createSession(fixture, { declare: [], model: [{ text: 'Hi.', usage: null }] });
    await turn(s, first(s, request('Hello')));
    const [t] = (await readTrace(s.client)).turns as [TurnView];
    const m = models(t)[0] as ModelStepView;
    expect(m.tokens?.estimated).toBe(true);
    expect(m.tokens?.input).toBeGreaterThan(0);
    expect(formatUsage(m.tokens as never)).toMatch(/^~/);
    expect(t.tokens.estimated).toBe(true);
  });
});

describe('readTrace: retries', () => {
  it('shows a retried model call as attempts: the failed one, the wait, the one that succeeded', async () => {
    const s = await createSession(fixture, {
      declare: [],
      model: [
        { fail: { transient: true, message: 'rate limited' } },
        { text: 'ok', usage: { inputTokens: 10, outputTokens: 2 } },
      ],
    });
    await turn(s, first(s, request('Hi')));
    const [t] = (await readTrace(s.client)).turns as [TurnView];
    expect(t.status).toBe('done');
    const [m] = models(t) as [ModelStepView];
    expect(m.attempts.map((a) => a.status)).toEqual(['failed', 'succeeded']);
    expect(m.attempts[0]?.failure).toBe('rate limited');
    expect(m.attempts[0]?.waitBefore).toBeUndefined();
    expect(m.attempts[1]?.waitBefore).toEqual({ reason: 'retry back-off', ms: 500 });
    // The call’s tokens are its last attempt’s.
    expect(m.tokens).toEqual({ input: 10, output: 2 });
    expect(m.reply).toBe('ok');
  });

  it('a turn that failed says why', async () => {
    const s = await createSession(fixture, { declare: [], model: [{ fail: { message: 'bad key' } }] });
    await turn(s, first(s, request('Hi')));
    const [t] = (await readTrace(s.client)).turns as [TurnView];
    expect(t.status).toBe('failed');
    expect(t.failure).toBe('bad key');
    expect(models(t)[0]?.attempts[0]).toMatchObject({ status: 'failed', failure: 'bad key' });
    expect(models(t)[0]?.reply).toBeUndefined();
    expect(models(t)[0]?.finish).toBeUndefined();
  });
});

describe('readTrace: nested runs', () => {
  it('records a nested run’s calls under the call that ran it, after it', async () => {
    const s = await createSession(fixture, {
      declare: [],
      tools: {
        helper: async (context) => {
          const result = await context.run({
            input: [user('work on this')],
            tools: { enabled: ['inner'], preauthorized: ['inner'], definitions: [toolDefinition('inner')] },
            model: { target: { provider: 'scripted', model: 'scripted-1' }, parameters: {} },
            limits: {},
          });
          return { text: `helper ${result.status}` };
        },
        inner: () => ({ text: 'inner done' }),
      },
      model: [
        { calls: [{ tool: 'helper', arguments: {} }] },
        { calls: [{ tool: 'inner', arguments: {} }] },
        { text: 'nested done' },
        { text: 'final' },
      ],
    });
    await turn(s, first(s, request('Use the helper')));
    const [t] = (await readTrace(s.client)).turns as [TurnView];
    expect(shape(t)).toEqual(['model', 'helper', 'model', 'inner', 'model', 'model']);
    const helper = t.steps[1] as ToolStepView;
    expect(helper.result).toBe('helper completed');
    expect(t.steps.map((x) => x.under)).toEqual([
      undefined,
      undefined,
      helper.callId,
      helper.callId,
      helper.callId,
      undefined,
    ]);
    // The turn’s tokens include the nested run’s.
    const total = models(t).reduce((n, m) => n + (m.tokens?.output ?? 0), 0);
    expect(t.tokens.output).toBe(total);
  });
});

describe('readTrace: what each model call saw (harness.md §2)', () => {
  it('lists the window in its parts, and marks what is new since the run’s previous call', async () => {
    const s = await basic();
    const [t] = (await readTrace(s.client)).turns as [TurnView];
    const [a, b] = models(t) as [ModelStepView, ModelStepView];
    const roles = (m: ModelStepView) => m.window?.entries.map((e) => `${e.role}:${e.part}`);
    expect(roles(a)).toEqual(['system:guidance', 'tools:tools', 'user:state']);
    expect(a.window?.entries.every((e) => e.fresh)).toBe(true);
    expect(roles(b)).toEqual([
      'system:guidance',
      'tools:tools',
      'user:state',
      'assistant:conversation',
      'tool:conversation',
    ]);
    expect(b.window?.entries.filter((e) => e.fresh).map((e) => e.role)).toEqual(['assistant', 'tool']);
    expect(b.window?.entries.map((e) => e.title)).toEqual([
      'Guidance',
      'echo, boom',
      'State + request',
      'Tool calls',
      'tool results',
    ]);
    const results = b.window?.entries[4];
    expect(results?.failed).toBe(true);
    expect(results?.text).toContain('echoed');
    expect(results?.tokens).toBeGreaterThan(0);
    expect(b.window?.limit).toBeUndefined();
  });

  it('marks as new, in a later turn, only what was appended since the last call', async () => {
    const s = await basic();
    s.model.push({ text: 'Again.' });
    await turn(s, 'Request: once more');
    const trace = await readTrace(s.client, { window: 200_000 });
    expect(trace.turns).toHaveLength(2);
    const [, second] = trace.turns as [TurnView, TurnView];
    const m = models(second)[0] as ModelStepView;
    expect(m.window?.limit).toBe(200_000);
    expect(m.window?.entries.filter((e) => e.fresh).map((e) => e.title)).toEqual(['Reply', 'Request']);
    expect(m.window?.entries.map((e) => e.part).slice(-2)).toEqual(['conversation', 'conversation']);
    expect(second.prompt).toBe('once more');
  });

  it('reads the preamble from the first call, and the tools loaded from the start', async () => {
    const s = await basic();
    const trace = await readTrace(s.client);
    expect(trace.preamble.guidance).toBe('Be brief.');
    expect(trace.preamble.guidanceTokens).toBeGreaterThan(0);
    expect(trace.preamble.tools.map((t) => t.id)).toEqual(['echo', 'boom']);
    expect(trace.preamble.toolsTokens).toBeGreaterThan(0);
    expect(trace.loaded).toEqual([
      { id: 'echo', name: 'echo' },
      { id: 'boom', name: 'boom' },
    ]);
    expect(trace.declared).toEqual([]);
  });

  it('gives the turn a tool was first seen loaded in, and what is declared but not loaded', async () => {
    const s = await createSession(fixture, {
      declare: [],
      tools: {
        echo: () => ({ text: 'echoed' }),
        disclose: () => ({ text: 'ok', contributions: [declareTools([toolDefinition('late')])] }),
        late: () => ({ text: 'late' }),
      },
      model: [
        { text: 'one' },
        { calls: [{ tool: 'disclose', arguments: {} }] },
        { calls: [{ tool: 'late', arguments: {} }] },
        { text: 'two' },
      ],
    });
    await turn(s, first(s, request('first'), { deferred: ['late'] }));
    const before = await readTrace(s.client);
    expect(before.declared).toEqual(['late']);
    expect(before.loaded.map((t) => t.id)).toEqual(['echo', 'disclose']);
    await turn(s, request('second'));
    const after = await readTrace(s.client);
    expect(after.declared).toEqual([]);
    expect(after.loaded.map((t) => [t.id, t.turn])).toEqual([
      ['echo', undefined],
      ['disclose', undefined],
      ['late', 2],
    ]);
    // The window shows a declared tool by name.
    const firstCall = models(after.turns[0] as TurnView)[0] as ModelStepView;
    const definitions = firstCall.window?.entries.find((e) => e.role === 'tools');
    expect(definitions?.title).toBe('2 tool definition messages');
    expect(definitions?.text).toContain('late');
  });
});

describe('readTrace: where the provider’s cache breaks (harness.md §11)', () => {
  const breaks = (t: TurnView) => models(t).map((m) => m.prefixBreak?.reason);

  it('says nothing for requests that extend each other, whatever the adapter reports', async () => {
    const s = await basic();
    s.model.push({ text: 'Again.' });
    await turn(s, 'Request: again');
    const trace = await readTrace(s.client);
    expect(trace.turns.flatMap(breaks)).toEqual([undefined, undefined, undefined]);
    // Without segments from the adapter there is nothing to compare.
    const t = await createSession(fixture, { declare: [], model: [{ text: 'a' }, { text: 'b' }] });
    await t.send(first(t, request('one')));
    await t.send('Request: two');
    expect((await readTrace(t.client)).turns.flatMap((x) => models(x).map((m) => m.prefixBreak))).toEqual([
      undefined,
      undefined,
    ]);
  });

  it('a changed tool list breaks it at the tools', async () => {
    const s = await createSession(fixture, {
      declare: [],
      tools: {
        disclose: () => ({ text: 'ok', contributions: [declareTools([toolDefinition('late')])] }),
        late: () => ({ text: 'late' }),
      },
      model: [{ calls: [{ tool: 'disclose', arguments: {} }] }, { text: 'one' }, { text: 'two' }],
    });
    await turn(s, first(s, request('first'), { deferred: ['late'] }));
    await turn(s, request('second'));
    const [t1, t2] = (await readTrace(s.client)).turns as [TurnView, TurnView];
    expect(breaks(t1)).toEqual([undefined, 'its tools changed']);
    expect(models(t1)[1]?.prefixBreak).toMatchObject({ kept: 0, of: 3 });
    expect(breaks(t2)).toEqual([undefined]);
  });

  it('a changed system text breaks it after the tools', async () => {
    adapterState.systemVersion = 0;
    const s = await createSession(fixture, { declare: [], model: [{ text: 'one' }, { text: 'two' }] });
    await turn(s, first(s, request('first')));
    adapterState.systemVersion = 1;
    await turn(s, 'Request: second');
    adapterState.systemVersion = 0;
    const [t1, t2] = (await readTrace(s.client)).turns as [TurnView, TurnView];
    expect(breaks(t1)).toEqual([undefined]);
    expect(breaks(t2)).toEqual(['its system text changed']);
    expect(models(t2)[0]?.prefixBreak).toMatchObject({ kept: 1, of: 3 });
  });

  it('an edit to the conversation (a compaction) breaks it at the first message that changed', async () => {
    const s = await createSession(fixture, {
      declare: [],
      model: [{ text: 'one' }, { text: 'two' }],
      drive: with_({
        'next-turn': ({ payload }, ctx) =>
          payload.size.inputTokens > 0 && ctx.inbox.length > 0
            ? {
                begin: {
                  keys: ctx.inbox.map((item) => item.key),
                  compact: { entries: [...ctx.head.entries.slice(0, 2), user('Summary of the work so far.')] },
                },
              }
            : 'pass',
      }),
    });
    await turn(s, first(s, request('first')));
    await turn(s, 'Request: second');
    const trace = await readTrace(s.client);
    const [t1, t2] = trace.turns as [TurnView, TurnView];
    expect(breaks(t1)).toEqual([undefined]);
    // Guidance and tools kept; the first message after them (the request) was replaced by a summary.
    expect(breaks(t2)).toEqual(['the conversation changed at message 1 (an edit or a compaction)']);
    expect(models(t2)[0]?.prefixBreak).toMatchObject({ kept: 2 });
    expect(t2.milestones.some((m) => m.detail?.startsWith('compacted main'))).toBe(true);
  });

  it('a changed model or setting breaks it whole', async () => {
    const s = await createSession(fixture, {
      declare: [],
      model: [{ text: 'one' }, { text: 'two' }, { text: 'three' }],
    });
    await turn(s, first(s, request('first')));
    await turn(s, 'Request: second', configured(s, { parameters: { temperature: 0.2 } }));
    await turn(s, 'Request: third', configured(s, { target: { provider: 'scripted', model: 'scripted-2' } }));
    const trace = await readTrace(s.client);
    expect(trace.turns.flatMap(breaks)).toEqual([
      undefined,
      "the model's settings changed",
      'another model (scripted-2)',
    ]);
    expect(models(trace.turns[2] as TurnView)[0]).toMatchObject({
      model: 'scripted-2',
      prefixBreak: { kept: 0, of: 5 },
    });
  });
});

describe('readTrace: a turn that is still going (harness.md §4, §6)', () => {
  it('is thinking while a model step is pending, and working while its tools run', async () => {
    const s = await createSession(fixture, {
      declare: [],
      tools: {
        slow: (context) =>
          context.signal.aborted
            ? { text: 'x' }
            : s.clock.sleep(60_000, context.signal).then(() => ({ text: 'slow done' })),
      },
      model: [{ calls: [{ tool: 'slow', arguments: {} }] }, { text: 'wait', waitMs: 60_000 }, { text: 'end' }],
    });
    const { launch } = await s.send(first(s, request('go')), { run: false });
    const running = runnerOf(s).run(launch as string);
    await until(() => s.model.invocations.length === 1 && s.tools.callsOf('slow').length === 1);
    const working = await readTrace(s.client, { now: () => s.clock.now() });
    expect(working.turns[0]?.status).toBe('working');
    expect(working.execution).toMatchObject({ runner: 'live', pendingStep: 'tools' });
    expect(working.execution.lease?.leaseMs).toBeGreaterThan(0);
    expect(working.turns[0]?.steps.map((x) => x.kind)).toEqual(['model']);
    s.clock.advance(60_000).catch(() => undefined);
    await until(() => s.model.invocations.length === 2);
    const thinking = await readTrace(s.client);
    expect(thinking.turns[0]?.status).toBe('thinking');
    expect(thinking.execution.pendingStep).toBe('model');
    expect(thinking.turns[0]?.end).toBeUndefined();
    expect(thinking.turns[0]?.steps.map((x) => x.kind)).toEqual(['model', 'tool']);
    await s.clock.runUntil(running);
    expect((await readTrace(s.client)).turns[0]?.status).toBe('done');
  });

  it('is waiting when it is parked for an answer', async () => {
    const s = await createSession(fixture, {
      declare: [],
      tools: { ask: () => ({ text: 'asked' }) },
      model: [{ calls: [{ tool: 'ask', arguments: {} }] }],
      configuration: { tools: { enabled: ['ask'], preauthorized: [] } },
      drive: with_({ permission: () => ({ park: {} }) }),
    });
    await turn(s, first(s, request('go')));
    const trace = await readTrace(s.client);
    expect(trace.turns[0]?.status).toBe('waiting');
    expect(trace.execution).toMatchObject({ parked: true, pendingStep: 'tools' });
  });

  it('is stopped when it was stopped', async () => {
    const s = await createSession(fixture, { declare: [], model: [{ text: 'slow', waitMs: 600_000 }] });
    const { launch } = await s.send(first(s, request('go')), { run: false });
    const running = runnerOf(s).run(launch as string);
    await until(() => s.model.invocations.length === 1);
    const turnId = (await s.client.execution()).turn?.id as string;
    await s.client.stop({ key: 'stop-1', turnId });
    await s.clock.runUntil(running);
    const trace = await readTrace(s.client);
    expect(trace.turns[0]?.status).toBe('stopped');
    expect(trace.turns[0]?.failure).toBeUndefined();
  });
});

describe('readTrace: the budget (harness.md §11)', () => {
  it('shows the kind of count nearest its ceiling, from the profile or from the limits the turn ran under', async () => {
    const s = await basic();
    const named = await readTrace(s.client, { ceilings: { uncachedInput: 1000, output: 30 } });
    // 300 of 1000 uncached input, 24 of 30 output.
    expect(named.budget).toEqual({
      kind: 'output',
      consumed: 24,
      ceiling: 30,
      dimension: 'output',
      scope: 'latest turn',
    });
    const unlimited = await readTrace(s.client);
    expect(unlimited.budget).toEqual({
      kind: 'uncachedInput',
      consumed: 0,
      ceiling: 0,
      dimension: 'uncached input',
      scope: 'latest turn',
    });

    const limited = await createSession(fixture, {
      declare: [],
      model: [{ text: 'ok', usage: { inputTokens: 500, outputTokens: 10 } }],
      configuration: {
        model: { target: { provider: 'scripted', model: 'scripted-1' }, parameters: {} },
        tools: { enabled: [], preauthorized: [] },
        limits: { uncachedInput: 1000, output: 100, modelCalls: 10 },
      },
    });
    await turn(limited, first(limited, request('hi')));
    expect((await readTrace(limited.client)).budget).toMatchObject({
      kind: 'uncachedInput',
      consumed: 500,
      ceiling: 1000,
    });
  });
});

describe('readTrace: the profile (the host’s part)', () => {
  it('reads the request with the host’s own marker', async () => {
    const s = await createSession(fixture, { declare: [], model: [{ text: 'ok' }] });
    await turn(s, first(s, 'Question: why?\nWhere: here'));
    const trace = await readTrace(s.client, {
      requestOf: (m) => {
        const text = m.parts.map((p) => (p.kind === 'text' ? p.text : '')).join('');
        return { prompt: /Question: (.*)/.exec(text)?.[1] ?? '', focus: /Where: (.*)/.exec(text)?.[1] ?? '' };
      },
    });
    expect(trace.turns[0]).toMatchObject({ prompt: 'why?', focus: 'here' });
  });

  it('classifies the window with the host’s own parts', async () => {
    const s = await basic();
    const trace = await readTrace(s.client, { windowPart: (m) => (m.role === 'user' ? 'state' : 'guidance') });
    const last = models(trace.turns[0] as TurnView)[1] as ModelStepView;
    expect(last.window?.entries.map((e) => e.part)).toEqual(['guidance', 'guidance', 'state', 'guidance', 'state']);
  });

  it('takes a tool call’s title, whether it read, and what it changed from the host', async () => {
    const s = await basic();
    const seen: string[] = [];
    const trace = await readTrace(s.client, {
      toolNames: new Map([['echo', 'Echo']]),
      toolCall: (call) => {
        seen.push(`${call.toolId}:${call.name}:${call.status}:${call.text}`);
        return call.toolId === 'echo'
          ? { title: 'Echo a', read: true, changed: [{ id: 'scene-1', label: '§1' }] }
          : undefined;
      },
    });
    expect(seen).toEqual(['echo:Echo:succeeded:echoed', 'boom:boom:failed:it broke']);
    const [echo, boom] = tools(trace.turns[0] as TurnView) as [ToolStepView, ToolStepView];
    expect(echo).toMatchObject({
      name: 'Echo',
      title: 'Echo a',
      read: true,
      changed: [{ id: 'scene-1', label: '§1' }],
    });
    expect(boom.title).toBeUndefined();
    expect(boom.read).toBeUndefined();
    expect(boom.changed).toBeUndefined();
    expect(models(trace.turns[0] as TurnView)[0]?.toolCalls).toEqual(['Echo(a=1)', 'boom()']);
  });

  it('prices model calls with the host’s credits, and sums them for the turn and the budget', async () => {
    const s = await basic();
    const trace = await readTrace(s.client, {
      credits: (model, usage) =>
        model === 'scripted-1' ? (usage.uncachedInput + usage.cacheRead * 0.1) / 100 : undefined,
      ceilings: { output: 100 },
    });
    const [t] = trace.turns as [TurnView];
    const [a, b] = models(t) as [ModelStepView, ModelStepView];
    expect(a.tokens?.credits).toBeCloseTo(2.8);
    expect(b.tokens?.credits).toBeCloseTo(2);
    expect(t.tokens.credits).toBeCloseTo(4.8);
    expect(trace.budget.credits).toBeCloseTo(4.8);
    expect(formatUsage(t.tokens)).toBe('2.2k in (82% cached, 100 written) · 24 out · 4.8 credits');
    const unpriced = await readTrace(s.client);
    expect(unpriced.turns[0]?.tokens.credits).toBeUndefined();
  });

  it('keeps what it read in a cache, so a re-read fetches only what is new', async () => {
    const s = await basic();
    const counts = { entries: 0, heads: 0, blobs: 0 };
    const reader: HarnessReader = {
      ...s.client,
      entries: (ids) => (counts.entries++, s.client.entries(ids)),
      heads: (ids) => (counts.heads++, s.client.heads(ids)),
      blobs: (refs) => (counts.blobs++, s.client.blobs(refs)),
    };
    const cache = { entries: new Map(), heads: new Map(), blobs: new Map() };
    const one = await readTrace(reader, { cache });
    const fetched = { ...counts };
    expect(fetched.entries).toBeGreaterThan(0);
    const two = await readTrace(reader, { cache });
    expect(counts).toEqual(fetched);
    expect(two).toEqual(one);
  });
});
