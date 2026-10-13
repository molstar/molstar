/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { createMemoryStorage } from '@molstar/ai-components';
import type { Decisions } from '@molstar/ai-harness';
import { AgentSession, createDrive, DENIED, type DriveHost, extendBudget, plainFailure, turnBudget } from '../index.js';
import {
  callOf,
  fixture,
  guidance,
  headOf,
  messagesOf,
  resultsIn,
  tool,
  type Toy,
  toy,
  toyToolkit,
} from './helpers.js';

const outcomes = (r: { tools: { name: string; status: string }[] }) => r.tools.map((t) => `${t.name}:${t.status}`);

describe('a tool call that fails transiently (tools-and-session.md §5)', () => {
  it('is retried when the attempt changed nothing; the model sees nothing of it', async () => {
    const f = fixture([{ calls: [callOf('fetch_read')] }, { text: 'Done.' }]);
    f.api.flaky.failures = 2;
    const r = await f.session.turn({ prompt: 'fetch' });
    expect(r.status).toBe('completed');
    expect(f.api.flaky.calls).toBe(3);
    expect(outcomes(r)).toEqual(['fetch_read:succeeded']);
    expect(resultsIn(await headOf(f.session))).toEqual(['Fetched.']);
  });

  it('gives up after the retries and tells the model; a write that changed nothing is retried too', async () => {
    const f = fixture([{ calls: [callOf('fetch_write')] }, { text: 'It is down.' }], { toolRetries: 2 });
    f.api.flaky.failures = 10;
    const r = await f.session.turn({ prompt: 'write' });
    expect(f.api.flaky.calls).toBe(3);
    expect(outcomes(r)).toEqual(['fetch_write:failed']);
    expect(r.tools[0]?.text).toContain('unavailable');
  });

  it('is not retried when it may have changed something, unless the tool is safe to repeat', async () => {
    // The attempt applied a change (or it is not known): running it again could repeat it.
    for (const error of ['applied', 'throw'] as const) {
      const f = fixture([{ calls: [callOf('fetch_write')] }, { text: 'Failed.' }]);
      f.api.flaky.failures = 1;
      f.api.flaky.error = error;
      const r = await f.session.turn({ prompt: 'write' });
      expect(f.api.flaky.calls).toBe(1);
      expect(outcomes(r)).toEqual(['fetch_write:failed']);
    }
    // The same failure of a tool marked idempotent is retried.
    const safe = fixture([{ calls: [callOf('fetch_safe')] }, { text: 'Done.' }]);
    safe.api.flaky.failures = 1;
    safe.api.flaky.error = 'applied';
    const r = await safe.session.turn({ prompt: 'write' });
    expect(safe.api.flaky.calls).toBe(2);
    expect(outcomes(r)).toEqual(['fetch_safe:succeeded']);
  });

  it('a failed write skips the calls after it in that response; a failed read does not', async () => {
    const write = fixture([{ calls: [callOf('fetch_write'), callOf('set_value')] }, { text: 'Failed.' }]);
    write.api.flaky.failures = 1;
    write.api.flaky.error = 'applied';
    const w = await write.session.turn({ prompt: 'write, then bump' });
    expect(outcomes(w)).toEqual(['fetch_write:failed', 'set_value:skipped']);
    expect(w.tools[1]?.text).toBe('Not run: fetch_write failed before it.');
    expect(write.api.state.value).toBe(0);
    expect(resultsIn(await headOf(write.session))[1]).toContain('Tool execution skipped: fetch_write failed');

    const read = fixture([{ calls: [callOf('fetch_read'), callOf('set_value')] }, { text: 'Done.' }]);
    read.api.flaky.failures = 10;
    read.api.flaky.error = 'applied';
    const r = await read.session.turn({ prompt: 'read, then bump' });
    expect(outcomes(r)).toEqual(['fetch_read:failed', 'set_value:succeeded']);
    expect(read.api.state.value).toBe(1);
  });

  it('a validation error is not retried', async () => {
    const f = fixture([{ calls: [callOf('set_value', { by: 'many' })] }, { text: 'Oops.' }]);
    const r = await f.session.turn({ prompt: 'bump' });
    expect(outcomes(r)).toEqual(['set_value:failed']);
  });
});

describe('permission', () => {
  it('a denial is told to the model and the turn goes on; the person sees a failed call', async () => {
    const f = fixture([{ calls: [callOf('ask_service')] }, { text: 'Not allowed, so I did not.' }], {
      approve: () => 'deny',
    });
    const r = await f.session.turn({ prompt: 'ask the service' });
    expect(r.status).toBe('completed');
    expect(r.text).toBe('Not allowed, so I did not.');
    expect(r.tools).toEqual([
      expect.objectContaining({ name: 'ask_service', status: 'failed', text: DENIED, read: false }),
    ]);
    expect(f.api.log.lines).toEqual([]);
    expect(resultsIn(await headOf(f.session))).toEqual([`Denied: ${DENIED}`]);
  });

  it('a rule allows a call without asking', async () => {
    const f = fixture([{ calls: [callOf('ask_service')] }, { text: 'Done.' }], { approve: () => 'allow' });
    const r = await f.session.turn({ prompt: 'ask' });
    expect(outcomes(r)).toEqual(['ask_service:succeeded']);
  });

  it('parks at the decision, and the answers given on resume decide it', async () => {
    for (const answer of ['allow', 'deny'] as const) {
      const f = fixture([{ calls: [callOf('ask_service')] }, { text: 'Finished.' }]);
      const statuses: string[] = [];
      const parked = await f.session.turn({ prompt: 'ask' }, { onStatus: (s) => statuses.push(s.kind) });
      expect(parked.status).toBe('awaiting-approval');
      expect(parked.approval).toEqual([expect.objectContaining({ tool: 'ask_service', args: {} })]);
      expect(statuses).toContain('permission');
      expect(await f.session.interrupted()).toBe(true);
      const callId = parked.approval?.[0]?.callId ?? '';
      const done = await f.session.resume({}, { [callId]: answer });
      expect(done?.status).toBe('completed');
      expect(outcomes(done as never)).toEqual([`ask_service:${answer === 'allow' ? 'succeeded' : 'failed'}`]);
      expect(f.api.log.lines).toEqual(answer === 'allow' ? ['ask_service'] : []);
      expect(await f.session.resume()).toBeUndefined();
    }
  });

  it('the host’s own decision comes first; its pass leaves it to the session’s', async () => {
    const f = fixture([{ calls: [callOf('ask_service')] }, { text: 'Done.' }], {
      drive: { permission: () => 'allow' },
    });
    expect((await f.session.turn({ prompt: 'ask' })).tools[0]?.status).toBe('succeeded');
    const passing = fixture([{ calls: [callOf('ask_service')] }, { text: 'Done.' }], {
      drive: { permission: () => 'pass' },
      approve: () => 'deny',
    });
    expect((await passing.session.turn({ prompt: 'ask' })).tools[0]?.status).toBe('failed');
  });
});

describe('a call the model should not have made (repair)', () => {
  it('a tool not loaded yet: answered failed and loaded, the rest of the response skipped', async () => {
    const f = fixture([
      { calls: [callOf('load_map'), callOf('set_value')] },
      { calls: [callOf('load_map')] },
      { text: 'Loaded the map.' },
    ]);
    const r = await f.session.turn({ prompt: 'load a map' });
    expect(outcomes(r)).toEqual(['load_map:failed', 'set_value:skipped', 'load_map:succeeded']);
    expect(r.tools[0]?.text).toBe(
      'load_map was not loaded; it is now. Call it again. It belongs to the maps skill: disclose_skills(["maps"]) gives its guidance.',
    );
    expect(r.tools[1]?.text).toContain('Not run: load_map was not loaded yet');
    expect(f.api.state.value).toBe(0);
    expect(f.api.log.lines).toEqual(['load_map']);
  });

  it('an unknown tool points to the disclosure tools; invalid arguments say how they come', async () => {
    const f = fixture([
      { calls: [callOf('fly'), { tool: 'set_value', arguments: 'by two' }, callOf('set_value')] },
      { text: 'Ok.' },
    ]);
    const r = await f.session.turn({ prompt: 'go' });
    expect(r.tools.map((t) => t.text)).toEqual([
      'There is no tool "fly". Load tools with disclose_skills or disclose_tools; their descriptions list what there is.',
      'set_value takes its arguments as one JSON object.',
      'Bumped.\nChanged: value.',
    ]);
    // Nothing is skipped: the call that was fine ran.
    expect(outcomes(r)).toEqual(['fly:failed', 'set_value:failed', 'set_value:succeeded']);
  });

  it('a tool that cannot run now (loaded earlier, no longer enabled): said so, nothing loaded, nothing skipped', async () => {
    // One release had a tool the next one has not: its definition is still in the conversation's head.
    const storage = createMemoryStorage();
    const legacy = tool({
      id: 'legacy',
      description: 'An old tool.',
      input: { kind: 'object', fields: {} },
      modules: [],
      run: () => ({ ok: true, did: 'Legacy.' }),
    });
    const before = fixture([{ calls: [callOf('disclose_tools', { tools: ['legacy'] })] }, { text: 'One.' }], {
      storage,
      harnessId: 'h-legacy',
      toolkit: toyToolkit({ tools: [legacy] }),
    });
    await before.session.turn({ prompt: 'first' });
    const f = fixture(
      [{ calls: [callOf('legacy'), callOf('set_value')] }, { text: 'Ok.' }],
      { storage, harnessId: 'h-legacy' },
      {},
      before.api,
    );
    const r = await f.session.turn({ prompt: 'second' });
    expect(outcomes(r)).toEqual(['legacy:failed', 'set_value:succeeded']);
    expect(r.tools[0]?.text).toBe('legacy cannot run now.');
    expect(r.status).toBe('completed');
    expect(f.model.invocations).toHaveLength(2);
  });
});

/** A drive over a stand-in host, for the decisions a session cannot be made to reach. */
function drive() {
  const toolkit = toyToolkit();
  const shown: unknown[] = [];
  const host: DriveHost<never> = {
    toolkit: toolkit as never,
    api: {} as never,
    policy: {},
    decided: new Map(),
    turn: () => ({}) as never,
    registration: (id) => toolkit.tools.find((t) => t.definition.tool.id === id) as never,
    budgetFor: () => turnBudget(100),
    outputCap: () => undefined,
    compaction: async () => undefined,
    status: () => undefined,
    shown: (e) => void shown.push(e),
  };
  const context = { depth: 0, inbox: [] } as never;
  return {
    shown,
    repair: (reason: Decisions['repair'][0]['reason'], toolId = 'set_value') =>
      createDrive(host).decide(
        {
          kind: 'repair',
          path: [],
          payload: { call: { callId: 'c1', toolId, arguments: {} }, reason, unstarted: [] },
        },
        context,
      ),
  };
}

describe('the repair decision for each reason', () => {
  it('not-enabled and no-handler say it cannot run now, load nothing and skip nothing', async () => {
    const d = drive();
    for (const reason of ['not-enabled', 'no-handler'] as const) {
      const answer = await d.repair(reason);
      expect(answer).toEqual({ report: { parts: [{ kind: 'text', text: 'set_value cannot run now.' }] } });
    }
    expect(d.shown).toHaveLength(2);
  });

  it('not-loaded loads the definition and skips the rest; unknown-tool and invalid-arguments load nothing', async () => {
    const d = drive();
    const loaded = (await d.repair('not-loaded', 'archive')) as {
      report: { contributions?: unknown[]; skip?: { text: string } };
    };
    expect(loaded.report.contributions).toHaveLength(1);
    expect(loaded.report.skip?.text).toContain('archive was not loaded yet');
    for (const reason of ['unknown-tool', 'invalid-arguments'] as const) {
      const answer = (await d.repair(reason)) as { report: { contributions?: unknown; skip?: unknown } };
      expect(answer.report.contributions).toBeUndefined();
      expect(answer.report.skip).toBeUndefined();
    }
  });
});

describe('the budget (tools-and-session.md §5)', () => {
  const ZERO = { modelCalls: 0, toolCalls: 0, uncachedInput: 0, cacheRead: 0, cacheWrite: 0, output: 0, reasoning: 0 };
  /** The budget decision when `toolCalls` is short by one call. */
  const short = (limit: number, consumed: number): Decisions['budget'][0] => ({
    dimension: 'toolCalls',
    budget: { limits: { toolCalls: limit }, consumed: { ...ZERO, toolCalls: consumed } },
    need: 1,
  });
  const policy = { initial: { toolCalls: 60 }, ceiling: { toolCalls: 200 } };

  it('extends a short budget by its initial amount, up to the ceiling', () => {
    expect(extendBudget(short(60, 60), policy)).toEqual({ toolCalls: 120 });
    expect(extendBudget(short(180, 180), policy)).toEqual({ toolCalls: 200 });
    expect(extendBudget(short(200, 200), policy)).toBeUndefined();
  });

  it('extends each kind of token on its own: cache reads past theirs, the rest unchanged', () => {
    const cacheShort: Decisions['budget'][0] = {
      dimension: 'cacheRead',
      budget: { limits: { cacheRead: 5_000_000, uncachedInput: 500_000 }, consumed: { ...ZERO, cacheRead: 4_990_000 } },
      need: 20_000,
    };
    const split = {
      initial: { cacheRead: 5_000_000, uncachedInput: 500_000 },
      ceiling: { cacheRead: 25_000_000, uncachedInput: 2_000_000 },
    };
    expect(extendBudget(cacheShort, split)).toEqual({ cacheRead: 10_000_000 });
  });

  it('scales with the compaction trigger, the same for every model', () => {
    expect(turnBudget(160_000)).toEqual({
      initial: { toolCalls: 60, uncachedInput: 640_000, cacheRead: 6_400_000, output: 200_000 },
      ceiling: { toolCalls: 200, uncachedInput: 3_200_000, cacheRead: 32_000_000, output: 1_000_000 },
    });
    expect(turnBudget(1_000, 5).initial.toolCalls).toBe(5);
  });

  it('extends a turn that used its tool calls, and fails one at its ceiling', async () => {
    const calls = Array.from({ length: 61 }, () => ({ calls: [callOf('get_state')] }));
    const extended = fixture([...calls, { text: 'Done.' }]);
    const statuses: { kind: string; detail?: string }[] = [];
    const r = await extended.session.turn({ prompt: 'read a lot' }, { onStatus: (s) => statuses.push(s) });
    expect(r.status).toBe('completed');
    expect(statuses.filter((s) => s.kind === 'budget-extended')).toEqual([
      { kind: 'budget-extended', detail: 'toolCalls → 120' },
    ]);

    const capped = fixture([...calls.slice(0, 5), { text: 'Done.' }], { maxToolCalls: 3 });
    const c = await capped.session.turn({ prompt: 'read a lot' });
    expect(c.status).toBe('failed');
    expect(c.error).toBe('The turn reached its budget.');
  });

  it('asks once more with twice the room when a reply is cut short at the output limit', async () => {
    const f = fixture([{ text: 'Let me think about this very long', finish: 'output-limit' }, { text: 'Done.' }], {
      maxOutputTokens: 1000,
    });
    const statuses: string[] = [];
    const r = await f.session.turn({ prompt: 'Hi.' }, { onStatus: (s) => statuses.push(s.kind) });
    expect(r.status).toBe('completed');
    expect(r.text).toBe('Done.');
    expect(f.model.projections.map((p) => p.configuration.parameters.maxOutputTokens)).toEqual([1000, 2000]);
    expect(statuses.filter((s) => s === 'widened')).toHaveLength(1);
    // Once per turn: a second cut-short reply is accepted as it is.
    const again = fixture(
      [
        { text: 'a', finish: 'output-limit' },
        { text: 'b', finish: 'output-limit' },
      ],
      { maxOutputTokens: 1000 },
    );
    const second = await again.session.turn({ prompt: 'Hi.' });
    expect(again.model.projections.map((p) => p.configuration.parameters.maxOutputTokens)).toEqual([1000, 2000]);
    expect(second.text).toBe('b');
  });
});

describe('a failing model call', () => {
  const names = { providerNames: { scripted: 'Scripted' } };
  const failing = (message: string, times = 1, transient = true) => {
    const steps = Array.from({ length: times }, () => ({ fail: { message, transient } }));
    return fixture([...steps, { text: 'Done.' }], names);
  };

  it('is retried after a wait, and the turn goes on', async () => {
    const f = failing('Your credit balance is too low.');
    const statuses: { kind: string; detail?: string }[] = [];
    const r = await f.session.turn({ prompt: 'Hi.' }, { onStatus: (s) => statuses.push(s) });
    expect(r).toMatchObject({ status: 'completed', text: 'Done.' });
    expect(statuses).toContainEqual({ kind: 'retrying', detail: 'Your credit balance is too low.' });
  });

  it('ends the turn with the model’s own error once the retries run out', async () => {
    const r = await failing('Your credit balance is too low.', 10).session.turn({ prompt: 'Hi.' });
    expect(r.status).toBe('failed');
    expect(r.error).toContain('credit balance is too low');
  });

  it('is not retried when it cannot succeed (a missing model), and says what to do', async () => {
    const f = failing('model: scripted (HTTP 404)', 10, false);
    const statuses: string[] = [];
    const r = await f.session.turn({ prompt: 'Hi.' }, { onStatus: (s) => statuses.push(s.kind) });
    expect(r.status).toBe('failed');
    expect(r.error).toBe('No model "scripted-1" at Scripted: check the model id in the options.');
    expect(f.model.invocations).toHaveLength(1);
    expect(statuses).not.toContain('retrying');
  });

  it('is not retried when the provider asks for a wait no turn can make (a daily limit)', async () => {
    const f = failing('Rate limit reached. Please try again in 63h29m22.56s. (HTTP 429)', 10);
    const r = await f.session.turn({ prompt: 'Hi.' });
    expect(r.status).toBe('failed');
    expect(r.error).toMatch(/^Scripted asks to wait 63 hours before the next call/);
    expect(f.model.invocations).toHaveLength(1);
  });

  it('reads milliseconds as milliseconds: a wait of 337 ms is retried, not 337 minutes', async () => {
    const r = await failing('Rate limit reached. Please try again in 337ms. (HTTP 429)').session.turn({
      prompt: 'Hi.',
    });
    expect(r).toMatchObject({ status: 'completed', text: 'Done.' });
  });

  it('says it in plain words, with the host’s provider names', () => {
    const model = { target: { provider: 'acme', model: 'm1' }, parameters: {} };
    expect(plainFailure({ code: 'provider', message: 'x (HTTP 401)' }, model, { acme: 'Acme' })).toBe(
      'Acme refused the key (x (HTTP 401)).',
    );
    expect(
      plainFailure(
        { code: 'provider', message: 'x (HTTP 403)' },
        { ...model, target: { provider: 'google', model: 'm' } },
      ),
    ).toBe('Google refused the key (x (HTTP 403)).');
    expect(plainFailure({ code: 'capacity', message: '' })).toContain('no longer fits');
  });
});

describe('the model per turn', () => {
  it('sizes each model by its own limits: the application’s window, the model’s output maximum', async () => {
    const f = fixture([{ text: 'One.' }, { text: 'Two.' }, { text: 'Three.' }], {}, {});
    // An adapter that knows one model and refuses the others.
    const inner = f.model.adapter;
    const seen: { model: string; window: number; output: number | undefined; reasoning: unknown }[] = [];
    const adapter = {
      ...inner,
      capabilitiesFor: async (c: { target: { model: string } }) => {
        if (c.target.model !== 'scripted-1') throw new Error(`Unknown model ${c.target.model}.`);
        return { ...(await inner.capabilitiesFor(c as never)), contextWindowTokens: 1_114_112, maxOutputTokens: 8192 };
      },
      project: async (input: Parameters<typeof inner.project>[0]) => {
        seen.push({
          model: input.configuration.target.model,
          window: input.capabilities.contextWindowTokens,
          output: input.configuration.parameters.maxOutputTokens,
          reasoning: input.configuration.parameters.reasoning,
        });
        return inner.project(input);
      },
    };
    const session = new AgentSession<Toy>({
      models: {
        adapterFor: () => adapter,
        default: { provider: 'google', model: 'scripted-1', contextWindow: 1_000_000 },
      },
      toolkit: toyToolkit(),
      api: toy(),
      guidance,
    });
    await session.turn({ prompt: 'One.' });
    // A model the adapter does not know runs with the window the application gives it.
    const custom = await session.turn({
      prompt: 'Two.',
      model: { provider: 'google', model: 'gemini-2.5-pro-custom', contextWindow: 200_000, maxOutput: 4096 },
    });
    expect(custom.status).toBe('completed');
    // Without a window it cannot be sized, and says so.
    const unknown = await session.turn({ prompt: 'Three.', model: { provider: 'google', model: 'gemini-9' } });
    expect(unknown.status).toBe('failed');
    expect(unknown.error).toContain('give its context window');
    expect(session.contextWindow()).toBe(1_000_000);
    expect(seen).toEqual([
      // The asked-for 16k clamped to the model's 8k; the application's window over the adapter's.
      { model: 'scripted-1', window: 1_000_000, output: 8192, reasoning: undefined },
      // No thinking asked for, none sent.
      { model: 'gemini-2.5-pro-custom', window: 200_000, output: 4096, reasoning: undefined },
    ]);
  });

  it('runs a turn on another provider’s model, with its own thinking', async () => {
    const make = (name: string) => {
      const m = fixture([{ text: `from ${name}` }, { text: `from ${name}` }]);
      return m.model;
    };
    const other = make('claude');
    const main = fixture([{ text: 'from main' }, { text: 'from main' }]);
    const session = new AgentSession<Toy>({
      models: {
        adapterFor: (p) => (p === 'scripted' ? main.model.adapter : p === 'anthropic' ? other.adapter : undefined),
        default: { provider: 'scripted', model: 'scripted-1' },
      },
      toolkit: toyToolkit(),
      api: toy(),
      guidance,
    });
    const first = await session.turn({ prompt: 'One.' });
    const second = await session.turn({
      prompt: 'Two.',
      model: { provider: 'anthropic', model: 'claude-x', thinking: { mode: 'adaptive', effort: 'high' } },
    });
    const third = await session.turn({ prompt: 'Three.' });
    expect([first.text, second.text, third.text]).toEqual(['from main', 'from claude', 'from main']);
    expect(other.projections.map((p) => p.configuration.parameters.reasoning)).toEqual([
      { mode: 'adaptive', effort: 'high' },
    ]);
    const unknown = await session.turn({ prompt: 'Four.', model: { provider: 'openai', model: 'x' } });
    expect(unknown).toMatchObject({ status: 'failed', error: 'No adapter for openai.' });
  });
  it('thinks per model as the host says, under the turn’s own choice', async () => {
    const main = fixture([{ text: 'a' }, { text: 'b' }, { text: 'c' }, { text: 'd' }]);
    const session = new AgentSession<Toy>({
      models: { adapterFor: () => main.model.adapter, default: { provider: 'scripted', model: 'fast' } },
      toolkit: toyToolkit(),
      api: toy(),
      guidance,
      // As a host passes the catalog's default (`defaultReasoning`).
      thinking: ({ model }) => (model === 'fast' ? { mode: 'adaptive', effort: 'low' } : undefined),
    });
    await session.turn({ prompt: 'One.' });
    await session.turn({ prompt: 'Two.', model: { provider: 'scripted', model: 'other' } });
    await session.turn({
      prompt: 'Three.',
      model: { provider: 'scripted', model: 'fast', thinking: { mode: 'adaptive', effort: 'high' } },
    });
    await session.turn({ prompt: 'Four.', model: { provider: 'scripted', model: 'fast', thinking: 'model' } });
    expect(main.model.projections.map((p) => p.configuration.parameters.reasoning)).toEqual([
      { mode: 'adaptive', effort: 'low' },
      undefined,
      { mode: 'adaptive', effort: 'high' },
      undefined,
    ]);
  });
});

describe('a conversation that began before a release', () => {
  it('learns tools added since, and guidance that changed, once', async () => {
    const storage = createMemoryStorage();
    const extra = tool({
      id: 'added_later',
      description: 'Came with a release.',
      input: { kind: 'object', fields: {} },
      modules: [],
      run: () => ({ ok: true, did: 'Added.' }),
    });
    const api = fixture().api;
    const session = (withExtra: boolean, guidanceText: string, changed = false) => {
      const base = toyToolkit(withExtra ? { tools: [extra] } : {});
      const toolkit = changed
        ? {
            ...base,
            allDefinitions: () =>
              base
                .allDefinitions()
                .map((d) => (d.name === 'get_state' ? { ...d, description: `${d.description} Changed.` } : d)),
            coreDefinitions: () =>
              base
                .coreDefinitions()
                .map((d) => (d.name === 'get_state' ? { ...d, description: `${d.description} Changed.` } : d)),
          }
        : base;
      return fixture(
        [{ text: 'Done.' }],
        { storage, harnessId: 'h-release', guidance: guidanceText, toolkit },
        {},
        api,
      );
    };
    const defs = (messages: Awaited<ReturnType<typeof headOf>>) =>
      messages.flatMap((m) => m.parts.flatMap((p) => (p.kind === 'tool-definition' ? [p.definition.name] : [])));
    const guidances = (messages: Awaited<ReturnType<typeof headOf>>) =>
      messages.flatMap((m) =>
        m.parts.flatMap((p) =>
          p.kind === 'text' && (m.role === 'system' || p.text.startsWith('Updated guidance')) ? [p.text] : [],
        ),
      );

    const one = session(false, 'Guidance one.');
    await one.session.turn({ prompt: 'first' });
    expect(defs(await headOf(one.session))).not.toContain('added_later');

    // A later release: one more tool (declared, not loaded), new guidance. The next turn adds exactly those.
    const two = session(true, 'Guidance two.');
    await two.session.turn({ prompt: 'second' });
    const after = await headOf(two.session);
    expect(defs(after).filter((n) => n === 'added_later')).toHaveLength(1);
    expect(guidances(after).at(-1)).toMatch(/^Updated guidance .*\n\nGuidance two\.$/s);

    // Nothing new: nothing added.
    const three = session(true, 'Guidance two.');
    await three.session.turn({ prompt: 'third' });
    const third = await headOf(three.session);
    expect(defs(third)).toHaveLength(defs(after).length);
    expect(guidances(third)).toHaveLength(2);

    // A changed definition is sent again, once.
    const four = session(true, 'Guidance two.', true);
    await four.session.turn({ prompt: 'fourth' });
    await session(true, 'Guidance two.', true).session.turn({ prompt: 'fifth' });
    expect(defs(await headOf(four.session)).filter((n) => n === 'get_state')).toHaveLength(2);
    expect(messagesOf(four.model, 0).at(-1)?.role).toBe('user');
  });
});
