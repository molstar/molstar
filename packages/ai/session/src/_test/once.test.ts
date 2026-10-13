/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type ModelStep, scriptedModel } from '@molstar/ai-conformance/testing';
import { createToolkit } from '@molstar/ai-tools';
import { runOnce, Settlement, terminalToolFor } from '../index.js';
import { callOf, modelsOf, realClock, tool, toy } from './helpers.js';

const terminal = terminalToolFor<ReturnType<typeof toy>>();

/** The one-off setting: a tool that reads, and a terminal one that proposes a number (the model may be asked again). */
function setup(script: ModelStep[], clock = false) {
  const settlement = new Settlement<{ value: number }>();
  const api = toy();
  const read = tool({
    id: 'get_state',
    description: 'Reads the value.',
    input: { kind: 'object', fields: {} },
    modules: ['state'],
    read: true,
    run: ({ api }) => ({ ok: true, did: `Value ${api.state.value}.` }),
  });
  const propose = terminal(settlement, {
    id: 'propose',
    description: 'Propose a value. Nothing runs until the person confirms.',
    input: { kind: 'object', fields: { value: { kind: 'number' } } },
    modules: [],
    run: ({ args }) =>
      args.value < 0
        ? { ok: false, code: 'propose.value', message: 'The value must not be negative.' }
        : { ok: true, value: { value: args.value }, did: 'Proposed.' },
  });
  const toolkit = createToolkit<typeof api>({ tools: [read, propose], core: ['get_state', 'propose'], groups: [] });
  const model = scriptedModel(script, clock ? { clock: realClock } : {});
  const run = (extra: Partial<Parameters<typeof runOnce<typeof api, { value: number }>>[0]> = {}) =>
    runOnce({
      models: modelsOf(model.adapter),
      guidance: 'You propose a value.',
      toolkit,
      api,
      prompt: 'Propose three.',
      settlement,
      ...extra,
    });
  return { run, model, api, settlement };
}

describe('one-off runs (tools-and-session.md §7)', () => {
  it('a terminal call settles the run with its value, and the run stops there', async () => {
    const { run, model, api } = setup([{ calls: [callOf('propose', { value: 3 })] }, { text: 'Never said.' }]);
    const r = await run();
    expect(r).toEqual({ status: 'settled', value: { value: 3 }, text: '' });
    // One model call: the proposal ends the run; the host's state is untouched.
    expect(model.invocations).toHaveLength(1);
    expect(model.remaining()).toBe(1);
    expect(api.state.value).toBe(0);
    // All its tools are offered from the start.
    const tools = (model.invocations[0]?.payload as unknown as { tools: string[] }).tools;
    expect(tools).toEqual(expect.arrayContaining(['get_state', 'propose']));
  });

  it('a proposal that does not fit goes back to the model within the run, which fixes it', async () => {
    const { run, model } = setup([
      { calls: [callOf('propose', { value: -1 })] },
      { calls: [callOf('get_state'), callOf('propose', { value: 2 })] },
    ]);
    const r = await run();
    expect(r.status).toBe('settled');
    expect(r.value).toEqual({ value: 2 });
    expect(model.invocations).toHaveLength(2);
  });

  it('the first terminal call settles; the value of a later one in the same response changes nothing', async () => {
    const { run } = setup([{ calls: [callOf('propose', { value: 1 }), callOf('propose', { value: 2 })] }]);
    expect((await run()).value).toEqual({ value: 1 });
  });

  it('an answer in words is the status completed, with its text', async () => {
    const { run } = setup([{ text: 'I cannot tell which value you mean.' }]);
    expect(await run()).toEqual({ status: 'completed', text: 'I cannot tell which value you mean.' });
  });

  it('a model that fails ends the run as failed, with its error', async () => {
    const { run } = setup([{ fail: { message: 'The service is down.', transient: false } }]);
    const r = await run();
    expect(r.status).toBe('failed');
    expect(r.error).toContain('The service is down.');
  });

  it('gives up after the timeout', async () => {
    const { run } = setup([{ text: 'Too late.', waitMs: 10_000 }], true);
    const started = Date.now();
    const r = await run({ timeoutMs: 150 });
    expect(r.status).toBe('timeout');
    expect(Date.now() - started).toBeLessThan(5_000);
  });

  it('stops when the signal aborts', async () => {
    const { run } = setup([{ text: 'Too late.', waitMs: 10_000 }], true);
    const controller = new AbortController();
    const running = run({ signal: controller.signal });
    setTimeout(() => controller.abort(), 100);
    expect((await running).status).toBe('aborted');
    const { run: again } = setup([{ text: 'Never.' }]);
    expect((await again({ signal: AbortSignal.abort() })).status).toBe('aborted');
  });

  it('a settlement is cleared for the next run', async () => {
    const { run, model, settlement } = setup([{ calls: [callOf('propose', { value: 3 })] }]);
    expect((await run()).status).toBe('settled');
    model.push({ text: 'Not this time.' });
    const second = await run();
    expect(second.status).toBe('completed');
    expect(settlement.settled).toBe(false);
  });
});
