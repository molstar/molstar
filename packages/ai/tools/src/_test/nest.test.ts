/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { ModelConfiguration } from '@molstar/ai-agent';
import type { RunResult, RunSpec } from '@molstar/ai-harness';
import { createToolBindings, NestedTools, toolFor } from '../index.js';
import { call, runContext, type Toy, toy } from './helpers.js';

interface Choice {
  name: string;
}
const tool = toolFor<Toy, Choice>();
const none = { kind: 'object', fields: {} } as const;
const done: RunResult = {
  status: 'completed',
  usage: { modelCalls: 1, toolCalls: 1, uncachedInput: 0, cacheRead: 0, cacheWrite: 0, output: 0, reasoning: 0 },
};

describe('nested tools (harness.md §9)', () => {
  const configure = (choice: Choice): ModelConfiguration => ({
    target: { provider: 'test', model: choice.name },
    parameters: {},
  });

  it('a tool runs a helper under its call: its tools, model and limits go to the run, and only its own call reaches them', async () => {
    const api = toy();
    const inner = tool({
      id: 'inner',
      description: 'The helper’s tool.',
      input: none,
      modules: ['state'],
      run: ({ api }) => {
        api.state.bump(5);
        return { ok: true, did: 'Inner done.' };
      },
    });
    const outer = tool({
      id: 'outer',
      description: 'Uses a helper.',
      input: none,
      modules: [],
      run: async ({ nest }) => {
        if (!nest) return { ok: false, code: 'no.nest', message: 'Called without a turn.' };
        const r = await nest({
          messages: [{ role: 'user', parts: [{ kind: 'text', text: 'Help.' }] }],
          tools: [inner],
          api,
          model: { name: 'small' },
          limits: { toolCalls: 3 },
        });
        return { ok: true, did: `Helper ${r.status}.` };
      },
    });
    const nested = new NestedTools<Toy, Choice>(configure);
    const bindings = createToolBindings([outer], api, { nested });
    const specs: RunSpec[] = [];
    const innerFrames: string[] = [];
    const context = runContext(
      outer.definition,
      {},
      {
        callId: 'outer-call',
        run: async (spec) => {
          specs.push(spec);
          // The harness would run the helper's calls; here one is run through the bindings, as the harness does,
          // with the nested call's path ending in the outer call.
          const handler = bindings.resolve(inner.definition);
          for await (const f of handler?.invoke(
            runContext(inner.definition, {}, { callId: 'inner-call', path: ['outer-call'] }),
          ) ?? [])
            if (f.kind === 'result') innerFrames.push(f.status);
          return done;
        },
      },
    );
    // Before the run, nothing resolves the helper's tool.
    expect(bindings.resolve(inner.definition)).toBeUndefined();
    const out: string[] = [];
    for await (const f of bindings.resolve(outer.definition)?.invoke(context) ?? [])
      if (f.kind === 'result') out.push(f.parts.map((p) => (p.kind === 'text' ? p.text : '')).join(''));
    expect(out).toEqual(['Helper completed.']);
    expect(innerFrames).toEqual(['succeeded']);
    expect(api.state.value).toBe(5);
    expect(specs).toHaveLength(1);
    expect(specs[0]).toMatchObject({
      tools: { enabled: ['inner'], preauthorized: ['inner'] },
      model: { target: { provider: 'test', model: 'small' } },
      limits: { toolCalls: 3 },
    });
    expect(specs[0]?.tools.definitions.map((d) => d.tool.id)).toEqual(['inner']);
    // After it, the helper's tool is gone.
    expect(bindings.resolve(inner.definition)).toBeUndefined();
  });

  it('a call to a helper’s tool from another call than the one that runs it is refused', async () => {
    const api = toy();
    const inner = tool({
      id: 'inner',
      description: 'Inner.',
      input: none,
      modules: [],
      run: () => ({ ok: true, did: 'x' }),
    });
    const nested = new NestedTools<Toy, Choice>(configure);
    const bindings = createToolBindings([], api, { nested });
    const refused = async (path: string[]) => {
      let message = '';
      await nested.run(
        runContext(
          inner.definition,
          {},
          {
            callId: 'c1',
            run: async () => {
              try {
                for await (const _ of bindings
                  .resolve(inner.definition)
                  ?.invoke(runContext(inner.definition, {}, { path })) ?? [])
                  void _;
              } catch (e) {
                message = (e as Error).message;
              }
              return done;
            },
          },
        ),
        { messages: [], tools: [inner], api, model: { name: 'm' }, limits: {} },
      );
      return message;
    };
    expect(await refused(['c1'])).toBe('');
    expect(await refused(['c2'])).toBe('inner belongs to a helper that is not running.');
  });

  it('a tool called without a turn has no nest', async () => {
    const outer = tool({
      id: 'outer',
      description: 'Uses a helper.',
      input: none,
      modules: [],
      run: ({ nest }) => ({ ok: true, did: nest ? 'nest' : 'none' }),
    });
    expect((await call(outer, toy(), {})).text).toBe('none');
  });
});
