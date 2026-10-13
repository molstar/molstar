/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { RepeatedReads, type ToolResult } from '../index.js';
import { call, toy, tool } from './helpers.js';

const none = { kind: 'object', fields: {} } as const;
const make = (id: string, read: boolean, run: () => ToolResult | Promise<ToolResult>, extra = {}) =>
  tool({ id, description: `${id}.`, input: none, modules: [], read, run, ...extra });

describe('the effect on the result frame (tools-and-session.md §1)', () => {
  const cases: [string, boolean, () => ToolResult, string][] = [
    ['success of a read', true, () => ({ ok: true, did: 'Read.' }), 'none'],
    ['success of a write', false, () => ({ ok: true, did: 'Wrote.' }), 'applied'],
    ['the tool’s own failure, of a read', true, () => ({ ok: false, code: 'x', message: 'No.' }), 'none'],
    ['the tool’s own failure, of a write', false, () => ({ ok: false, code: 'x', message: 'No.' }), 'none'],
    ['a result that says it', false, () => ({ ok: true, did: 'Nothing to do.', effect: 'none' }), 'none'],
    [
      'a failure that says it (the change was made, the rest was not)',
      false,
      () => ({ ok: false, code: 'x', message: 'Half.', effect: 'applied' }),
      'applied',
    ],
    ['a read that says it changed something', true, () => ({ ok: true, did: 'Hm.', effect: 'applied' }), 'applied'],
  ];
  for (const [name, read, run, effect] of cases)
    it(`${name}: ${effect}`, async () => {
      expect((await call(make('t', read, run), toy(), {})).result.effect).toBe(effect);
    });

  it('a throw: none for a read, unknown for a write', async () => {
    const boom = () => {
      throw new Error('Boom.');
    };
    expect((await call(make('r', true, boom), toy(), {})).result).toMatchObject({ status: 'failed', effect: 'none' });
    expect((await call(make('w', false, boom), toy(), {})).result).toMatchObject({
      status: 'failed',
      effect: 'unknown',
    });
  });

  it('a failed argument check is a failure before any change: none', async () => {
    const t = tool({
      id: 'needs',
      description: 'Needs a name.',
      input: { kind: 'object', fields: { name: { kind: 'string' } } },
      modules: [],
      run: () => ({ ok: true, did: 'Named.' }),
    });
    expect((await call(t, toy(), {})).result).toMatchObject({ status: 'failed', effect: 'none' });
  });

  it('a repeated read answered in one line is a read: none', async () => {
    const api = toy();
    const repeats = new RepeatedReads({ version: () => api.state.version });
    const t = make('get', true, () => ({ ok: true, did: 'Value.' }));
    await call(t, api, {}, { repeats });
    expect((await call(t, api, {}, { repeats })).result).toMatchObject({ status: 'succeeded', effect: 'none' });
  });
});

describe('afterWrite (tools-and-session.md §1)', () => {
  it('runs after a write, before the result, with the tool and the call', async () => {
    const order: string[] = [];
    const afterWrite = async (c: { toolId: string; callId: string }) => {
      order.push(`saved ${c.toolId} ${c.callId}`);
    };
    const t = make('w', false, () => {
      order.push('ran');
      return { ok: true, did: 'Wrote.' };
    });
    const r = await call(t, toy(), {}, { afterWrite }, { callId: 'c-9' });
    expect(order).toEqual(['ran', 'saved w c-9']);
    expect(r.result).toMatchObject({ status: 'succeeded', effect: 'applied' });
  });

  it('does not run after a read, a failure, or a write that says it changed nothing', async () => {
    const saved: string[] = [];
    const afterWrite = async (c: { toolId: string }) => void saved.push(c.toolId);
    await call(
      make('r', true, () => ({ ok: true, did: 'Read.' })),
      toy(),
      {},
      { afterWrite },
    );
    await call(
      make('f', false, () => ({ ok: false, code: 'x', message: 'No.' })),
      toy(),
      {},
      { afterWrite },
    );
    await call(
      make('z', false, () => ({ ok: true, did: 'Same.', effect: 'none' })),
      toy(),
      {},
      { afterWrite },
    );
    await call(
      make('t', false, () => {
        throw new Error('Boom.');
      }),
      toy(),
      {},
      { afterWrite },
    );
    expect(saved).toEqual([]);
  });

  it('when it throws, the call fails transient with effect applied', async () => {
    const afterWrite = async () => {
      throw new Error('disk full');
    };
    const r = await call(
      make('w', false, () => ({ ok: true, did: 'Wrote.' })),
      toy(),
      {},
      { afterWrite },
    );
    expect(r.result).toMatchObject({ status: 'failed', transient: true, effect: 'applied' });
    expect(r.text).toBe('The change could not be saved: disk full');
  });
});

describe('the call’s identity and the registration’s flags', () => {
  it('run gets the call id, the key and whether an earlier attempt may have run', async () => {
    const seen: unknown[] = [];
    const t = tool({
      id: 'who',
      description: 'Who.',
      input: none,
      modules: [],
      run: ({ call }) => {
        seen.push(call);
        return { ok: true, did: 'Done.' };
      },
    });
    await call(t, toy(), {}, {}, { callId: 'c1', key: 'k1', prior: 'none' });
    await call(t, toy(), {}, {}, { callId: 'c1', key: 'k1', prior: 'possible' });
    expect(seen).toEqual([
      { callId: 'c1', key: 'k1', prior: 'none' },
      { callId: 'c1', key: 'k1', prior: 'possible' },
    ]);
  });

  it('carries idempotent, read, approval and optional on the registration', () => {
    const t = tool({
      id: 'flags',
      description: 'Flags.',
      input: none,
      modules: ['state', 'extra'],
      optional: ['extra'],
      idempotent: true,
      approval: true,
      approvedBy: (_args, api) => api.state.value > 0,
      run: () => ({ ok: true, did: 'Done.' }),
    });
    expect(t).toMatchObject({ idempotent: true, approval: true, modules: ['state', 'extra'], optional: ['extra'] });
    expect(t.read).toBeUndefined();
    const api = toy();
    expect(t.approvedBy?.({}, api)).toBe(false);
    api.state.bump();
    expect(t.approvedBy?.({}, api)).toBe(true);
    expect(make('r', true, () => ({ ok: true, did: 'x' })).read).toBe(true);
  });

  it('a module the host lacks is undefined for a tool that does without it', async () => {
    let extra: unknown = 'unset';
    const t = tool({
      id: 'optional',
      description: 'Optional.',
      input: none,
      modules: ['extra'],
      optional: ['extra'],
      read: true,
      run: ({ api }) => {
        extra = api.extra;
        return { ok: true, did: 'Done.' };
      },
    });
    await call(t, toy(), {});
    expect(extra).toBeUndefined();
  });
});
