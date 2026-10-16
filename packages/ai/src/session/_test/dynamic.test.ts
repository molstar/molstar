/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { callOf, fixture, tool, toyToolkit } from './helpers.js';

// A dynamic toolkit in a session (decision 25 of the prism-ui spec, `.v6/prism-ui/spec/README.md`; plan track A6):
// tools added between turns are offered at the next; a lazy skill disclosed during a turn brings its tools in that
// turn.

const none = { kind: 'object', fields: {} } as const;
const doubler = (id: string) =>
  tool({
    id,
    description: `Doubles the value (${id}).`,
    input: none,
    modules: ['state'],
    run: ({ api }) => {
      api.state.bump(api.state.value);
      return { ok: true, did: 'Doubled.', changed: ['value'] };
    },
  });

describe('a dynamic toolkit in a session', () => {
  it('offers a tool added between turns at the next turn', async () => {
    const toolkit = toyToolkit();
    const { session, api } = fixture(
      [{ text: 'Hello.' }, { calls: [callOf('double_value')] }, { text: 'Doubled it.' }],
      { toolkit },
    );
    api.state.bump(3);
    expect((await session.turn({ prompt: 'hi' }, {})).status).toBe('completed');
    toolkit.add({ tools: [doubler('double_value')], core: ['double_value'] });
    const r = await session.turn({ prompt: 'double it' }, {});
    expect(r.status).toBe('completed');
    expect(r.tools.map((t) => [t.name, t.status])).toEqual([['double_value', 'succeeded']]);
    expect(api.state.value).toBe(6);
  });

  it('loads a lazy skill disclosed during a turn, and runs its tool in that turn', async () => {
    const toolkit = toyToolkit();
    let removeStub = () => {};
    removeStub = toolkit.add({
      groups: [
        {
          id: 'math',
          title: 'Math',
          summary: 'arithmetic on the value',
          guidance: '',
          tools: [],
          pending: [{ id: 'triple_value' }],
          load: async () => {
            removeStub();
            toolkit.add({
              tools: [doubler('triple_value')],
              groups: [
                {
                  id: 'math',
                  title: 'Math',
                  summary: 'arithmetic on the value',
                  guidance: 'Use triple_value.',
                  tools: ['triple_value'],
                },
              ],
            });
          },
        },
      ],
    });
    const { session, api } = fixture(
      [
        { calls: [callOf('disclose_skills', { skills: ['math'] })] },
        { calls: [callOf('triple_value')] },
        { text: 'Done.' },
      ],
      { toolkit },
    );
    api.state.bump(2);
    const r = await session.turn({ prompt: 'triple it' }, {});
    expect(r.status).toBe('completed');
    expect(r.tools.map((t) => [t.name, t.status])).toEqual([
      ['disclose_skills', 'succeeded'],
      ['triple_value', 'succeeded'],
    ]);
    expect(api.state.value).toBe(4);
  });
});
