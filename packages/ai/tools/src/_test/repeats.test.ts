/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { createToolBindings, RepeatedReads } from '../index.js';
import { call, runContext, toy, tool } from './helpers.js';

const none = { kind: 'object', fields: { q: { kind: 'string', optional: true } } } as const;
const get = tool({
  id: 'get',
  description: 'Reads.',
  input: none,
  modules: ['state'],
  read: true,
  run: ({ api }) => ({ ok: true, did: `Value ${api.state.value}.` }),
});
const bump = tool({
  id: 'bump',
  description: 'Writes.',
  input: none,
  modules: ['state'],
  run: ({ api }) => {
    api.state.bump();
    return { ok: true, did: 'Bumped.' };
  },
});

describe('repeated reads (tools-and-session.md §1)', () => {
  it('the same read with the same arguments while the version is unchanged answers in one line', async () => {
    const api = toy();
    const repeats = new RepeatedReads({ version: () => api.state.version });
    expect((await call(get, api, { q: 'a' }, { repeats })).text).toBe('Value 0.');
    expect((await call(get, api, { q: 'a' }, { repeats })).text).toBe(
      'Unchanged: nothing changed since your earlier get call with these arguments, so its result above is current. Do not read it again; act on it.',
    );
    // Other arguments are another read.
    expect((await call(get, api, { q: 'b' }, { repeats })).text).toBe('Value 0.');
  });

  it('a changed version reads fresh', async () => {
    const api = toy();
    const repeats = new RepeatedReads({ version: () => api.state.version });
    await call(get, api, {}, { repeats });
    await call(bump, api, {}, { repeats });
    expect((await call(get, api, {}, { repeats })).text).toBe('Value 1.');
  });

  it('a version that cannot be read is no reason to repeat', async () => {
    const api = toy();
    const repeats = new RepeatedReads({
      version: () => {
        throw new Error('no state');
      },
    });
    await call(get, api, {}, { repeats });
    expect((await call(get, api, {}, { repeats })).text).toBe('Value 0.');
  });

  it('a version of undefined is still a version', async () => {
    const api = toy();
    const repeats = new RepeatedReads({ version: () => undefined });
    expect((await call(get, api, {}, { repeats })).text).toBe('Value 0.');
    expect((await call(get, api, {}, { repeats })).text).toContain('Unchanged');
  });

  it('after the streak of reads in a row, each result says to act; an action starts it again', async () => {
    const api = toy();
    const repeats = new RepeatedReads({ version: () => api.state.version, streak: 3 });
    const texts: string[] = [];
    for (const q of ['a', 'b', 'c', 'd']) texts.push((await call(get, api, { q }, { repeats })).text);
    expect(texts.map((t) => t.includes('reads in a row'))).toEqual([false, false, true, true]);
    expect(texts[2]).toBe(
      'Value 0.\n\n(3 reads in a row without a change: you have what you need. Act now with what you know; read more only for one specific missing fact.)',
    );
    await call(bump, api, {}, { repeats });
    expect((await call(get, api, { q: 'e' }, { repeats })).text).toBe('Value 1.');
  });

  it('the default streak is six', async () => {
    const api = toy();
    const repeats = new RepeatedReads({ version: () => api.state.version });
    const noted: boolean[] = [];
    for (let i = 0; i < 6; i++)
      noted.push((await call(get, api, { q: String(i) }, { repeats })).text.includes('in a row'));
    expect(noted).toEqual([false, false, false, false, false, true]);
  });

  it('reset starts a turn over', async () => {
    const api = toy();
    const repeats = new RepeatedReads({ version: () => api.state.version, streak: 2 });
    await call(get, api, {}, { repeats });
    repeats.reset();
    expect((await call(get, api, {}, { repeats })).text).toBe('Value 0.');
  });

  it('the exempt tools count as neither a read nor an action (by default the disclosure tools)', async () => {
    const api = toy();
    const disclose = tool({
      id: 'disclose_skills',
      description: 'D.',
      input: none,
      modules: [],
      run: () => ({ ok: true, did: 'Loaded.' }),
    });
    const other = tool({
      id: 'other',
      description: 'O.',
      input: none,
      modules: [],
      run: () => ({ ok: true, did: 'Done.' }),
    });
    const repeats = new RepeatedReads({ version: () => api.state.version, streak: 3 });
    await call(get, api, { q: 'a' }, { repeats });
    await call(get, api, { q: 'b' }, { repeats });
    await call(disclose, api, {}, { repeats });
    expect((await call(get, api, { q: 'c' }, { repeats })).text).toContain('3 reads in a row');
    const custom = new RepeatedReads({ version: () => 0, streak: 3, exempt: ['other'] });
    await call(get, api, { q: 'a' }, { repeats: custom });
    await call(get, api, { q: 'b' }, { repeats: custom });
    await call(other, api, {}, { repeats: custom });
    expect((await call(get, api, { q: 'c' }, { repeats: custom })).text).toContain('3 reads in a row');
    await call(disclose, api, {}, { repeats: custom });
    expect((await call(get, api, { q: 'd' }, { repeats: custom })).text).toBe('Value 0.');
  });

  it('is bound with the tools', async () => {
    const api = toy();
    const bindings = createToolBindings([get], api, {
      repeats: new RepeatedReads({ version: () => api.state.version }),
    });
    const handler = bindings.resolve(get.definition);
    const run = async () => {
      let text = '';
      for await (const f of handler?.invoke(runContext(get.definition, {})) ?? [])
        if (f.kind === 'result') text = f.parts.map((p) => (p.kind === 'text' ? p.text : '')).join('');
      return text;
    };
    expect(await run()).toBe('Value 0.');
    expect(await run()).toContain('Unchanged');
    expect(bindings.resolve({ ...get.definition, tool: { id: 'unknown' } })).toBeUndefined();
  });
});
