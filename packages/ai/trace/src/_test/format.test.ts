/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Message } from '@molstar/ai-agent';
import {
  cachedShare,
  credits,
  defaultRequestOf,
  defaultWindowPart,
  formatUsage,
  k,
  ms,
  tokensOf,
  turnStatusOf,
} from '../index.js';

const text = (t: string, role: 'user' | 'assistant' = 'user'): Message => ({
  role,
  parts: [{ kind: 'text', text: t }],
});

describe('formatting', () => {
  it('k: counts of tokens', () => {
    expect([0, 812, 999, 1000, 3210, 9999, 10_000, 41_300, 999_999, 1_000_000, 1_500_000, 12_000_000].map(k)).toEqual([
      '0',
      '812',
      '999',
      '1.0k',
      '3.2k',
      '10k',
      '10k',
      '41k',
      '1M',
      '1M',
      '1.5M',
      '12M',
    ]);
  });

  it('ms: durations', () => {
    expect([0, 340, 999.6, 1000, 1840, 9990, 10_000, 42_300].map(ms)).toEqual([
      '0 ms',
      '340 ms',
      '1.0 s',
      '1.0 s',
      '1.8 s',
      '10 s',
      '10 s',
      '42 s',
    ]);
  });

  it('credits: as people read them', () => {
    expect([0.04, 0.4, 1, 9.96, 41.2, 1200].map(credits)).toEqual([
      '0 credits',
      '0.4 credits',
      '1 credit',
      '10 credits',
      '41 credits',
      '1200 credits',
    ]);
  });

  it('formatUsage: tokens in (cached share), out, credits', () => {
    expect(formatUsage({ input: 3200, output: 22 })).toBe('3.2k in · 22 out');
    expect(formatUsage({ input: 3200, output: 22, cached: 2720 })).toBe('3.2k in (85% cached) · 22 out');
    expect(formatUsage({ input: 3200, output: 22, cached: 2720, written: 400, credits: 0.4 })).toBe(
      '3.2k in (85% cached, 400 written) · 22 out · 0.4 credits',
    );
    expect(formatUsage({ input: 100, output: 0, written: 100 })).toBe('100 in (100 written) · 0 out');
    expect(formatUsage({ input: 100, output: 5, estimated: true })).toBe('~100 in · ~5 out');
    expect(cachedShare({ input: 0, cached: 5 })).toBeUndefined();
    expect(cachedShare({ input: 10 })).toBeUndefined();
  });

  it('tokensOf: the input includes what the cache served and wrote', () => {
    const usage = { uncachedInput: 100, cacheRead: 800, cacheWrite: 50, output: 20, reasoning: 0 };
    expect(tokensOf(usage)).toEqual({ input: 950, output: 20, cached: 800, written: 50 });
    expect(tokensOf({ ...usage, reasoning: 7, estimated: true }, 2)).toEqual({
      input: 950,
      output: 20,
      cached: 800,
      written: 50,
      reasoning: 7,
      estimated: true,
      credits: 2,
    });
  });
});

describe('the default profile', () => {
  it('requestOf: the text after the request marker, and the Focus line before it', () => {
    expect(defaultRequestOf(text('Request: make it red'))).toEqual({ prompt: 'make it red' });
    expect(defaultRequestOf(text('Counter at 0.\n\nFocus: the display\n\nRequest: make it red'))).toEqual({
      prompt: 'make it red',
      focus: 'the display',
    });
    expect(
      defaultRequestOf(
        text('State.\n\nFocus: a\n\nNotes since your last turn:\n- x\n\nRequest: do\nit\n\nRequest: twice'),
      ),
    ).toEqual({ prompt: 'do\nit\n\nRequest: twice', focus: 'a' });
    // A Focus line in the prompt is the person’s words, not the host’s.
    expect(defaultRequestOf(text('Request: a\nFocus: b'))).toEqual({ prompt: 'a\nFocus: b' });
    expect(defaultRequestOf(text('just words'))).toEqual({ prompt: 'just words' });
    expect(defaultRequestOf(text(''))).toEqual({ prompt: '(no request recorded)' });
    expect(defaultRequestOf(text('Focus: \n\nRequest: x'))).toEqual({ prompt: 'x' });
  });

  it('windowPart: guidance, tools, state or conversation', () => {
    const definition = { tool: { id: 't' }, name: 't', description: 'd', inputSchema: {} };
    expect(defaultWindowPart({ role: 'system', parts: [{ kind: 'text', text: 'Be brief.' }] })).toBe('guidance');
    expect(
      defaultWindowPart({ role: 'system', parts: [{ kind: 'tool-definition', definition, deferred: false }] }),
    ).toBe('tools');
    expect(defaultWindowPart(text('Counter at 0.\n\nRequest: go'))).toBe('state');
    expect(defaultWindowPart(text('Request: go'))).toBe('conversation');
    expect(defaultWindowPart(text('plain'))).toBe('conversation');
    expect(defaultWindowPart(text('Counter at 0.\n\nRequest: go', 'assistant'))).toBe('conversation');
    expect(
      defaultWindowPart({
        role: 'user',
        parts: [{ kind: 'tool-result', requestId: 'r', tool: { id: 't' }, status: 'succeeded', parts: [] }],
      }),
    ).toBe('conversation');
  });

  it('turnStatusOf: the harness’s statuses as Prism’s', () => {
    const execution = (turn: object | undefined, phase: 'none' | 'live' = 'live') =>
      ({ runner: { phase }, ...(turn ? { turn } : {}) }) as never;
    const active = (pending: 'model' | 'tools', extra: object = {}) =>
      execution({ id: 't', pending: { kind: pending }, ...extra });
    expect(turnStatusOf('completed', execution(undefined), 't')).toBe('done');
    expect(turnStatusOf('failed', execution(undefined), 't')).toBe('failed');
    expect(turnStatusOf('stopped', execution(undefined), 't')).toBe('stopped');
    expect(turnStatusOf('reset', execution(undefined), 't')).toBe('stopped');
    expect(turnStatusOf('active', active('model'), 't')).toBe('thinking');
    expect(turnStatusOf('active', active('tools'), 't')).toBe('working');
    expect(turnStatusOf('active', active('tools', { park: { path: [] } }), 't')).toBe('waiting');
    expect(turnStatusOf('active', active('model', { blocked: 'stale-schema' }), 't')).toBe('waiting');
    expect(turnStatusOf('active', execution({ id: 't', pending: { kind: 'model' } }, 'none'), 't')).toBe('waiting');
    // A turn the Execution does not hold (its summary is not final yet).
    expect(turnStatusOf('active', execution(undefined), 't')).toBe('working');
    expect(turnStatusOf(undefined, active('model'), 'other')).toBe('working');
  });
});
