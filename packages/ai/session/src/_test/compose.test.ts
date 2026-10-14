/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import type { ToolDefinition } from '@molstar/ai-agent';
import { fingerprint, LABEL_DEFINITIONS } from '@molstar/ai-tools';
import {
  FOCUS_MARKER,
  GUIDANCE_UPDATE,
  LABEL_GUIDANCE,
  LABEL_GUIDANCE_UPDATE,
  LABEL_REQUEST,
  NOTES_HEADING,
  REQUEST_MARKER,
  requestOf,
  turnInput,
} from '../index.js';
import { textOf } from './helpers.js';

const def = (id: string): ToolDefinition => ({
  tool: { id },
  name: id,
  description: id,
  inputSchema: { type: 'object' },
});

describe('the request text (tools-and-session.md §4)', () => {
  it('is [preamble][state], the Focus, the notes, the request: a part absent leaves out its separator', () => {
    expect(
      requestOf({
        prompt: 'hello',
        state: 'State.',
        focus: 'the display',
        notes: ['a', 'b'],
        preamble: 'Memory.\n\n',
      }),
    ).toBe('Memory.\n\nState.\n\nFocus: the display\n\nNotes since your last turn:\n- a\n- b\n\nRequest: hello');
    expect(requestOf({ prompt: 'hello' })).toBe('Request: hello');
    expect(requestOf({ prompt: 'hello', state: 'State.' })).toBe('State.\n\nRequest: hello');
    // The Focus line is there when a focus is given, even an empty one ("none").
    expect(requestOf({ prompt: 'hello', focus: 'none' })).toBe('Focus: none\n\nRequest: hello');
    expect(requestOf({ prompt: 'hello', notes: [] })).toBe('Request: hello');
    expect(requestOf({ prompt: 'hello', notes: ['a'], notesHeading: 'The person changed:' })).toBe(
      'The person changed:\n- a\n\nRequest: hello',
    );
    expect(NOTES_HEADING).toBe('Notes since your last turn:');
    expect([REQUEST_MARKER, FOCUS_MARKER]).toEqual(['Request: ', 'Focus: ']);
  });
});

describe('the messages of a turn', () => {
  const base = { prompt: 'go', guidance: 'Guidance.', tools: [def('a')], deferred: [def('b')] };

  it('first: the guidance, one message of definitions (loaded, then declared), the request', () => {
    const [guidance, tools, request, ...rest] = turnInput({ ...base, first: true });
    expect(rest).toEqual([]);
    expect(guidance).toEqual({
      role: 'system',
      parts: [{ kind: 'text', text: 'Guidance.' }],
      label: { kind: LABEL_GUIDANCE, key: fingerprint('Guidance.') },
    });
    expect(tools?.label).toEqual({ kind: LABEL_DEFINITIONS });
    expect(tools?.parts.map((p) => (p.kind === 'tool-definition' ? [p.definition.tool.id, p.deferred] : []))).toEqual([
      ['a', false],
      ['b', true],
    ]);
    expect(textOf(request as never)).toBe('Request: go');
    expect(request?.label).toEqual({ kind: LABEL_REQUEST, data: { prompt: 'go' } });
  });

  it('labels the request with the person’s words and focus, whatever the text around them', () => {
    const [request] = turnInput({ ...base, first: false, prompt: 'Request: twice', focus: 'A:45', state: 'S' });
    expect(request?.label).toEqual({ kind: LABEL_REQUEST, data: { prompt: 'Request: twice', focus: 'A:45' } });
  });

  it('later: the request alone; or what changed, as definitions (system) and guidance (user), then the request', () => {
    expect(turnInput({ ...base, first: false })).toHaveLength(1);
    const messages = turnInput({
      ...base,
      first: false,
      catchUp: { tools: [def('a')], deferred: [def('c')], guidance: true },
    });
    expect(messages.map((m) => m.role)).toEqual(['system', 'user', 'user']);
    expect(textOf(messages[1] as never)).toBe(`${GUIDANCE_UPDATE}\n\nGuidance.`);
    expect(messages.map((m) => m.label?.kind)).toEqual([LABEL_DEFINITIONS, LABEL_GUIDANCE_UPDATE, LABEL_REQUEST]);
    expect(turnInput({ ...base, first: false, catchUp: { tools: [], guidance: true } }).map((m) => m.role)).toEqual([
      'user',
      'user',
    ]);
  });

  it('carries attachments after the text', () => {
    const reference = { id: 'f', mediaType: 'application/pdf' } as never;
    const [message] = turnInput({ ...base, first: false, attachments: [reference] });
    expect(message?.parts.map((p) => p.kind)).toEqual(['text', 'reference']);
  });
});
