/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import type { Message, Provider } from '../../../messages/index.js';
import { type ModelStep, scriptedModel } from '../../../testing/index.js';
import type { HeadEdit } from '../../../harness/index.js';
import { createToolkit } from '../../../tools/index.js';
import { AgentSession } from '../../index.js';
import { callOf, guidance, modelsOf, tool, type Toy, toy } from '../../_test/helpers.js';
import { GUIDANCE_UPDATE, guidanceLabel, LABEL_CARRY_OVER, LABEL_MEMORY, LABEL_REQUEST } from '../../compose.js';
import {
  CARRY_OVER_HEADER,
  type CompactionDomain,
  checkedCompaction,
  checkWishes,
  compactionFor,
  DEFAULT_COMPACTION,
  type HeadEntry,
  parseCarryOver,
  planCompaction,
  SIZE_NOTE,
  toolkitDomain,
} from '../index.js';

const none = { kind: 'object', fields: {} } as const;
const quiet = (id: string, read = false) =>
  tool({
    id,
    description: `The ${id} tool.`,
    input: none,
    modules: [],
    read,
    run: () => ({ ok: true, did: `${id} done.` }),
  });

/** A toolkit with a core, two skills with a tool each, and two tools loaded by name. */
function kit() {
  return createToolkit<Toy>({
    tools: [
      quiet('get_state', true),
      quiet('set_look'),
      quiet('volume_load'),
      quiet('shape_add'),
      quiet('remove_scene'),
      quiet('move_scene'),
    ],
    core: ['get_state', 'set_look'],
    groups: [
      {
        id: 'volumes',
        title: 'Volumes',
        summary: 'maps',
        guidance: 'Skill volumes: load maps.',
        tools: ['volume_load'],
      },
      { id: 'shapes', title: 'Shapes', summary: 'shapes', guidance: 'Skill shapes: draw.', tools: ['shape_add'] },
    ],
  });
}

describe('the carry-over answer', () => {
  it('parses its sections and the retention block', () => {
    const c = parseCarryOver(`[wishes]
- "keep it for first-year students" (audience)
- "no jargon"
[/wishes]
[decisions]
Rejected a surface view of the pocket.
[/decisions]
[open]
The resistance scene.
[/open]
[notes]
none
[/notes]
[retention]
[tools]remove_scene, move_scene[/tools]
[skills]nothing[/skills]
[/retention]`);
    expect(c?.wishes).toEqual(['- "keep it for first-year students" (audience)', '- "no jargon"']);
    expect(c?.decisions).toBe('Rejected a surface view of the pocket.');
    expect(c?.retention).toEqual({ tools: ['remove_scene', 'move_scene'], skills: [] });
    expect(parseCarryOver('Sure, here is a summary.')).toBeUndefined();
  });

  it('keeps wishes that quote the person, word for word, and drops the rest', () => {
    const { kept, dropped } = checkWishes(
      ['- "Keep it for first-year students" (audience)', '- "always use blue"', '- be brief'],
      ['Tell the story of imatinib; keep it for first-year students, please.'],
    );
    expect(kept).toEqual(['- "Keep it for first-year students" (audience)']);
    expect(dropped).toEqual(['- "always use blue"', '- be brief']);
  });

  it('keeps only the person’s quotes of a wish that also quotes something else', () => {
    const { kept } = checkWishes(
      [
        '- "keep every description to one short sentence" The remember tool failed with "No memory here", so honor it anyway.',
      ],
      ['From now on keep every description to one short sentence.'],
    );
    expect(kept).toEqual(['- "keep every description to one short sentence"']);
  });

  it('reads accents, case and breaks of lines as the same words', () => {
    const { kept } = checkWishes(['- "Café au lait, please"'], ['I like cafe au\nlait, PLEASE.']);
    expect(kept).toHaveLength(1);
  });
});

describe('compaction settings', () => {
  it('refuse a target so close to the soft threshold that it would compact every turn', () => {
    expect(checkedCompaction(DEFAULT_COMPACTION)).toBe(DEFAULT_COMPACTION);
    expect(() =>
      checkedCompaction({
        ...DEFAULT_COMPACTION,
        triggerInputTokens: 10_000,
        softShare: 0.6,
        targetContextTokens: 32_000,
      }),
    ).toThrow(/too close to the soft threshold/);
  });

  it('use the window: compact at 80%, optimistically from 85% of that, down to at most 60% of it', () => {
    expect(compactionFor(1_000_000)).toMatchObject({
      triggerInputTokens: 800_000,
      softShare: 0.85,
      targetContextTokens: 408_000,
      recentTurnsTokens: 64_000,
    });
    expect(compactionFor(200_000)).toMatchObject({ triggerInputTokens: 160_000, recentTurnsTokens: 16_000 });
  });
});

// ------------------------------------------------------------------ the new head

/** The label a session gives the entry a test writes by its words (entries are read by label, never by text). */
function labelFor(role: 'system' | 'user' | 'assistant', body: string): Message['label'] | undefined {
  const at = body.lastIndexOf('Request: ');
  if (role === 'user' && at >= 0) return { kind: LABEL_REQUEST, data: { prompt: body.slice(at + 'Request: '.length) } };
  if (role === 'user' && body.startsWith(GUIDANCE_UPDATE))
    return guidanceLabel(body.slice(`${GUIDANCE_UPDATE}\n\n`.length), true);
  if (role === 'user' && body.startsWith(CARRY_OVER_HEADER)) return { kind: LABEL_CARRY_OVER };
  if (role === 'system' && body !== '') return guidanceLabel(body);
  return undefined;
}

/** Head entries for the plan: ids are the messages' labels. */
const text = (id: string, role: 'system' | 'user' | 'assistant', body: string): HeadEntry => {
  const label = labelFor(role, body);
  return { id, message: { role, parts: [{ kind: 'text', text: body }], ...(label ? { label } : {}) } };
};
const definition = (id: string, entryId: string, deferred = false): HeadEntry => ({
  id: entryId,
  message: {
    role: 'system',
    parts: [
      {
        kind: 'tool-definition',
        definition: { tool: { id }, name: id, description: id, inputSchema: { type: 'object' } },
        deferred,
      },
    ],
  },
});
const kept = (edit: HeadEdit) => edit.entries.filter((e): e is string => typeof e === 'string');
const added = (edit: HeadEdit) => edit.entries.filter((e): e is Message => typeof e !== 'string');
const textOf = (m: Message) => m.parts.map((p) => (p.kind === 'text' ? p.text : '')).join('');
const domain = toolkitDomain(kit());
const plain = { domain, pins: { skills: [], tools: [] }, settings: { recentTurnsTokens: 0, maxSummaryTokens: 4_000 } };

describe('the toolkit’s domain', () => {
  it('knows a request by its label, and its words by the label’s data', () => {
    const labelled = (kind: string, body: string, data?: Message['label']): HeadEntry => ({
      ...text('x', 'user', body),
      message: { role: 'user', parts: [{ kind: 'text', text: body }], label: data ?? { kind } },
    });
    const request = labelled('', 'State.\n\nRequest: go', { kind: LABEL_REQUEST, data: { prompt: 'Request: go' } });
    expect(domain.isRequest(request.message)).toBe(true);
    expect(domain.requestText?.(request.message)).toBe('Request: go');
    // A labelled entry of another kind is no request, whatever its text.
    expect(domain.isRequest(labelled(LABEL_CARRY_OVER, 'Request: go').message)).toBe(false);
  });

  it('knows the core, the groups and a request by its label', () => {
    expect(domain.coreToolIds).toEqual(['get_state', 'set_look', 'disclose_skills', 'disclose_tools']);
    expect(domain.groups.map((g) => g.id)).toEqual(['volumes', 'shapes']);
    expect(domain.isRequest(text('u', 'user', 'State.\n\nRequest: go').message)).toBe(true);
    expect(domain.isRequest(text('a', 'assistant', 'Request: go').message)).toBe(false);
    expect(domain.isRequest(text('u', 'user', 'Skill volumes: load maps.').message)).toBe(false);
    // Not by its text: without the label, a message with the marker is no request.
    const { label: _label, ...unlabelled } = text('u', 'user', 'State.\n\nRequest: go').message;
    expect(domain.isRequest(unlabelled)).toBe(false);
  });
});

describe('the new head', () => {
  it('keeps one system message, the latest guidance: older versions do not pile up', () => {
    const plan = planCompaction({
      entries: [
        text('g1', 'system', 'Guidance one.'),
        text('u1', 'user', 'Request: first'),
        text('a1', 'assistant', 'Done.'),
        text('g2', 'user', 'Updated guidance (it replaces the earlier version):\n\nGuidance two.'),
        text('u2', 'user', 'Request: second'),
      ],
      answer: undefined,
      ...plain,
    });
    const systems = added(plan.edit).filter((m) => m.role === 'system');
    expect(systems.map(textOf)).toEqual(['Guidance two.']);
    expect(kept(plan.edit).at(-1)).toBe('u2');
  });

  it('labelled: the latest guidance is the last update, given again as the system message, labelled', () => {
    const labelled = (id: string, role: 'system' | 'user', body: string, label: Message['label']): HeadEntry => {
      const entry = text(id, role, body);
      return { ...entry, message: { ...entry.message, label } as Message };
    };
    const plan = planCompaction({
      entries: [
        labelled('g1', 'system', 'Guidance one.', guidanceLabel('Guidance one.')),
        text('u1', 'user', 'Request: first'),
        text('a1', 'assistant', 'Done.'),
        labelled('g2', 'user', `${GUIDANCE_UPDATE}\n\nGuidance two.`, guidanceLabel('Guidance two.', true)),
        text('u2', 'user', 'Request: second'),
      ],
      answer: undefined,
      ...plain,
      memory: 'The memory.',
    });
    const systems = added(plan.edit).filter((m) => m.role === 'system');
    expect(systems.map(textOf)).toEqual(['Guidance two.']);
    expect(systems[0]?.label).toEqual(guidanceLabel('Guidance two.'));
    const notes = added(plan.edit).filter((m) => m.role === 'user');
    expect(notes.map((m) => m.label?.kind)).toEqual([LABEL_CARRY_OVER, LABEL_MEMORY]);
  });

  it('starts from the memory as it now stands (its events go with the old context)', () => {
    const plan = planCompaction({
      entries: [
        text('u1', 'user', 'The memory:\n## Notes\n- Old. [m1]\n\nRequest: first'),
        text('a1', 'assistant', 'memory: forget m1 (note): Old.'),
        text('u2', 'user', 'Request: second'),
      ],
      answer: undefined,
      ...plain,
      memory: 'The memory:\n## Notes\n- New. [m2]',
    });
    expect(added(plan.edit).some((m) => textOf(m).includes('New. [m2]'))).toBe(true);
    expect(kept(plan.edit)).not.toContain('u1');
  });

  it('may quote the host’s own sources as the person’s words', () => {
    const plan = planCompaction({
      entries: [
        text('u1', 'user', 'Request: first'),
        text('a1', 'assistant', 'Done.'),
        text('u2', 'user', 'Request: second'),
      ],
      answer: '[wishes]\n- "always label the chains"\n[/wishes]',
      ...plain,
      domain: { ...domain, wishSources: () => ['The person’s note: always label the chains.'] },
    });
    expect(plan.report.wishesKept).toBe(1);
  });
});

describe('skills in the new head', () => {
  const k = kit();
  const volumes = k.groups.find((x) => x.id === 'volumes');
  it('counts a skill loaded by its guidance, keeps it as guidance, and not one whose tools only were loaded', () => {
    if (!volumes) throw new Error('no skills');
    const plan = planCompaction({
      entries: [
        text('g', 'system', 'Guidance.'),
        text('u1', 'user', 'Request: first'),
        { id: 'skill-volumes', message: k.groupEntry(volumes) },
        definition('shape_add', 'def-shape_add'),
        text('a1', 'assistant', 'Done.'),
        text('u2', 'user', 'Request: second'),
      ],
      answer:
        '[decisions]\nnone\n[/decisions]\n[retention]\n[tools][/tools]\n[skills]volumes, shapes[/skills]\n[/retention]',
      domain,
      pins: { skills: [], tools: [] },
      settings: { recentTurnsTokens: 0, maxSummaryTokens: 8_000 },
    });
    expect(plan.report.skillsKept).toEqual(['volumes']);
    expect(plan.report.unknownRetention).toContain('shapes');
    expect(added(plan.edit)).toContainEqual(k.groupEntry(volumes));
  });

  it('keeps a pinned skill whatever the answer says', () => {
    if (!volumes) throw new Error('no skills');
    const plan = planCompaction({
      entries: [
        text('g', 'system', 'Guidance.'),
        text('u1', 'user', 'Request: first'),
        { id: 'skill-volumes', message: k.groupEntry(volumes) },
        text('u2', 'user', 'Request: second'),
      ],
      answer: undefined,
      ...plain,
      pins: { skills: ['volumes'], tools: [] },
    });
    expect(plan.report.skillsKept).toEqual(['volumes']);
  });
});

describe('declared tools in the new head', () => {
  const entries = [
    text('g', 'system', 'Guidance.'),
    definition('set_look', 'core'),
    definition('remove_scene', 'decl-remove', true),
    definition('move_scene', 'decl-move', true),
    text('u1', 'user', 'Request: first'),
    definition('remove_scene', 'load-remove'),
    text('a1', 'assistant', 'Loaded.'),
    text('u2', 'user', 'Request: second'),
  ];
  const tools = (pinned: string[]) => {
    const plan = planCompaction({ entries, answer: undefined, ...plain, pins: { skills: [], tools: pinned } });
    const parts = added(plan.edit).flatMap((m) => m.parts.flatMap((p) => (p.kind === 'tool-definition' ? [p] : [])));
    return {
      plan,
      loaded: parts.filter((p) => !p.deferred).map((p) => p.definition.tool.id),
      declared: parts.filter((p) => p.deferred).map((p) => p.definition.tool.id),
    };
  };

  it('keeps a retained tool loaded (up front now) and the others declared', () => {
    const retained = tools(['remove_scene']);
    expect(retained.loaded).toEqual(['set_look', 'remove_scene']);
    expect(retained.declared).toEqual(['move_scene']);
    expect(retained.plan.report.toolsKept).toEqual(['set_look', 'remove_scene']);
    const dropped = tools([]);
    expect(dropped.loaded).toEqual(['set_look']);
    expect(dropped.declared).toEqual(['remove_scene', 'move_scene']);
  });
});

// ------------------------------------------------------------------ in a session

/**
 * A model whose requests grow with the conversation: the main script answers the turns, and a compaction call (the
 * carry-over instruction last) is answered by `carryOver`. Stories' scripted model chose the answer by the request;
 * the conformance one plays a sequence, so two of them are told apart here by that instruction.
 */
function growingModel(script: ModelStep[], carryOver: string) {
  const main = scriptedModel([...script, ...Array.from({ length: 30 }, () => ({ text: 'Done.' }))], {
    capabilities: { contextWindowTokens: 10_000_000, maxOutputTokens: 16_384 },
  });
  const carry = scriptedModel(Array.from({ length: 20 }, () => ({ text: carryOver })));
  const requests: string[] = [];
  const asks = (messages: Message[]) =>
    messages.at(-1)?.parts.some((p) => p.kind === 'text' && p.text.startsWith('The conversation is being condensed'));
  const provider: Provider = {
    ...main.provider,
    invoke(input) {
      const messages = (input.request.payload as unknown as { messages: Message[] }).messages;
      const compaction = asks(messages);
      requests.push(`${compaction ? 'COMPACTION ' : ''}${JSON.stringify(messages)}`);
      return (compaction ? carry : main).provider.invoke(input);
    },
  };
  return { provider, requests, main, carry };
}

const pad = (label: string) => `${label} ${'lorem ipsum dolor sit amet '.repeat(400)}`;
const settingsFor = (s: { trigger: number; soft?: number; target: number; summary?: number }) => ({
  triggerInputTokens: s.trigger,
  softShare: s.soft ?? 0.85,
  targetContextTokens: s.target,
  maxSummaryTokens: s.summary ?? 4_000,
  instructionReserveTokens: 0,
  recentTurnsTokens: 0,
});

function sessionOver(model: ReturnType<typeof growingModel>, compaction: ReturnType<typeof settingsFor>, extra = {}) {
  return new AgentSession<Toy>({
    models: modelsOf(model.provider),
    toolkit: kit(),
    api: toy(),
    guidance,
    compaction,
    ...extra,
  });
}

describe('compaction in a session (scripted model, the real harness loop)', () => {
  it('condenses at a turn’s start: wishes quoted from the person survive, invented ones do not', async () => {
    const m = growingModel(
      [
        { text: 'Done.' },
        { calls: [callOf('disclose_tools', { tools: ['remove_scene'] })] },
        { text: 'Loaded.' },
        { text: 'Done.' },
        { text: 'Condensed and done.' },
      ],
      `[wishes]
- "keep it for first-year students" (audience)
- "always color the drug blue" (invented)
[/wishes]
[decisions]
The person preferred cartoons over surfaces.
[/decisions]
[open]
nothing
[/open]
[notes]
none
[/notes]
[retention]
[tools]remove_scene[/tools]
[skills][/skills]
[/retention]`,
    );
    const session = sessionOver(m, settingsFor({ trigger: 10_000, target: 5_000 }));
    const prompts = [
      `${pad('FIRST')} Keep it for first-year students.`,
      `${pad('SECOND')} Load what you need to remove scenes.`,
      `${pad('THIRD')} Anything else?`,
      'Now, the fourth request.',
    ];
    const results = [];
    for (const prompt of prompts) results.push(await session.turn({ prompt }));
    expect(results.map((r) => r.status)).toEqual(['completed', 'completed', 'completed', 'completed']);
    const reports = results.flatMap((r) => r.compactions ?? []);
    expect(reports.length).toBeGreaterThan(0);
    for (const report of reports) {
      expect(report).toMatchObject({ fallback: false, wishesKept: 1 });
      expect(report.wishesDropped).toEqual(['- "always color the drug blue" (invented)']);
    }
    expect(reports.at(-1)?.toolsKept).toContain('remove_scene');
    // The compaction call saw the old context; the request after it carries the notes, not the turns.
    expect(m.requests.some((r) => r.startsWith('COMPACTION') && r.includes('FIRST'))).toBe(true);
    const last = m.requests.at(-1) ?? '';
    expect(last).toContain(CARRY_OVER_HEADER);
    // The carried-over context is conversation, not a system message; the system messages are the guidance and the
    // tools.
    expect(last).toContain(`{"role":"user","parts":[{"kind":"text","text":"${CARRY_OVER_HEADER}`);
    expect(last.match(/"role":"system"/g)?.length).toBe(2);
    expect(last).toContain('keep it for first-year students');
    expect(last).not.toContain('always color the drug blue');
    expect(last).not.toContain('FIRST lorem');
    expect(last).toContain('Now, the fourth request.');
    expect(last).toContain('remove_scene');
  });

  it('accepts a compaction above its target, notes it for the next one, and does not compact again at once', async () => {
    const m = growingModel(
      [],
      '[wishes]\n[/wishes]\n[decisions]\nnone\n[/decisions]\n[open]\nnothing\n[/open]\n[notes]\nnone\n[/notes]',
    );
    // The target is below the fixed prefix (a long guidance and the core tools): no compaction can reach it.
    const session = sessionOver(m, settingsFor({ trigger: 14_000, target: 100 }), {
      guidance: `${guidance} ${'Be careful. '.repeat(3000)}`,
    });
    const results = [];
    for (const label of ['FIRST', 'SECOND', 'THIRD', 'FOURTH', 'FIFTH'])
      results.push(await session.turn({ prompt: pad(label) }));
    expect(results.every((r) => r.status === 'completed')).toBe(true);
    // One compaction (at the second turn's start, soft), accepted above its target. The later requests pass the
    // trigger again (~16k, 19k, 21k), but not the headroom past what the compaction left, so it does not compact again.
    const reports = results.flatMap((r) => r.compactions ?? []);
    expect(results[1]?.compactions).toHaveLength(1);
    expect(reports).toHaveLength(1);
    expect(reports[0]).toMatchObject({ overTarget: true, targetOutOfReach: true });
    // What cannot be condensed is over the target by itself: no note asks to condense harder.
    expect(m.requests.at(-1)).not.toContain(SIZE_NOTE);
    expect(m.requests.filter((r) => r.startsWith('COMPACTION'))).toHaveLength(1);
  });

  it('asks the next compaction to condense harder when the carried-over part is what is over the target', async () => {
    const m = growingModel(
      [],
      `[wishes]\n[/wishes]\n[decisions]\nnone\n[/decisions]\n[open]\nnothing\n[/open]\n[notes]\n${'A long note. '.repeat(1500)}\n[/notes]`,
    );
    const session = sessionOver(m, settingsFor({ trigger: 6_500, target: 3_000, summary: 8_000 }));
    const results = [];
    for (const label of ['FIRST', 'SECOND', 'THIRD', 'FOURTH', 'FIFTH'])
      results.push(await session.turn({ prompt: pad(label) }));
    expect(results.every((r) => r.status === 'completed')).toBe(true);
    const report = results.flatMap((r) => r.compactions ?? [])[0];
    expect(report).toMatchObject({ overTarget: true });
    expect(report?.targetOutOfReach).toBeUndefined();
    expect(report?.fixedTokens).toBeLessThan(3_000);
    const after = m.requests.filter((r) => !r.startsWith('COMPACTION')).at(-1) ?? '';
    expect(after).toContain(SIZE_NOTE);
    expect(after).toMatch(/what is carried over should be at most about [\d,]+/);
  });

  it('condenses mid-turn: the request and the last step stay, the turn so far is carried, the turn finishes', async () => {
    const reads = Array.from({ length: 6 }, (_, look) => ({ calls: [callOf('get_state')], text: pad(`READ${look}`) }));
    const m = growingModel(
      [...reads, { text: 'Read them all.' }],
      `[wishes]
[/wishes]
[this turn]
Reading the state again and again; five reads done, the last next.
[/this turn]
[retention]
[tools][/tools]
[skills][/skills]
[/retention]`,
    );
    const session = sessionOver(m, settingsFor({ trigger: 9_000, soft: 10, target: 4_000 }));
    const r = await session.turn({ prompt: 'Read the state, one by one.' });
    expect(r.status).toBe('completed');
    expect(r.text).toBe('Read them all.');
    expect(r.compactions?.length).toBeGreaterThan(0);
    const after = m.requests.filter((x) => !x.startsWith('COMPACTION')).at(-1) ?? '';
    expect(after).toContain('This turn so far:');
    expect(after).toContain('Read the state, one by one.');
  });

  it('carries a digest when the compaction call gives no usable answer, keeping earlier wishes', async () => {
    const m = growingModel(
      [{ text: 'One.' }, { text: 'Two.' }, { text: 'Three.' }, { text: 'Four.' }],
      'Sorry, I cannot do that.',
    );
    const session = sessionOver(m, settingsFor({ trigger: 8_000, soft: 0.7, target: 3_000, summary: 8_000 }));
    await session.turn({ prompt: `${pad('FIRST')} Make scene two brighter.` });
    await session.turn({ prompt: `${pad('SECOND')} And the third?` });
    const r = await session.turn({ prompt: 'Third request.' });
    const after = await session.turn({ prompt: 'Fourth request.' });
    expect([r.status, after.status]).toEqual(['completed', 'completed']);
    const reports = [...(r.compactions ?? []), ...(after.compactions ?? [])];
    expect(reports[0]).toMatchObject({ fallback: true });
    const withDigest =
      m.requests.filter((x) => !x.startsWith('COMPACTION') && x.includes('Earlier turns:')).at(-1) ?? '';
    expect(withDigest).toContain('Asked: FIRST lorem');
    expect(withDigest).toContain('Replied: One.');
  });

  it('follows each turn’s model: its window sets the turn’s budget', async () => {
    const m = growingModel([{ calls: [callOf('get_state')] }, { text: 'Done.' }], '');
    // The provider reports a 1M window for the model "big" and 10M for the others.
    const provider: Provider = {
      ...m.provider,
      capabilities: async (c) => ({
        ...(await m.provider.capabilities(c)),
        contextWindowTokens: c.target.model === 'big' ? 1_000_000 : 10_000_000,
      }),
    };
    const session = new AgentSession<Toy>({
      models: modelsOf(provider),
      toolkit: kit(),
      api: toy(),
      guidance,
    });
    const client = await session.client();
    // While the turn runs, the limit on its uncached input: 4 × the trigger (80% of a 1M window).
    let limit: number | undefined;
    await session.turn(
      { prompt: 'hello', model: { provider: 'scripted', model: 'big' } },
      {
        onTool: (e) => {
          if (e.status === 'started')
            void client.execution().then((x) => {
              limit = x.turn?.budget.limits.uncachedInput;
            });
        },
      },
    );
    expect(limit).toBe(3_200_000);
    expect(session.budget({ provider: 'scripted', model: 'big' }).initial.uncachedInput).toBe(3_200_000);
  });
});

// A host's domain: it pins a skill and gives the memory.
describe('the host’s part of the domain', () => {
  it('gives pins and the memory to a compaction in a session', async () => {
    const m = growingModel(
      [{ text: 'One.' }, { text: 'Two.' }, { text: 'Three.' }, { text: 'Four.' }],
      '[decisions]\nnone\n[/decisions]\n[retention]\n[tools][/tools]\n[skills][/skills]\n[/retention]',
    );
    const memory: CompactionDomain['memory'] = () => 'The memory:\n- Keeps it simple. [m1]';
    const session = sessionOver(m, settingsFor({ trigger: 9_000, target: 4_000, summary: 8_000 }), {
      compactionDomain: { memory },
    });
    const results = [];
    for (const label of ['FIRST', 'SECOND', 'THIRD', 'FOURTH'])
      results.push(await session.turn({ prompt: pad(label) }));
    expect(results.flatMap((r) => r.compactions ?? []).length).toBeGreaterThan(0);
    expect(m.requests.at(-1)).toContain('Keeps it simple. [m1]');
  });
});
