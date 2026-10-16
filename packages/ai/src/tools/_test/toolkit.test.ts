/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import type { Message } from '../../messages/index.js';
import {
  closestHint,
  createToolkit,
  definitionEntries,
  fingerprint,
  LABEL_DEFINITIONS,
  LABEL_ENTRY,
  LABEL_SKILL,
  type ToolGroup,
  type ToolkitSpec,
  type ToolRegistration,
} from '../index.js';
import { call, definedIn, type Toy, textsIn, toy, tool } from './helpers.js';

const none = { kind: 'object', fields: {} } as const;
const make = (id: string, description = `${id}. More about ${id}.`): ToolRegistration<Toy> =>
  tool({ id, description, input: none, modules: [], run: () => ({ ok: true, did: `${id} done.` }) });

const VOLUMES: ToolGroup<Toy> = {
  id: 'volumes',
  title: 'Density maps',
  summary: 'maps drawn with the structure',
  tags: ['maps'],
  guidance: 'Skill volumes: maps. Use load_volume, then set_volume_look; remove_volume_look takes one away.',
  tools: ['load_volume', 'set_volume_look'],
  extras: ['remove_volume_look'],
};
const SELECTIONS: ToolGroup<Toy> = {
  id: 'selections',
  title: 'Selections',
  summary: 'select by distance and logic',
  guidance: 'Skill selections: query_reference gives symbols.',
  tools: ['query_reference'],
  modules: ['extra'],
  unavailable: (api) =>
    api.extra ? undefined : { message: 'Selections cannot be evaluated here (no analysis).', hint: 'Use targets.' },
};
const PDBE = { id: 'pdbe', summary: 'PDBe data notes', text: 'Source pdbe: structures and maps.', tools: ['invoke'] };
const NOTES = { id: 'notes', summary: 'house notes', text: 'Notes: be brief.' };

const ids = [
  'get_state',
  'set_value',
  'remove_scene',
  'move_scene',
  'load_volume',
  'set_volume_look',
  'remove_volume_look',
  'query_reference',
  'invoke',
];
const spec = (over: Partial<ToolkitSpec<Toy>> = {}): ToolkitSpec<Toy> => ({
  tools: ids.map((id) => make(id)),
  core: ['get_state', 'set_value'],
  groups: [VOLUMES, SELECTIONS],
  entries: [PDBE, NOTES],
  ...over,
});
const kit = createToolkit<Toy>(spec());

/** Calls a disclosure tool against a conversation that holds these tools and keys. */
async function disclose(name: string, args: unknown, held: { loaded?: string[]; keys?: string[] } = {}, api = toy()) {
  api.loaded.push(...(held.loaded ?? []));
  api.keys.push(...(held.keys ?? []));
  const t = kit.disclosureTools().find((d) => d.definition.tool.id === name);
  if (!t) throw new Error(name);
  return call(t, api, args as never);
}

describe('a toolkit: one place for each tool (tools-and-session.md §3)', () => {
  it('places the core, the groups’ tools and extras, an entry’s tools and every other tool singly', () => {
    expect(kit.placementOf('get_state')).toEqual({ kind: 'core' });
    expect(kit.placementOf('disclose_skills')).toEqual({ kind: 'core' });
    expect(kit.placementOf('disclose_tools')).toEqual({ kind: 'core' });
    expect(kit.placementOf('load_volume')).toMatchObject({ kind: 'group', group: { id: 'volumes' }, extra: false });
    expect(kit.placementOf('remove_volume_look')).toMatchObject({
      kind: 'group',
      group: { id: 'volumes' },
      extra: true,
    });
    expect(kit.placementOf('invoke')).toMatchObject({ kind: 'entry', entry: { id: 'pdbe' } });
    expect(kit.placementOf('remove_scene')).toEqual({ kind: 'single' });
    expect(kit.placementOf('nothing')).toEqual({ kind: 'single' });
  });

  it('the core is what is offered from the first turn, the disclosure tools with it', () => {
    expect(kit.coreDefinitions().map((d) => d.tool.id)).toEqual([
      'get_state',
      'set_value',
      'disclose_skills',
      'disclose_tools',
    ]);
    expect(kit.allDefinitions().map((d) => d.tool.id)).toEqual([...ids, 'disclose_skills', 'disclose_tools']);
    expect(kit.tools.map((t) => t.definition.tool.id)).toEqual(kit.allDefinitions().map((d) => d.tool.id));
    expect(kit.disclosureTools().map((t) => t.definition.tool.id)).toEqual(['disclose_skills', 'disclose_tools']);
  });

  it('a core that lists the disclosure tools is the same core', () => {
    const k = createToolkit<Toy>(spec({ core: ['get_state', 'disclose_skills', 'set_value', 'disclose_tools'] }));
    expect(k.coreDefinitions().map((d) => d.tool.id)).toEqual([
      'get_state',
      'disclose_skills',
      'set_value',
      'disclose_tools',
    ]);
  });

  it('throws, naming the tool, when it has two places', () => {
    expect(() => createToolkit<Toy>(spec({ core: ['get_state', 'load_volume'] }))).toThrow(
      'Tool "load_volume" has more than one place: the core and group "volumes".',
    );
    expect(() => createToolkit<Toy>(spec({ groups: [VOLUMES, { ...SELECTIONS, extras: ['load_volume'] }] }))).toThrow(
      'Tool "load_volume" has more than one place: group "volumes" and group "selections" (extras).',
    );
    expect(() => createToolkit<Toy>(spec({ entries: [{ ...PDBE, tools: ['remove_volume_look'] }] }))).toThrow(
      'Tool "remove_volume_look" has more than one place: group "volumes" (extras) and entry "pdbe".',
    );
    expect(() => createToolkit<Toy>(spec({ groups: [{ ...VOLUMES, extras: ['load_volume'] }] }))).toThrow(
      'Tool "load_volume" has more than one place',
    );
  });

  it('throws, naming the tool, when a declaration names one that is not registered, or one is registered twice', () => {
    expect(() => createToolkit<Toy>(spec({ core: ['get_state', 'ghost'] }))).toThrow(
      'Tool "ghost" is listed in the core but is not registered.',
    );
    expect(() => createToolkit<Toy>(spec({ groups: [{ ...VOLUMES, tools: ['load_volume', 'ghost'] }] }))).toThrow(
      'Tool "ghost" is listed in group "volumes" but is not registered.',
    );
    expect(() => createToolkit<Toy>(spec({ tools: [...ids, 'get_state'].map((id) => make(id)) }))).toThrow(
      'Tool "get_state" is registered twice.',
    );
    expect(() => createToolkit<Toy>(spec({ tools: [...ids, 'disclose_skills'].map((id) => make(id)) }))).toThrow(
      'Tool "disclose_skills" is a disclosure tool',
    );
    expect(() => createToolkit<Toy>(spec({ groups: [VOLUMES, VOLUMES] }))).toThrow('Group "volumes" is declared twice');
    expect(() => createToolkit<Toy>(spec({ entries: [{ ...NOTES, id: 'volumes' }] }))).toThrow(
      '"volumes" is declared twice',
    );
  });
});

describe('the disclosure tools’ catalog is generated from the toolkit', () => {
  const [skills, tools] = kit.disclosureTools();

  it('disclose_skills lists the skills and the other texts', () => {
    expect(skills?.definition.description).toBe(
      'Load skills: their guidance and their tools, all a request needs in one call, once per conversation. Skills: "volumes" (maps drawn with the structure); "selections" (select by distance and logic). Other texts, loaded the same way: "pdbe" (PDBe data notes); "notes" (house notes).',
    );
    expect(skills?.definition.inputSchema).toMatchObject({
      properties: {
        skills: {
          type: 'array',
          items: { type: 'string', enum: ['volumes', 'selections'] },
          description: 'Skill names.',
        },
      },
      required: ['skills'],
    });
  });

  it('disclose_tools lists the tools that load by name, by their first sentence', () => {
    expect(tools?.definition.description).toBe(
      "Load tools by name, all you need in one call; they stay available for the rest of the conversation. Tools: remove_scene (remove_scene.); move_scene (move_scene.). Tools that belong to a skill come with disclose_skills; a skill's guidance names the ones it loads by name.",
    );
    expect(tools?.definition.inputSchema).toMatchObject({
      properties: { tools: { items: { enum: ['remove_scene', 'move_scene'] } } },
    });
  });

  it('their ids can be named, and the texts that name them follow', () => {
    const k = createToolkit<Toy>(spec({ names: { skills: 'load_skill', tools: 'load_tool' } }));
    expect(k.disclosureTools().map((t) => t.definition.tool.id)).toEqual(['load_skill', 'load_tool']);
    expect(k.disclosureTools()[1]?.definition.description).toContain('come with load_skill;');
    expect(k.groupText(VOLUMES)).toContain('By name when needed (load_tool, or call it)');
    expect(k.placementOf('load_skill')).toEqual({ kind: 'core' });
  });

  it('no catalog, no enum when there is nothing to list', () => {
    const k = createToolkit<Toy>({ tools: [make('a')], core: ['a'], groups: [] });
    const [s, t] = k.disclosureTools();
    expect(s?.definition.inputSchema).toMatchObject({ properties: { skills: { items: { type: 'string' } } } });
    expect(JSON.stringify(s?.definition.inputSchema)).not.toContain('enum');
    expect(JSON.stringify(t?.definition.inputSchema)).not.toContain('enum');
    expect(t?.definition.description).not.toContain('Tools:');
  });

  it('their JSON Schemas are exactly these, byte for byte (what the model sees does not change with the input)', () => {
    const k = createToolkit<Toy>({
      tools: ['a', 'b', 'c', 'd'].map((id) => make(id)),
      core: ['a'],
      groups: [
        { id: 'g1', title: 'G', summary: 's', guidance: 'g', tools: ['b'] },
        { id: 'g2', title: 'G2', summary: 's2', guidance: 'g', tools: ['c'] },
      ],
    });
    const [s, t] = k.disclosureTools();
    expect(JSON.stringify(s?.definition.inputSchema)).toBe(
      '{"type":"object","properties":{"skills":{"type":"array","items":{"type":"string","enum":["g1","g2"]},"description":"Skill names."}},"required":["skills"]}',
    );
    expect(JSON.stringify(t?.definition.inputSchema)).toBe(
      '{"type":"object","properties":{"tools":{"type":"array","items":{"type":"string","enum":["d"]},"description":"Tool names."}},"required":["tools"]}',
    );
    const [bareSkills, bareTools] = createToolkit<Toy>({
      tools: [make('a')],
      core: ['a'],
      groups: [],
    }).disclosureTools();
    expect(JSON.stringify(bareSkills?.definition.inputSchema)).toBe(
      '{"type":"object","properties":{"skills":{"type":"array","items":{"type":"string"},"description":"Skill names."}},"required":["skills"]}',
    );
    expect(JSON.stringify(bareTools?.definition.inputSchema)).toBe(
      '{"type":"object","properties":{"tools":{"type":"array","items":{"type":"string"},"description":"Tool names."}},"required":["tools"]}',
    );
  });

  it('the disclosure tools take the modules the groups’ availability reads', () => {
    expect(skills?.modules).toEqual(['disclosure', 'extra']);
    expect(tools?.modules).toEqual(['disclosure']);
  });
});

describe('disclosing a group', () => {
  it('contributes its tools’ definitions (not its extras) and its guidance, and says so', async () => {
    const r = await disclose('disclose_skills', { skills: ['volumes'] });
    expect(r.status).toBe('succeeded');
    expect(r.text).toBe(
      'Skill volumes loaded; its tools: load_volume, set_volume_look; by name when needed: remove_volume_look. Its guidance follows.',
    );
    expect(definedIn(r.result.contributions)).toEqual(['load_volume', 'set_volume_look']);
    expect(textsIn(r.result.contributions)).toEqual([kit.groupText(VOLUMES)]);
    expect(kit.groupText(VOLUMES)).toBe(
      `${VOLUMES.guidance}\n\nLoaded with this skill: load_volume, set_volume_look. By name when needed (disclose_tools, or call it): remove_volume_look.`,
    );
    // Definitions are loaded, not deferred.
    const first = r.result.contributions?.[0];
    expect(first?.role === 'system' && first.parts.every((p) => p.kind === 'tool-definition' && !p.deferred)).toBe(
      true,
    );
  });

  it('several at once, any case, and a bare string for a list', async () => {
    const r = await disclose('disclose_skills', { skills: ['Volumes', ' PDBE '] });
    expect(r.text).toContain('Skill volumes loaded');
    expect(r.text).toContain('pdbe loaded; its tools: invoke.');
    expect(definedIn(r.result.contributions)).toEqual(['load_volume', 'set_volume_look', 'invoke']);
    expect((await disclose('disclose_skills', { skills: 'Volumes' })).text).toContain('Skill volumes loaded');
  });

  it('once: a group whose guidance is already in the conversation is said to be loaded, and not given again', async () => {
    const r = await disclose('disclose_skills', { skills: ['volumes', 'notes'] }, { keys: [kit.groupKey(VOLUMES)] });
    expect(r.text).toBe(
      'notes loaded. Skill volumes is already loaded: the guidance is earlier in this conversation, and the tools are available. Its guidance follows.',
    );
    expect(textsIn(r.result.contributions)).toEqual(['Notes: be brief.']);
    const both = await disclose(
      'disclose_skills',
      { skills: ['volumes', 'notes'] },
      { keys: [kit.groupKey(VOLUMES), kit.groupKey(NOTES)] },
    );
    expect(both.text).toBe(
      'Skills volumes, notes are already loaded: the guidance is earlier in this conversation, and the tools are available.',
    );
    expect(both.result.contributions).toBeUndefined();
  });

  it('a tool the conversation already has is not defined again', async () => {
    const r = await disclose('disclose_skills', { skills: ['volumes'] }, { loaded: ['load_volume'] });
    expect(definedIn(r.result.contributions)).toEqual(['set_volume_look']);
  });

  it('a changed text is given again: the key follows the text', () => {
    const changed = { ...VOLUMES, guidance: `${VOLUMES.guidance} Newly added.` };
    expect(kit.groupKey(changed)).not.toBe(kit.groupKey(VOLUMES));
    expect(kit.groupKey(VOLUMES)).toBe(`skill:volumes@${fingerprint(kit.groupText(VOLUMES))}`);
    expect(kit.groupKey(PDBE)).toBe(`entry:pdbe@${fingerprint(PDBE.text)}`);
  });

  it('an unknown id is answered with the closest ones', async () => {
    const r = await disclose('disclose_skills', { skills: ['volume'] });
    expect(r.status).toBe('failed');
    expect(r.text).toBe('No skill "volume". Skills and other texts, closest first: volumes, notes, pdbe, selections.');
    const mixed = await disclose('disclose_skills', { skills: ['notes', 'selection-x'] });
    expect(mixed.status).toBe('succeeded');
    expect(mixed.text).toContain('notes loaded.');
    expect(mixed.text).toContain('No skill "selection-x". Skills and other texts, closest first: selections,');
  });

  it('a group that cannot be used here says why, and what to do instead', async () => {
    const r = await disclose('disclose_skills', { skills: ['selections'] });
    expect(r).toMatchObject({
      status: 'failed',
      text: 'Selections cannot be evaluated here (no analysis). Use targets.',
    });
    const ok = await disclose('disclose_skills', { skills: ['selections'] }, {}, toy({ extra: { name: 'analysis' } }));
    expect(ok.status).toBe('succeeded');
    expect(definedIn(ok.result.contributions)).toEqual(['query_reference']);
  });

  it('what is refused rides as a warning beside what loaded', async () => {
    const r = await disclose('disclose_skills', { skills: ['volumes', 'selections', 'ghost'] });
    expect(r.status).toBe('succeeded');
    expect(r.text).toContain('Skill volumes loaded');
    expect(r.text).toContain(
      'Warnings (fix if it affects the request, else tell the user): No skill "ghost". Skills and other texts, closest first:',
    );
    expect(r.text).toContain('Selections cannot be evaluated here (no analysis). Use targets.');
  });

  it('reads what is loaded from the disclosure module of the API', async () => {
    const api = toy();
    const asked: string[] = [];
    api.disclosure = {
      loaded: async () => (asked.push('loaded'), []),
      skills: async () => (asked.push('skills'), []),
    };
    await call(kit.disclosureTools()[0] as ToolRegistration<Toy>, api, { skills: ['volumes'] });
    expect(asked.sort()).toEqual(['loaded', 'skills']);
  });
});

describe('disclosing tools by name', () => {
  it('loads singles and extras, several at once; the core and the loaded are already available', async () => {
    const r = await disclose(
      'disclose_tools',
      { tools: ['remove_scene', 'move_scene', 'remove_volume_look', 'get_state', 'remove_scene'] },
      { loaded: ['move_scene'] },
    );
    expect(r.text).toBe('Loaded: remove_scene, remove_volume_look. Already available: move_scene, get_state.');
    expect(definedIn(r.result.contributions)).toEqual(['remove_scene', 'remove_volume_look']);
  });

  it('says what cannot be loaded that way, and where it is', async () => {
    const r = await disclose('disclose_tools', {
      tools: ['remove_scene', 'load_volume', 'invoke', 'move_sene', 'zzzz'],
    });
    expect(r.status).toBe('succeeded');
    expect(r.text).toContain('Loaded: remove_scene.');
    expect(r.text).toContain('load_volume belongs to the volumes skill: disclose_skills(["volumes"])');
    expect(r.text).toContain('invoke comes with pdbe: disclose_skills(["pdbe"])');
    expect(r.text).toMatch(/no tool "move_sene" \(closest: move_scene,/);
    expect(r.text).toContain('no tool "zzzz"');
    expect(definedIn(r.result.contributions)).toEqual(['remove_scene']);
  });

  it('an entry’s tool that the conversation has is available', async () => {
    const r = await disclose('disclose_tools', { tools: ['invoke'] }, { loaded: ['invoke'] });
    expect(r.text).toBe('Already available: invoke.');
  });

  it('nothing loaded is a failure that says what to name; none named too', async () => {
    const r = await disclose('disclose_tools', { tools: ['load_volume'] });
    expect(r).toMatchObject({
      status: 'failed',
      text: 'load_volume belongs to the volumes skill: disclose_skills(["volumes"]).',
    });
    expect((await disclose('disclose_tools', { tools: [] })).text).toBe(
      'Name the tools to load, e.g. tools: ["remove_scene"].',
    );
    expect((await disclose('disclose_tools', { tools: 'a, b' })).status).toBe('failed');
  });
});

describe('what of the conversation is a group’s text', () => {
  it('finds the keys of the groups and entries whose text is among the messages, as it is now', () => {
    const text = (t: string): Message => ({ role: 'user', parts: [{ kind: 'text', text: t }] });
    const messages: Message[] = [
      text('hello'),
      kit.groupEntry(VOLUMES),
      { role: 'assistant', parts: [{ kind: 'text', text: kit.groupText(SELECTIONS) }] },
      kit.groupEntry(PDBE),
      text(`${NOTES.text} Old.`),
    ];
    expect(kit.groupsIn(messages)).toEqual([kit.groupKey(VOLUMES), kit.groupKey(PDBE)]);
    const k = createToolkit<Toy>(spec({ groups: [{ ...VOLUMES, guidance: 'Changed.' }, SELECTIONS] }));
    // The guidance given earlier is not the guidance now.
    expect(k.groupsIn(messages)).toEqual([k.groupKey(PDBE)]);
  });

  it('a group or an entry as a context entry is its text as a user message', () => {
    expect(kit.groupEntry(NOTES)).toEqual({
      role: 'user',
      parts: [{ kind: 'text', text: 'Notes: be brief.' }],
      label: { kind: LABEL_ENTRY, key: kit.groupKey(NOTES) },
    });
    expect(kit.groupEntry(VOLUMES)).toEqual({
      role: 'user',
      parts: [{ kind: 'text', text: kit.groupText(VOLUMES) }],
      label: { kind: LABEL_SKILL, key: kit.groupKey(VOLUMES) },
    });
  });

  it('reads a labelled entry by its label and one written before labels by its text', () => {
    const unlabelled = ({ role, parts }: Message): Message => ({ role, parts }) as Message;
    // By its label: the text is not what is matched.
    const relabelled: Message = { ...kit.groupEntry(VOLUMES), parts: [{ kind: 'text', text: 'reworded' }] } as Message;
    expect(kit.groupsIn([relabelled])).toEqual([kit.groupKey(VOLUMES)]);
    // A label of another kind, or with another key, is not a group's, even with the group's text.
    const other: Message = { ...kit.groupEntry(VOLUMES), label: { kind: 'host:note' } } as Message;
    expect(kit.groupsIn([other])).toEqual([]);
    // Before labels: by the exact text.
    expect(kit.groupsIn([unlabelled(kit.groupEntry(VOLUMES))])).toEqual([kit.groupKey(VOLUMES)]);
  });

  it('definitionEntries is one system message of tool-definition parts', () => {
    const [a, b] = kit.allDefinitions();
    expect(definitionEntries([])).toEqual([]);
    expect(definitionEntries([a, b] as never, true)).toEqual([
      {
        role: 'system',
        parts: [
          { kind: 'tool-definition', definition: a, deferred: true },
          { kind: 'tool-definition', definition: b, deferred: true },
        ],
        label: { kind: LABEL_DEFINITIONS },
      },
    ]);
    expect(kit.definitionEntries([a] as never)).toEqual(definitionEntries([a] as never));
  });
});

describe('closest ids', () => {
  it('lists the known ones closest to what was asked first, a few at most', () => {
    expect(closestHint('skills', 'volume', ['selections', 'volumes', 'shapes'])).toBe(
      'Skills, closest first: volumes, shapes, selections.',
    );
    expect(closestHint('skills', 'x', ['a', 'b', 'c'], 2)).toBe('Skills, closest first: a, b (and 1 more).');
    expect(closestHint('skills', 'x', [])).toBe('There are no skills yet.');
  });
});
