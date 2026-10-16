/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { createToolBindings, createToolkit, type ToolGroup, type ToolRegistration } from '../index.js';
import { call, definedIn, type Toy, toy, tool } from './helpers.js';

// A dynamic toolkit (decision 25 of the prism-ui spec, `.v6/prism-ui/spec/README.md`; plan track A6): parts added and
// removed after creation, lazy skills that bring their tools when disclosed, the turn's tool set, and the eviction
// hint.

const none = { kind: 'object', fields: {} } as const;
const make = (id: string, approval = false): ToolRegistration<Toy> =>
  tool({
    id,
    description: `${id}. More.`,
    input: none,
    modules: [],
    approval,
    run: () => ({ ok: true, did: `${id} done.` }),
  });

const base = () => createToolkit<Toy>({ tools: [make('get_state')], core: ['get_state'], groups: [] });

describe('a dynamic toolkit', () => {
  it('adds a part, says so, and removes it; the removed tools are the eviction hint', () => {
    const kit = base();
    let changes = 0;
    kit.onChange(() => changes++);
    const group: ToolGroup<Toy> = {
      id: 'anim',
      title: 'Animation',
      summary: 'moves',
      guidance: 'Use spin.',
      tools: ['spin'],
    };
    const remove = kit.add({ tools: [make('spin'), make('zoom')], groups: [group] });
    expect(kit.version).toBe(1);
    expect(changes).toBe(1);
    expect(kit.registration('spin')).toBeDefined();
    expect(kit.placementOf('spin')).toMatchObject({ kind: 'group', group: { id: 'anim' } });
    expect(kit.placementOf('zoom')).toEqual({ kind: 'single' });
    // The catalog in the disclosure tool's description follows.
    expect(kit.disclosureTools()[0].definition.description).toMatch(/"anim" \(moves\)/);
    remove();
    expect(kit.registration('spin')).toBeUndefined();
    expect([...kit.removed()].sort()).toEqual(['spin', 'zoom']);
    expect(kit.version).toBe(2);
    // Added again: no longer evicted.
    kit.add({ tools: [make('spin')] });
    expect(kit.removed()).toEqual(['zoom']);
  });

  it('refuses a part that conflicts, and stays as it was', () => {
    const kit = base();
    expect(() => kit.add({ tools: [make('get_state')] })).toThrow(/registered twice/);
    expect(() => kit.add({ core: ['nope'] })).toThrow(/not registered/);
    expect(kit.version).toBe(0);
    expect(kit.tools.map((t) => t.definition.tool.id)).toEqual(['get_state', 'disclose_skills', 'disclose_tools']);
  });

  it('enables and preauthorizes what a turn may call, the tools lazy skills will bring included', () => {
    const kit = base();
    kit.add({ tools: [make('erase', true)] });
    kit.add({
      groups: [
        {
          id: 'export',
          title: 'Export',
          summary: 'files out',
          guidance: '',
          tools: [],
          pending: [{ id: 'save_image' }, { id: 'upload', approval: true }],
          load: async () => {},
        },
      ],
    });
    const set = kit.toolSet();
    expect(set.enabled.sort()).toEqual([
      'disclose_skills',
      'disclose_tools',
      'erase',
      'get_state',
      'save_image',
      'upload',
    ]);
    expect(set.preauthorized.sort()).toEqual(['disclose_skills', 'disclose_tools', 'get_state', 'save_image']);
  });

  it('loads a lazy skill when it is disclosed: its real group and tools replace the stub, and come with the result', async () => {
    const kit = base();
    let loads = 0;
    const real: ToolGroup<Toy> = {
      id: 'export',
      title: 'Export',
      summary: 'files out',
      guidance: 'Use save_image.',
      tools: ['save_image'],
    };
    let removeStub = () => {};
    removeStub = kit.add({
      groups: [
        {
          id: 'export',
          title: 'Export',
          summary: 'files out',
          guidance: '',
          tools: [],
          pending: [{ id: 'save_image' }],
          load: async () => {
            loads++;
            removeStub();
            kit.add({ tools: [make('save_image')], groups: [real] });
          },
        },
      ],
    });
    const api = toy();
    const disclose = kit.disclosureTools()[0];
    const r = await call(disclose, api, { skills: ['export'] });
    expect(loads).toBe(1);
    expect(r.status).toBe('succeeded');
    expect(r.text).toMatch(/Skill export loaded; its tools: save_image/);
    expect(definedIn(r.result.contributions ?? [])).toEqual(['save_image']);
    // The bindings find the tool that was not there when they were made.
    const bindings = createToolBindings(kit, api);
    expect(bindings.resolve(kit.registration('save_image')!.definition)).toBeDefined();
  });

  it('says when a lazy skill cannot be loaded', async () => {
    const kit = base();
    kit.add({
      groups: [
        {
          id: 'broken',
          title: 'Broken',
          summary: 'fails',
          guidance: '',
          tools: [],
          load: async () => {
            throw new Error('offline');
          },
        },
      ],
    });
    const r = await call(kit.disclosureTools()[0], toy(), { skills: ['broken'] });
    expect(r.status).not.toBe('succeeded');
    expect(r.text).toMatch(/could not be loaded: offline/);
  });
});
