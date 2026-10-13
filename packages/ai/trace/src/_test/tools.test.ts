/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import type { ToolDefinition } from '@molstar/ai-agent';
import { definitionTokens, estimatedTokens, type ToolEntry, toolsView } from '../index.js';

const def = (id: string, description = `The ${id} tool.`): ToolDefinition => ({
  tool: { id },
  name: id,
  description,
  inputSchema: { type: 'object', properties: { target: { type: 'string' } } },
});

const entries: ToolEntry[] = [
  { definition: def('get_state'), placement: { kind: 'core' }, read: true, modules: ['state'] },
  { definition: def('set_value'), placement: { kind: 'core' }, approval: true, modules: ['state'] },
  { definition: def('disclose_skills'), placement: { kind: 'core' } },
  { definition: def('load_volume'), placement: { kind: 'skill', skill: 'volumes' } },
  { definition: def('volume_info'), placement: { kind: 'skill', skill: 'volumes', extra: true }, read: true },
  { definition: def('select_near'), placement: { kind: 'skill', skill: 'selections' }, enabled: false },
  { definition: def('web_get'), placement: { kind: 'single' }, unavailable: 'The host has no web access.' },
];

describe('toolsView (analysis/viewer-mode.md §5)', () => {
  it('groups tools as the model is offered them: core, each skill, loaded by name', () => {
    const view = toolsView({
      tools: entries,
      skills: [
        { id: 'volumes', title: 'Density maps', summary: 'maps', guidance: 'Skill volumes: load a map.' },
        { id: 'selections', title: 'Selections', unavailable: 'Selections cannot be evaluated here.' },
      ],
    });
    expect(view.groups.map((g) => [g.id, g.kind, g.title])).toEqual([
      ['core', 'core', 'Core tools'],
      ['volumes', 'skill', 'Density maps'],
      ['selections', 'skill', 'Selections'],
      ['single', 'single', 'Tools loaded by name'],
    ]);
    expect(view.groups[0]?.tools.map((t) => t.id)).toEqual(['get_state', 'set_value', 'disclose_skills']);
    const volumes = view.groups[1];
    expect(volumes).toMatchObject({ summary: 'maps', guidance: 'Skill volumes: load a map.' });
    expect(volumes?.guidanceTokens).toBe(estimatedTokens('Skill volumes: load a map.'));
    // A skill’s extras come after its tools.
    expect(volumes?.tools.map((t) => [t.id, t.extra])).toEqual([
      ['load_volume', undefined],
      ['volume_info', true],
    ]);
    expect(view.groups[2]?.unavailable).toBe('Selections cannot be evaluated here.');
  });

  it('carries what the host knows of each tool: reads, approval, enabled, modules, why unavailable', () => {
    const view = toolsView({
      tools: entries,
      skills: [{ id: 'selections', title: 'Selections', unavailable: 'Selections cannot be evaluated here.' }],
    });
    const byId = new Map(view.groups.flatMap((g) => g.tools).map((t) => [t.id, t]));
    expect(byId.get('get_state')).toMatchObject({ read: true, approval: false, enabled: true, modules: ['state'] });
    expect(byId.get('set_value')).toMatchObject({ read: false, approval: true });
    expect(byId.get('select_near')).toMatchObject({
      enabled: false,
      unavailable: 'Selections cannot be evaluated here.',
    });
    expect(byId.get('web_get')?.unavailable).toBe('The host has no web access.');
    expect(byId.get('get_state')?.unavailable).toBeUndefined();
    expect(byId.get('load_volume')?.modules).toBeUndefined();
    expect(byId.get('get_state')?.input).toEqual({ type: 'object', properties: { target: { type: 'string' } } });
  });

  it('sizes a definition in tokens as the trace does (characters / 4 of its JSON)', () => {
    const view = toolsView({ tools: entries });
    const tool = view.groups[0]?.tools[0];
    expect(tool?.tokens).toBe(definitionTokens(def('get_state')));
    expect(tool?.tokens).toBe(Math.round(JSON.stringify(def('get_state')).length / 4));
  });

  it('says whether each tool is loaded now (and since when), declared, or only offered', () => {
    const view = toolsView({
      tools: entries,
      skills: [{ id: 'volumes', title: 'Density maps' }],
      loaded: [{ id: 'get_state' }, { id: 'set_value' }, { id: 'load_volume', turn: 3 }],
      declared: ['volume_info'],
    });
    const byId = new Map(view.groups.flatMap((g) => g.tools).map((t) => [t.id, t]));
    expect(byId.get('get_state')).toMatchObject({ state: 'loaded' });
    expect(byId.get('get_state')?.loadedAtTurn).toBeUndefined();
    expect(byId.get('load_volume')).toMatchObject({ state: 'loaded', loadedAtTurn: 3 });
    expect(byId.get('volume_info')).toMatchObject({ state: 'deferred' });
    expect(byId.get('disclose_skills')).toMatchObject({ state: 'offered' });
  });

  it('shows a skill a tool names even when the host gave no group for it, and leaves out empty groups', () => {
    const view = toolsView({ tools: [{ definition: def('a'), placement: { kind: 'skill', skill: 'mystery' } }] });
    expect(view.groups.map((g) => [g.id, g.title, g.tools.length])).toEqual([['mystery', 'mystery', 1]]);
    expect(toolsView({ tools: [] }).groups).toEqual([]);
    expect(
      toolsView({
        tools: [{ definition: def('a'), placement: { kind: 'core' } }],
        titles: { core: 'Always on' },
      }).groups[0]?.title,
    ).toBe('Always on');
  });

  it('takes a trace’s loaded list as it is', () => {
    const view = toolsView({
      tools: entries,
      loaded: [{ id: 'get_state', turn: undefined }],
    });
    expect(view.groups[0]?.tools[0]?.state).toBe('loaded');
  });
});
