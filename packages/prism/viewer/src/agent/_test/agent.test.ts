/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type ModelStep, scriptedModel } from '@molstar/ai-conformance/testing';
import { describe, expect, it } from 'vitest';
import { cbs, expectOk, headlessViewer } from '../../_test/harness.js';
import type { PrismViewer } from '../../viewer.js';
import { createViewerAgent } from '../session.js';

// The viewer's agent (agent.md, plan 4.1) on a headless viewer with a scripted model: no network, no real model.

const call = (tool: string, args: Record<string, unknown> = {}) => ({ tool, arguments: args as never });

function agentOn(v: PrismViewer, script: ModelStep[]) {
  const model = scriptedModel(script);
  const agent = createViewerAgent({
    viewer: v,
    models: {
      adapterFor: (p) => (p === 'scripted' ? model.adapter : undefined),
      default: { provider: 'scripted', model: 'scripted' },
    },
  });
  return { agent, model };
}

const representationTypes = (v: PrismViewer) =>
  v.plugin.managers.structure.hierarchy.current.structures[0].components.flatMap((c) =>
    c.representations.map((r) => (r.cell.transform.params as { type: { name: string } }).type.name),
  );

describe('the viewer agent', () => {
  it('runs core tools as the agent, the turn as one history step that Undo turn takes back', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    const before = representationTypes(v);
    const { agent, model } = agentOn(v, [
      { calls: [call('apply_style', { style: 'surface' }), call('color', { target: 'ligand', color: 'red' })] },
      { text: 'Done: a surface, the ligand red.' },
    ]);
    const r = await agent.turn('show it as a surface, the ligand red');
    expect(r.status).toBe('completed');
    expect(
      r.tools.map((t) => [t.name, t.status]),
      JSON.stringify(r.tools.map((t) => t.text)),
    ).toEqual([
      ['apply_style', 'succeeded'],
      ['color', 'succeeded'],
    ]);
    expect(representationTypes(v)).not.toEqual(before);
    // One step for the turn, with the agent's origin.
    const steps = v.managers.history.steps.value;
    expect(steps[0].origin).toMatchObject({ kind: 'agent', turn: agent.lastTurn });
    // The request carried the scene and the focus.
    const sent = JSON.stringify(model.projections[0]);
    expect(sent).toMatch(/structure \\"1CBS\\"/);
    expect(sent).toMatch(/Selected: nothing/);
    // The turn is in the session's record, settled, with its reply (journal.md §2).
    await v.record.flush();
    const turnItem = (await v.record.read(10)).items.find((i) => i.entry.kind === 'turn');
    expect(turnItem?.entry).toMatchObject({ prompt: 'show it as a surface, the ligand red', outcome: 'done' });
    expectOk(await v.api.history.undoTurn({ turn: agent.lastTurn! }));
    expect(representationTypes(v)).toEqual(before);
    agent.dispose();
  });

  it('tells the agent what the person did since its last turn', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    const { agent, model } = agentOn(v, [{ text: 'Hi.' }, { text: 'I see.' }]);
    await agent.turn('hello');
    expectOk(await v.api.as({ kind: 'person' }).styles.apply({ style: 'spacefill' }));
    await agent.turn('what changed?');
    expect(JSON.stringify(model.projections.at(-1))).toMatch(/Notes since your last turn/);
    agent.dispose();
  });

  it('loads a lazy feature when its skill is disclosed, and runs its tool in the same turn', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    expect(v.extensions.status('superposition')?.value).toBe('stub');
    const { agent } = agentOn(v, [
      { calls: [call('disclose_skills', { skills: ['superpose'] })] },
      { calls: [call('superpose_chains', {})] },
      { text: 'Tried.' },
    ]);
    const r = await agent.turn('superpose them');
    expect(v.extensions.status('superposition')?.value).toBe('loaded');
    expect(r.tools.map((t) => t.name)).toEqual(['disclose_skills', 'superpose_chains']);
    // One structure: it says why it cannot, as a tool result (not a crash).
    expect(r.tools[1].status).toBe('failed');
    agent.dispose();
  });

  it('asks before a destructive call', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    const { agent } = agentOn(v, [
      { calls: [call('disclose_skills', { skills: ['styles'] })] },
      { calls: [call('remove_component', { component: 'Ligand' })] },
      { text: 'Removed.' },
    ]);
    const r = await agent.turn('remove the ligand');
    expect(r.status).toBe('awaiting-approval');
    const pending = r.approval![0];
    const done = await agent.resume({ [pending.callId]: 'allow' });
    expect(done?.status).toBe('completed');
    agent.dispose();
  });
});

describe('the viewer agent while undo is off', () => {
  it('asks before the first change of the turn, then lets the rest go', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    v.managers.history.setCheckpointing(false);
    const { agent } = agentOn(v, [
      { calls: [call('apply_style', { style: 'surface' })] },
      { calls: [call('color', { target: 'ligand', color: 'red' })] },
      { text: 'Done.' },
    ]);
    const r = await agent.turn('surface, ligand red');
    expect(r.status).toBe('awaiting-approval');
    expect(r.approval?.[0].tool).toBe('apply_style');
    const done = await agent.resume({ [r.approval![0].callId]: 'allow' });
    expect(done?.status).toBe('completed');
    expect(done?.tools.map((t) => [t.name, t.status])).toContainEqual(['color', 'succeeded']);
    agent.dispose();
  });
});

describe('the viewer agent at G6', () => {
  it('ends a parked turn on stop: settled in the record, its group closed; and refuses a new turn while one waits', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    const { agent } = agentOn(v, [
      { calls: [call('disclose_skills', { skills: ['styles'] })] },
      { calls: [call('remove_component', { component: 'Ligand' })] },
      { text: 'Removed.' },
    ]);
    const r = await agent.turn('remove the ligand');
    expect(r.status).toBe('awaiting-approval');
    expect(agent.parked).toBe(true);
    await expect(agent.turn('something else')).rejects.toThrow(/waits for your approval/);
    await agent.stop();
    expect(agent.parked).toBe(false);
    await v.record.flush();
    const turn = (await v.record.read(20)).items.find((i) => i.entry.kind === 'turn');
    expect(turn?.entry).toMatchObject({ outcome: 'stopped' });
    expect(turn?.open).toBe(false);
    // Nothing left waiting: resuming has nothing to go on with.
    expect(await agent.resume({ [r.approval![0].callId]: 'allow' })).toBeUndefined();
    agent.dispose();
  });

  it('asks before loading from an address it cannot check, whatever its form', async () => {
    const v = await headlessViewer();
    for (const url of [' https://evil.example/x.cif', '//evil.example/x.cif', 'evil.example/x.cif']) {
      const { agent } = agentOn(v, [{ calls: [call('load_url', { url })] }, { text: 'Loaded.' }]);
      const r = await agent.turn('load it');
      expect(r.status, url).toBe('awaiting-approval');
      await agent.stop();
      agent.dispose();
    }
  });

  it('describes its tools by what this viewer can do, not by the scene of the moment', async () => {
    const v = await headlessViewer();
    const { agent } = agentOn(v, []);
    const description = (id: string) => agent.toolkit.registration(id)?.definition.description ?? '';
    // Nothing is loaded yet: the tools are offered as they are (the call says when it cannot run).
    for (const id of ['reset_camera', 'select', 'focus']) {
      expect(description(id), id).toMatch(/\w/);
      expect(description(id), id).not.toMatch(/Not now|Nothing is loaded/);
    }
    agent.dispose();
  });

  it("takes an extension's skills, over its own operations", async () => {
    const v = await headlessViewer({
      extensions: [
        {
          id: 'probe',
          title: 'Probe',
          skills: [
            {
              id: 'probe',
              title: 'Probe',
              summary: 'a probe',
              guidance: 'Skill probe.',
              tools: [{ id: 'probe_selection', operation: 'select.get' }],
            },
          ],
        },
      ],
    });
    const { agent } = agentOn(v, []);
    expect(agent.toolkit.groups.map((g) => g.id)).toContain('probe');
    agent.dispose();
  });

  it("undoes its previous turn with undo_turn, never the person's later steps", async () => {
    const v = await headlessViewer({ files: [cbs()] });
    const before = representationTypes(v);
    const { agent } = agentOn(v, [
      { calls: [call('apply_style', { style: 'surface' })] },
      { text: 'A surface.' },
      { calls: [call('undo_turn')] },
      { text: 'Undone.' },
    ]);
    await agent.turn('surface');
    expect(representationTypes(v)).not.toEqual(before);
    const r = await agent.turn('undo that');
    expect(r.tools.map((t) => [t.name, t.status])).toEqual([['undo_turn', 'succeeded']]);
    expect(representationTypes(v)).toEqual(before);
    agent.dispose();
  });
});
