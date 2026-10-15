/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type ModelStep, scriptedModel } from '@molstar/ai-conformance/testing';
import * as fs from 'node:fs';
import * as path from 'node:path';
import { observableState } from '@molstar/prism-ui/core';
import { createViewerAgent } from '@molstar/prism-viewer/agent/session';
import { DefaultPrismViewerSpec } from '@molstar/prism-viewer/default-spec';
import { createPrismViewerHeadless } from '@molstar/prism-viewer/headless';
import type { PrismViewer } from '@molstar/prism-viewer/viewer';
import { afterEach, describe, expect, it } from 'vitest';
import type { AgentTurn, HistoryItem } from '../../model.js';
import { createViewerAgentModel } from '../agent-model.js';
import type { ViewerAISettings } from '../settings-model.js';

// The assistant's model over the viewer's agent (plan 4.4) on a headless viewer with a scripted model: no network.

const cbs = () =>
  new File(
    [fs.readFileSync(path.resolve(import.meta.dirname, '../../../../../../data/examples/1cbs_updated.cif'))],
    '1cbs.cif',
  );

const created: PrismViewer[] = [];
afterEach(async () => {
  for (const v of created.splice(0)) await v.dispose();
});

async function headlessViewer(o: { files: File[] }): Promise<PrismViewer> {
  const v = await createPrismViewerHeadless({ spec: DefaultPrismViewerSpec() });
  created.push(v);
  for (const f of o.files) {
    const r = await v.api.load.files({ files: [f] });
    if (!r.ok) throw new Error(r.message);
  }
  return v;
}

async function until(check: () => boolean, ms = 3000): Promise<void> {
  const t0 = Date.now();
  while (!check()) {
    if (Date.now() - t0 > ms) throw new Error(`timed out after ${ms} ms`);
    await new Promise((r) => setTimeout(r, 5));
  }
}

const call = (tool: string, args: Record<string, unknown> = {}) => ({ tool, arguments: args as never });

function modelOn(v: PrismViewer, script: ModelStep[], chosen = true) {
  const model = scriptedModel(script);
  const setup = observableState(chosen ? { provider: 'scripted', model: 'scripted' } : undefined);
  const settings = { setup: setup.asReadonly } as unknown as ViewerAISettings;
  const agent = createViewerAgent({
    viewer: v,
    models: {
      adapterFor: (p) => (p === 'scripted' ? model.adapter : undefined),
      default: { provider: 'scripted', model: 'scripted' },
    },
  });
  let openedSetup = 0;
  const m = createViewerAgentModel(v, agent, settings, { openSetup: () => openedSetup++ });
  return { m, agent, openedSetup: () => openedSetup };
}

/** What the viewer said to the person since `from` (newest first). */
const said = (v: PrismViewer) => v.managers.notifications.history.value.map((n) => n.text);

const removeLigand: ModelStep[] = [
  { calls: [call('disclose_skills', { skills: ['styles'] })] },
  { calls: [call('remove_component', { component: 'Ligand' })] },
  { text: 'Removed.' },
];

/** The record's turn entries, newest first. */
async function recordedTurns(v: PrismViewer) {
  await v.record.flush();
  return (await v.record.read(50)).items.filter((i) => i.entry.kind === 'turn');
}

const settled = (t: AgentTurn | undefined) => !!t && ['done', 'failed', 'stopped', 'waiting'].includes(t.status);
const turnsIn = (items: readonly HistoryItem[]) =>
  items.filter((i): i is Extract<HistoryItem, { kind: 'turn' }> => i.kind === 'turn').map((i) => i.turn);

describe('the assistant over the viewer agent', () => {
  it('shows a turn as it runs, ends it done with its checkpoint, and lists it in the history', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    const { m, agent } = modelOn(v, [{ calls: [call('apply_style', { style: 'surface' })] }, { text: 'A surface.' }]);
    m.ask('show a surface');
    expect(m.view.value).toBe('card');
    await until(() => settled(m.turn.value));
    const t = m.turn.value!;
    expect(t).toMatchObject({ status: 'done', text: 'A surface.', checkpoint: true, id: agent.lastTurn });
    expect(t.runs.map((r) => [r.command, r.status])).toEqual([['apply_style', 'done']]);
    // The history (the session's record, as History reads it): read while watched.
    const watching = m.history.subscribe();
    await v.record.flush();
    await until(() => turnsIn(m.history.value).length === 1);
    const [listed] = turnsIn(m.history.value);
    expect(listed).toMatchObject({ prompt: 'show a surface', status: 'done' });
    expect(listed.runs.map((r) => r.command)).toContain('styles.apply');
    // Undo turn: through the person's history.
    m.undo(t.id);
    await until(() => !!m.turn.value?.undone);
    expect(m.turn.value?.checkpoint).toBe(false);
    await until(() => turnsIn(m.history.value)[0]?.undone === true);
    watching.unsubscribe();
    m.dispose();
    agent.dispose();
  });

  it('waits for the Allow of a destructive call, then finishes the turn', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    const { m, agent } = modelOn(v, [
      { calls: [call('disclose_skills', { skills: ['styles'] })] },
      { calls: [call('remove_component', { component: 'Ligand' })] },
      { text: 'Removed.' },
    ]);
    m.ask('remove the ligand');
    await until(() => settled(m.turn.value));
    const waiting = m.turn.value!.runs.find((r) => r.status === 'waiting');
    expect(m.turn.value?.status).toBe('waiting');
    expect(waiting).toMatchObject({ command: 'remove_component', approval: true });
    m.approve(waiting!.id, true);
    await until(() => m.turn.value?.status === 'done');
    expect(m.turn.value?.text).toBe('Removed.');
    m.dispose();
    agent.dispose();
  });

  it('opens the setup when no model is chosen, and asks nothing', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    const { m, agent, openedSetup } = modelOn(v, [], false);
    m.ask('hello');
    expect(openedSetup()).toBe(1);
    expect(m.turn.value).toBeUndefined();
    m.dispose();
    agent.dispose();
  });

  it('stops a parked turn: the card, the record and the history group end, and a late Allow does nothing', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    // The scripted model goes on from where the stopped turn left it: the next turn gets the next step.
    const { m, agent } = modelOn(v, [...removeLigand.slice(0, 2), { text: 'Hello again.' }]);
    m.ask('remove the ligand');
    await until(() => m.turn.value?.status === 'waiting');
    const waiting = m.turn.value!.runs.find((r) => r.status === 'waiting')!;
    expect(agent.parked).toBe(true);

    m.stop();
    // At once on the card: stopped, nothing left to answer.
    expect(m.turn.value?.status).toBe('stopped');
    expect(m.turn.value?.runs.some((r) => r.status === 'waiting')).toBe(false);
    await until(() => !agent.parked);
    const [entry] = await recordedTurns(v);
    expect(entry).toMatchObject({ open: false, entry: { outcome: 'stopped' } });

    // Approving after the stop is ignored: the card stays stopped, nothing resumes.
    m.approve(waiting.id, true);
    expect(m.turn.value?.status).toBe('stopped');
    expect(agent.session.busy).toBe(false);

    // The group is closed: the next turn runs and ends on its own.
    m.ask('say hello');
    await until(() => m.turn.value?.status === 'done');
    expect(m.turn.value?.text).toBe('Hello again.');
    expect((await recordedTurns(v)).map((i) => i.open)).toEqual([false, false]);
    m.dispose();
    agent.dispose();
  });

  it('refuses a second ask while a turn waits for an approval, in words, and keeps the waiting card', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    const { m, agent } = modelOn(v, removeLigand);
    m.ask('remove the ligand');
    await until(() => m.turn.value?.status === 'waiting');
    const before = m.turn.value!;
    m.dismiss();

    m.ask('and color it red');
    expect(said(v)[0]).toMatch(/Answer or stop the waiting step first/);
    // The waiting turn stays on the card (shown again), its approval still live.
    expect(m.turn.value).toBe(before);
    expect(m.view.value).toBe('card');
    expect((await recordedTurns(v)).map((i) => i.entry.kind === 'turn' && i.entry.prompt)).toEqual([
      'remove the ligand',
    ]);

    m.approve(before.runs.find((r) => r.status === 'waiting')!.id, true);
    await until(() => m.turn.value?.status === 'done');
    expect(m.turn.value?.text).toBe('Removed.');
    m.dispose();
    agent.dispose();
  });

  it('refuses an ask while a turn runs', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    const { m, agent } = modelOn(v, [{ calls: [call('apply_style', { style: 'surface' })] }, { text: 'A surface.' }]);
    m.ask('show a surface');
    m.ask('and a cartoon');
    expect(said(v)[0]).toMatch(/A turn is running/);
    expect(m.turn.value?.prompt).toBe('show a surface');
    await until(() => m.turn.value?.status === 'done');
    m.dispose();
    agent.dispose();
  });

  it('answers several pending calls one by one, each shown as answered at once', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    const { m, agent } = modelOn(v, [
      { calls: [call('disclose_skills', { skills: ['styles'] })] },
      {
        calls: [call('remove_component', { component: 'Ligand' }), call('remove_component', { component: 'Polymer' })],
      },
      { text: 'Removed one.' },
    ]);
    m.ask('remove the ligand and the polymer');
    await until(() => m.turn.value?.status === 'waiting');
    const asked = m.turn.value!.runs.filter((r) => r.status === 'waiting');
    expect(asked).toHaveLength(2);

    m.approve(asked[0].id, true);
    let runs = m.turn.value!.runs;
    expect(runs.find((r) => r.id === asked[0].id)).toMatchObject({ status: 'running', detail: 'Allowed.' });
    expect(runs.find((r) => r.id === asked[1].id)?.status).toBe('waiting');
    // Still waiting for the other; answering the same one again does nothing.
    expect(m.turn.value?.status).toBe('waiting');
    m.approve(asked[0].id, false);
    expect(m.turn.value!.runs.find((r) => r.id === asked[0].id)?.status).toBe('running');

    m.approve(asked[1].id, false);
    runs = m.turn.value!.runs;
    expect(runs.find((r) => r.id === asked[1].id)?.status).not.toBe('waiting');
    await until(() => m.turn.value?.status === 'done');
    expect(m.turn.value?.text).toBe('Removed one.');
    m.dispose();
    agent.dispose();
  });

  it('holds the checkpoint only while the turn is the newest step, and says why Undo turn cannot run', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    const { m, agent } = modelOn(v, [{ calls: [call('apply_style', { style: 'surface' })] }, { text: 'A surface.' }]);
    m.ask('show a surface');
    await until(() => m.turn.value?.status === 'done');
    expect(m.turn.value?.checkpoint).toBe(true);
    // The person's next step releases it.
    const r = await v.api.as({ kind: 'person' }).styles.apply({ style: 'cartoon' });
    expect(r.ok).toBe(true);
    await until(() => m.turn.value?.checkpoint === false);
    m.undo(m.turn.value!.id);
    await until(() => /later step/.test(said(v)[0] ?? ''));
    expect(m.turn.value?.undone).toBeUndefined();
    m.dispose();
    agent.dispose();
  });

  it('stops a parked turn when it is disposed', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    const { m, agent } = modelOn(v, removeLigand);
    m.ask('remove the ligand');
    await until(() => m.turn.value?.status === 'waiting');
    m.dispose();
    await until(() => !agent.parked);
    const [entry] = await recordedTurns(v);
    expect(entry).toMatchObject({ open: false, entry: { outcome: 'stopped' } });
    agent.dispose();
  });

  it('says a failed turn on its card', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    const { m, agent } = modelOn(v, []);
    const turn = agent.turn.bind(agent);
    (agent as { turn: typeof turn }).turn = () => Promise.reject(new Error('The model is away.'));
    m.ask('hello');
    await until(() => m.turn.value?.status === 'failed');
    expect(m.turn.value?.text).toBe('The model is away.');
    m.dispose();
    agent.dispose();
  });
});
