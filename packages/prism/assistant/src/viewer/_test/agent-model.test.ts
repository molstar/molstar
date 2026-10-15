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
    // The history (the session's record): the turn with its run.
    await v.record.flush();
    await until(() => turnsIn(m.history.value).length === 1);
    const [listed] = turnsIn(m.history.value);
    expect(listed).toMatchObject({ prompt: 'show a surface', status: 'done' });
    expect(listed.runs.map((r) => r.command)).toContain('styles.apply');
    // Undo turn: through the person's history.
    m.undo(t.id);
    await until(() => !!m.turn.value?.undone);
    expect(m.turn.value?.checkpoint).toBe(false);
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
});
