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
import { createViewerAgent } from '@molstar/prism-viewer/agent/session';
import { DefaultPrismViewerSpec } from '@molstar/prism-viewer/default-spec';
import { createPrismViewerHeadless } from '@molstar/prism-viewer/headless';
import type { PrismViewer } from '@molstar/prism-viewer/viewer';
import { afterEach, describe, expect, it } from 'vitest';
import { findLimit, InspectorModel } from '../../inspector/model.js';
import { createInspectorSources } from '../inspector-sources.js';

// The Inspector over the viewer's agent (agent.md §1) on a headless viewer with a scripted model: no network.

const cbs = () =>
  new File(
    [fs.readFileSync(path.resolve(import.meta.dirname, '../../../../../../data/examples/1cbs_updated.cif'))],
    '1cbs.cif',
  );

const created: PrismViewer[] = [];
afterEach(async () => {
  for (const v of created.splice(0)) await v.dispose();
});

async function until(check: () => boolean, ms = 3000): Promise<void> {
  const t0 = Date.now();
  while (!check()) {
    if (Date.now() - t0 > ms) throw new Error(`timed out after ${ms} ms`);
    await new Promise((r) => setTimeout(r, 5));
  }
}

const call = (tool: string, args: Record<string, unknown> = {}) => ({ tool, arguments: args as never });

async function agentOn(script: ModelStep[]) {
  const v = await createPrismViewerHeadless({ spec: DefaultPrismViewerSpec() });
  created.push(v);
  const r = await v.api.load.files({ files: [cbs()] });
  if (!r.ok) throw new Error(r.message);
  const model = scriptedModel(script);
  const agent = createViewerAgent({
    viewer: v,
    models: {
      adapterFor: (p) => (p === 'scripted' ? model.adapter : undefined),
      default: { provider: 'scripted', model: 'scripted' },
    },
  });
  return { v, agent };
}

describe('the Inspector over the viewer agent', () => {
  it('reads the trace of a turn with its calls, lists the tools and skills with what is loaded, shows the limits', async () => {
    const { agent } = await agentOn([
      { calls: [call('disclose_skills', { skills: ['styles'] })] },
      { calls: [call('apply_style', { style: 'surface' })] },
      { text: 'A surface.' },
    ]);
    const sources = createInspectorSources(agent);
    const inspector = new InspectorModel(sources);

    // Before a turn: no trace (the Trace's empty state), the toolkit's tools, the session's limits.
    await sources.refresh();
    expect(sources.trace.value).toBeUndefined();
    const before = sources.tools.value;
    const core = before.groups.find((g) => g.kind === 'core')!;
    expect(core.tools.map((t) => t.name)).toEqual(expect.arrayContaining(['load_entry', 'disclose_skills']));
    const styles = () => sources.tools.value.groups.find((g) => g.id === 'styles')!;
    expect(styles()).toMatchObject({ kind: 'skill', title: 'Styles' });
    expect(styles().guidance).toBeTruthy();
    expect(styles().tools.every((t) => t.state !== 'loaded')).toBe(true);

    const result = await agent.turn('show a surface');
    expect(result.status).toBe('completed');
    // The harness's changes are read as they come (no refresh asked).
    await until(() => sources.trace.value?.turns.length === 1);
    await sources.refresh();
    const turn = sources.trace.value!.turns[0];
    expect(turn).toMatchObject({ prompt: 'show a surface', status: 'done' });
    const calls = turn.steps.flatMap((s) => (s.kind === 'tool' ? [[s.name, s.status]] : []));
    expect(calls).toEqual([
      ['disclose_skills', 'succeeded'],
      ['apply_style', 'succeeded'],
    ]);
    expect(turn.steps.filter((s) => s.kind === 'model')).toHaveLength(3);

    // The skill the turn disclosed is loaded now, its tools from turn 1; the core from the start.
    expect(styles().tools.find((t) => t.name === 'remove_component')).toMatchObject({
      state: 'loaded',
      loadedAtTurn: 1,
    });
    const coreNow = sources.tools.value.groups.find((g) => g.kind === 'core')!;
    expect(coreNow.tools.find((t) => t.name === 'apply_style')).toMatchObject({ state: 'loaded' });
    expect(coreNow.tools.find((t) => t.name === 'apply_style')?.loadedAtTurn).toBeUndefined();

    // The limits: the turn budget with the latest turn's use, the window, read-only.
    const toolCalls = findLimit(sources.limits.value, 'toolCalls')!;
    expect(toolCalls).toMatchObject({ editable: false, used: 2 });
    expect(toolCalls.ceiling).toBeGreaterThanOrEqual(toolCalls.value);
    expect(findLimit(sources.limits.value, 'contextWindow')?.value).toBeGreaterThan(0);
    expect(sources.setLimit).toBeUndefined();
    expect(inspector.editable).toBe(false);
    expect(inspector.setLimit('toolCalls', 10)).toBe(false);

    // The Inspector's rows follow: the turn, its request, its steps.
    const rows: string[] = [];
    const sub = inspector.rows.subscribe((r) => rows.splice(0, rows.length, ...r.map((x) => x.kind)));
    await until(() => rows.length > 0);
    expect(rows.slice(0, 2)).toEqual(['turn', 'request']);
    expect(rows.filter((k) => k === 'step')).toHaveLength(5);
    sub.unsubscribe();

    sources.dispose();
    agent.dispose();
  });
});
