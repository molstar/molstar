/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type ModelStep, scriptedModel } from '@molstar/ai/testing';
import { describe, expect, it } from 'vitest';
import { exampleFile, headlessViewer } from '../../_test/harness.js';
import { type EvalTask, evalReport, runEval } from '../eval.js';
import { FirstScope } from './first-scope.js';

// The first eval (agent.md §6) with a scripted model: the runner, the checks and Undo turn, offline. Each task is run by a
// script that does it right (it passes) and one that does it wrong (its check fails).

const call = (tool: string, args: Record<string, unknown> = {}) => ({ tool, arguments: args as never });
const task = (id: string) => FirstScope.find((t) => t.id === id)!;

async function run(t: EvalTask, script: ModelStep[]) {
  const v = await headlessViewer();
  const model = scriptedModel(script);
  return runEval(v, t, {
    models: {
      providerFor: (p) => (p === 'scripted' ? model.provider : undefined),
      default: { provider: 'scripted', model: 'scripted' },
    },
    file: (name) => exampleFile(name),
  });
}

const Right: Record<string, ModelStep[]> = {
  surface: [{ calls: [call('apply_style', { style: 'surface' })] }, { text: 'Done.' }],
  'ligand-red': [{ calls: [call('color', { target: 'ligand', color: 'red' })] }, { text: 'Done.' }],
  'select-range': [{ calls: [call('select', { target: 'A:10-20' })] }, { text: 'Done.' }],
  'remove-water': [
    { calls: [call('disclose_skills', { skills: ['styles'] })] },
    { calls: [call('remove_component', { component: 'Water' })] },
    { text: 'Done.' },
  ],
};

const Wrong: Record<string, ModelStep[]> = {
  surface: [{ calls: [call('apply_style', { style: 'spacefill' })] }, { text: 'Done.' }],
  'ligand-red': [{ calls: [call('color', { target: 'ligand', color: 'blue' })] }, { text: 'Done.' }],
  'select-range': [{ calls: [call('select', { target: 'A:10-25' })] }, { text: 'Done.' }],
  'remove-water': [{ text: 'I would rather not.' }],
};

describe('the first eval', () => {
  for (const t of FirstScope) {
    it(`passes "${t.id}" done right, the run log with the agent's runs`, async () => {
      const r = await run(t, Right[t.id]);
      expect(r.pass, evalReport(r)).toBe(true);
      expect(r.runs.some((x) => x.origin.kind === 'agent')).toBe(true);
      // The file load is the person's (the host's), first.
      expect(r.runs[0].operation).toBe('load.files');
    });

    it(`fails "${t.id}" done wrong`, async () => {
      const r = await run(t, Wrong[t.id]);
      expect(r.pass, evalReport(r)).toBe(false);
      expect(r.checks.some((c) => !c.pass && !c.name.startsWith('turn'))).toBe(true);
    });
  }

  it('answers approvals as the task says, and counts them', async () => {
    const r = await run(task('remove-water'), Right['remove-water']);
    expect(r.turns[0].approvals).toBe(1);
    const denied = await run({ ...task('remove-water'), approve: 'deny' }, Right['remove-water']);
    expect(denied.pass).toBe(false);
  });
});
