/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { ModelChoice, ModelSource, TurnResult } from '@molstar/ai';
import type { RunRecord } from '../api/run.js';
import type { PrismViewer } from '../viewer.js';
import { createAccount } from './account.js';
import { createViewerAgent, type ViewerAgent } from './session.js';

// Evals (api.md §11, agent.md §6): a task runs the agent over a headless viewer and is checked on the scene (the state
// the tools changed), not on the model's words. The result carries the run log (the API's runs with their arguments and
// origins), which replays in the app to show what the agent did. DOM-free: the host makes the viewer and gives the
// models (a scripted one in unit tests, real ones in `*.live.ts`).

/** One check of the scene after the task's turns. */
export interface EvalCheck {
  readonly name: string;
  readonly pass: boolean;
  /** What was found, when it failed. */
  readonly detail?: string;
}

export interface EvalTask {
  readonly id: string;
  /** What it tests, for the report. */
  readonly title: string;
  /** Example files the host loads before the turns (by name; the host reads them). */
  readonly files?: readonly string[];
  /** The person's requests, one turn each. */
  readonly prompts: readonly string[];
  /** How the person answers an approval (Allow by default). */
  readonly approve?: 'allow' | 'deny';
  /** The checks on the scene after the turns. */
  check(viewer: PrismViewer): readonly EvalCheck[] | Promise<readonly EvalCheck[]>;
  /** A turn takes back as one step: after the checks, Undo turn of each turn (newest first) restores the outline. */
  readonly undoes?: boolean;
}

export interface EvalTurn {
  readonly prompt: string;
  readonly status: TurnResult['status'];
  readonly tools: readonly { readonly name: string; readonly status: string }[];
  /** Approvals asked (and answered as the task says). */
  readonly approvals: number;
  readonly usage?: { readonly input: number; readonly output: number };
  readonly error?: string;
}

export interface EvalResult {
  readonly task: string;
  readonly pass: boolean;
  readonly checks: readonly EvalCheck[];
  readonly turns: readonly EvalTurn[];
  /** The API's runs during the task (the person's file loads first), in order. */
  readonly runs: readonly RunRecord[];
  readonly ms: number;
}

export interface EvalOptions<M extends ModelChoice = ModelChoice> {
  readonly models: ModelSource<M>;
  /** The model to ask (the source's default when not given). */
  readonly model?: M;
  /** An example file by name, for the task's `files`. */
  readonly file?: (name: string) => File | Promise<File>;
  /** Approval rounds a turn may ask before the eval gives up on it (5). */
  readonly maxApprovals?: number;
}

/** Runs a task on a viewer (made for it, headless) and checks the scene. */
export async function runEval<M extends ModelChoice>(
  viewer: PrismViewer,
  task: EvalTask,
  o: EvalOptions<M>,
): Promise<EvalResult> {
  const started = Date.now();
  const runs: RunRecord[] = [];
  const sub = viewer.api.events.finished.subscribe((r) => runs.push(r));
  const turns: EvalTurn[] = [];
  const turnIds: string[] = [];
  const checks: EvalCheck[] = [];
  let agent: ViewerAgent<M> | undefined;
  try {
    for (const name of task.files ?? []) {
      if (!o.file) throw new Error(`The task '${task.id}' needs files, and no file reader was given.`);
      const r = await viewer.api.load.files({ files: [await o.file(name)] });
      if (!r.ok) throw new Error(`Could not load ${name}: ${r.message}`);
    }
    // After the files, and the sequence's document for them (its runtime loads on its own: the structure tree's
    // Sequence rows read it, sequence.md §8.1), the agent starts on the scene the task gives it.
    await sequenceSettled(viewer);
    agent = createViewerAgent({ viewer, models: o.models });
    const before = sceneOf(viewer);
    for (const prompt of task.prompts) {
      const tools: { name: string; status: string }[] = [];
      const hooks = { onTool: (e: { name: string; status: string }) => void tools.push(e) };
      let r: TurnResult | undefined = await agent.turn(prompt, hooks, o.model);
      let approvals = 0;
      let rounds = 0;
      while (r?.status === 'awaiting-approval' && r.approval?.length && rounds < (o.maxApprovals ?? 5)) {
        rounds++;
        approvals += r.approval.length;
        const decision = task.approve ?? 'allow';
        r = await agent.resume({ decisions: Object.fromEntries(r.approval.map((p) => [p.callId, decision])), hooks });
      }
      // Still waiting (too many rounds, or nothing to answer): the turn is stopped, and fails.
      const gaveUp = r?.status === 'awaiting-approval';
      if (gaveUp) await agent.stop();
      const status = gaveUp ? 'stopped' : (r?.status ?? 'failed');
      turns.push({
        prompt,
        status,
        // Each call's end (a call also reports its start and progress).
        tools: tools.filter((t) => t.status !== 'started' && t.status !== 'progress'),
        approvals,
        ...(r?.usage ? { usage: { input: r.usage.input, output: r.usage.output } } : {}),
        ...(gaveUp
          ? { error: `still waiting for an approval after ${rounds} rounds` }
          : r?.error
            ? { error: r.error }
            : {}),
      });
      if (agent.lastTurn && !turnIds.includes(agent.lastTurn)) turnIds.push(agent.lastTurn);
    }
    for (const t of turns) {
      checks.push({
        name: `turn "${t.prompt}" completed`,
        pass: t.status === 'completed',
        ...(t.status === 'completed' ? {} : { detail: t.error ?? t.status }),
      });
    }
    try {
      checks.push(...(await task.check(viewer)));
    } catch (e) {
      checks.push({ name: 'the checks ran', pass: false, detail: e instanceof Error ? e.message : String(e) });
    }
    if (task.undoes) {
      for (const id of turnIds.slice().reverse()) {
        const r = await viewer.api.history.undoTurn({ turn: id });
        if (!r.ok) checks.push({ name: `Undo turn ${id}`, pass: false, detail: r.message });
      }
      const after = sceneOf(viewer);
      checks.push({
        name: 'Undo turn restores the scene',
        pass: after === before,
        ...(after === before ? {} : { detail: `before:\n${before}\nafter:\n${after}` }),
      });
    }
    return { task: task.id, pass: checks.every((c) => c.pass), checks, turns, runs, ms: Date.now() - started };
  } finally {
    sub.unsubscribe();
    agent?.dispose();
  }
}

/**
 * The scene as Undo turn must restore it: the agent's outline (agent.md §5), and the data tree's transforms with their
 * params (a color layer, which the outline does not show).
 */
/** The sequence's runtime in and its document made for what is loaded (bounded), then the tree's rows told of it. */
async function sequenceSettled(viewer: PrismViewer, ms = 10_000): Promise<void> {
  const sequence = viewer.views.sequence;
  if (!sequence?.source.structure.value) return;
  const viewers = (await sequence.load()).viewers.get();
  const t0 = Date.now();
  while (!viewers.structure.value && Date.now() - t0 < ms) await new Promise((r) => setTimeout(r, 10));
  // The tree's sequence rows follow once a burst (a microtask), then the scene's sync.
  await new Promise((r) => setTimeout(r, 0));
}

function sceneOf(viewer: PrismViewer): string {
  const transforms = [...viewer.plugin.state.data.cells.values()]
    .map(({ transform: t }) => `${t.ref} ${t.transformer.id} ${JSON.stringify(t.params ?? {})}`)
    .sort();
  return `${createAccount(viewer).outline(1000)}\n${transforms.join('\n')}`;
}

/** A report in a few lines: the task, its checks, its turns' tools and tokens. */
export function evalReport(r: EvalResult): string {
  const lines = [`${r.pass ? 'PASS' : 'FAIL'} ${r.task} (${(r.ms / 1000).toFixed(1)} s)`];
  for (const c of r.checks) lines.push(`  ${c.pass ? '✓' : '✗'} ${c.name}${c.detail ? `: ${c.detail}` : ''}`);
  for (const t of r.turns) {
    const usage = t.usage ? `, ${t.usage.input} in / ${t.usage.output} out` : '';
    lines.push(`  · ${t.tools.map((x) => `${x.name}:${x.status}`).join(' ') || 'no tools'}${usage}`);
  }
  return lines.join('\n');
}
