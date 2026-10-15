/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { PendingApproval, ToolEvent, TurnHooks, TurnResult } from '@molstar/ai-session';
import { observableState } from '@molstar/prism-ui/core';
import type { ViewerAgent } from '@molstar/prism-viewer/agent/session';
import type { RecordItem } from '@molstar/prism-viewer/journal/record';
import type { PrismViewer } from '@molstar/prism-viewer/viewer';
import type { AgentModel, AgentTurn, CommandRun, HistoryItem, RunStatus } from '../model.js';
import type { ViewerAISettings } from './settings-model.js';

// The assistant's model over the viewer's agent (plan 4.4, agent.md §1): the latest turn from the session's hooks (its
// text, its calls, waiting for the person), asking with the model chosen in the setup, approvals answered as the turn
// asks them, Undo turn through the viewer's history (the person's API), and the history read from the session's record
// (journal.md §3).

export interface ViewerAgentModel extends AgentModel {
  dispose(): void;
}

export function createViewerAgentModel(
  viewer: PrismViewer,
  agent: ViewerAgent,
  settings: ViewerAISettings,
  options: { openSetup?: () => void } = {},
): ViewerAgentModel {
  const turn = observableState<AgentTurn | undefined>(undefined);
  const view = observableState<'card' | 'line' | 'none'>('none');
  const history = observableState<readonly HistoryItem[]>([]);
  const person = viewer.api.as({ kind: 'person' });
  let pending: PendingApproval[] = [];
  const decided = new Map<string, 'allow' | 'deny'>();

  const update = (patch: Partial<AgentTurn>) => {
    if (turn.value) turn.set({ ...turn.value, ...patch });
  };
  const runOf = (e: ToolEvent): CommandRun => ({
    id: e.id ?? `${e.name}-${e.startedAt}`,
    command: e.name,
    title: e.name.replaceAll('_', ' '),
    source: { kind: 'assistant' },
    status: statusOf(e.status),
    ...(e.read ? { read: true } : {}),
    ...(e.text || e.progress ? { detail: e.text ?? e.progress } : {}),
    time: new Date(e.startedAt),
    ...(e.endedAt !== undefined ? { durationMs: e.endedAt - e.startedAt } : {}),
  });
  const hooks: TurnHooks = {
    onTurn: (e) => {
      if (e.kind === 'begun' && turn.value) update({ id: e.turnId });
    },
    onText: (text) => update({ text }),
    onTool: (e) => {
      const t = turn.value;
      if (!t) return;
      const run = runOf(e);
      const at = t.runs.findIndex((r) => r.id === run.id);
      update({ runs: at < 0 ? [...t.runs, run] : t.runs.map((r, i) => (i === at ? run : r)), status: 'working' });
    },
    onStatus: (e) => {
      if (e.kind === 'model') update({ status: 'thinking' });
      else if (e.kind === 'tools') update({ status: 'working' });
      // `permission` opens for every call with a rule, the policy's yes too: waiting is the turn's `awaiting-approval`.
      else if (e.kind === 'retrying') update({ retries: (turn.value?.retries ?? 0) + 1 });
    },
  };
  const settle = (r: TurnResult | undefined, started: number) => {
    if (!r || !turn.value) return;
    if (r.status === 'awaiting-approval') {
      pending = r.approval ?? [];
      decided.clear();
      const waiting = new Set(pending.map((p) => p.callId));
      update({
        status: 'waiting',
        runs: [
          ...turn.value.runs,
          ...pending
            .filter((p) => !turn.value!.runs.some((x) => x.id === p.callId))
            .map(
              (p): CommandRun => ({
                id: p.callId,
                command: p.tool,
                title: p.tool.replaceAll('_', ' '),
                source: { kind: 'assistant' },
                status: 'waiting',
                approval: true,
                detail: 'Waiting for your Allow.',
                time: new Date(),
              }),
            ),
        ].map((x) => (waiting.has(x.id) ? { ...x, status: 'waiting' as RunStatus, approval: true } : x)),
      });
      return;
    }
    pending = [];
    update({
      status: r.status === 'completed' ? 'done' : r.status === 'stopped' ? 'stopped' : 'failed',
      text: r.text || turn.value.text || (r.error ?? ''),
      checkpoint: viewer.managers.history.steps.value.some((s) => s.origin.turn === turn.value?.id),
      ...(r.usage ? { usage: { input: r.usage.input, output: r.usage.output } } : {}),
      durationMs: Date.now() - started,
    });
  };

  // The history: the session's record, newest last (the assistant's view lists it so).
  const readHistory = async () => {
    const { items } = await viewer.record.read(200);
    history.set(historyOf(items.slice().reverse()));
  };
  const off = viewer.record.onChange(() => void readHistory());
  void readHistory();

  return {
    settings,
    turn: turn.asReadonly,
    view: view.asReadonly,
    history: history.asReadonly,
    ask(prompt) {
      const setup = settings.setup.value;
      if (!setup) {
        options.openSetup?.();
        return;
      }
      if (agent.session.busy) return;
      const started = Date.now();
      turn.set({
        id: `pending-${started}`,
        prompt,
        source: { kind: 'assistant' },
        status: 'thinking',
        text: '',
        runs: [],
        checkpoint: false,
        time: new Date(started),
      });
      view.set('card');
      void agent
        .turn(prompt, hooks, { provider: setup.provider, model: setup.model })
        .then((r) => settle(r, started))
        .catch((e) => update({ status: 'failed', text: e instanceof Error ? e.message : String(e) }));
    },
    stop: () => void agent.stop(),
    approve(runId, allow) {
      if (!pending.some((p) => p.callId === runId)) return;
      decided.set(runId, allow ? 'allow' : 'deny');
      if (!pending.every((p) => decided.has(p.callId))) return;
      const started = turn.value?.time.getTime() ?? Date.now();
      update({ status: 'working' });
      void agent.resume(Object.fromEntries(decided), hooks).then((r) => settle(r, started));
    },
    undo(turnId) {
      void person.history.undoTurn({ turn: turnId }).then((r) => {
        if (r.ok && turn.value?.id === turnId) update({ undone: true, checkpoint: false });
      });
    },
    dismiss: () => view.set('line'),
    hide: () => view.set('none'),
    reopen: () => view.set('card'),
    dispose() {
      off();
    },
  };
}

function statusOf(s: ToolEvent['status']): RunStatus {
  return s === 'started' || s === 'progress'
    ? 'running'
    : s === 'succeeded'
      ? 'done'
      : s === 'skipped'
        ? 'stopped'
        : 'failed';
}

/** The record's items as the assistant's history: the person's runs one by one, an agent's turn with its runs. */
export function historyOf(items: readonly RecordItem[]): HistoryItem[] {
  const out: HistoryItem[] = [];
  const turns = new Map<string, AgentTurn>();
  for (const { id, entry, open } of items) {
    if (entry.kind === 'turn') {
      const t: AgentTurn = {
        id: entry.turn,
        prompt: entry.prompt,
        source: { kind: 'assistant' },
        status: open
          ? 'working'
          : entry.outcome === 'done'
            ? 'done'
            : entry.outcome === 'stopped'
              ? 'stopped'
              : 'failed',
        text: entry.reply ?? '',
        runs: [],
        checkpoint: false,
        time: new Date(entry.at),
      };
      turns.set(entry.turn, t);
      out.push({ kind: 'turn', turn: t });
    } else if (entry.kind === 'run') {
      const run: CommandRun = {
        id,
        command: entry.operation,
        title: entry.result.ok ? entry.result.did : entry.title,
        source:
          entry.origin.kind === 'agent'
            ? { kind: 'assistant' }
            : { kind: entry.origin.kind === 'person' ? 'person' : 'script' },
        status: entry.result.ok ? 'done' : 'failed',
        ...(entry.result.ok ? {} : { detail: entry.result.message }),
        time: new Date(entry.at),
        durationMs: entry.ms,
      };
      const t = entry.origin.turn ? turns.get(entry.origin.turn) : undefined;
      if (t) (t.runs as CommandRun[]).push(run);
      else out.push({ kind: 'run', run });
    } else if (entry.kind === 'undo' && entry.turn) {
      const t = turns.get(entry.turn);
      if (t) (t as { undone?: boolean }).undone = true;
    }
  }
  return out;
}
