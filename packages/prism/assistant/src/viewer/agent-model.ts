/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { PendingApproval, ToolEvent, TurnHooks, TurnResult } from '@molstar/ai-session';
import { observableState, type ReadonlyState } from '@molstar/prism-ui/core';
import type { ViewerAgent } from '@molstar/prism-viewer/agent/session';
import { createHistoryModel, type HistoryRow, type HistoryRunRow } from '@molstar/prism-viewer/journal/history-model';
import type { PrismViewer } from '@molstar/prism-viewer/viewer';
import { map, Observable } from 'rxjs';
import type { AgentModel, AgentTurn, CommandRun, HistoryItem, RunStatus, TurnStatus } from '../model.js';
import type { ViewerAISettings } from './settings-model.js';

// The assistant's model over the viewer's agent (plan 4.4, agent.md §1): the latest turn from the session's hooks (its
// text, its calls, waiting for the person), asking with the model chosen in the setup, approvals answered as the turn
// asks them, Undo turn through the viewer's history (the person's API), and the history as the viewer's History reads
// the session's record (journal.md §3), read only while something watches it.
//
// One turn at a time: asking while one runs, or while one waits for an approval, is refused in words (the card shows
// the turn that is there). Stop ends the card's turn at once, parked or running. The card's checkpoint follows the
// viewer's history: held while the newest step is this turn's and can be restored.

export interface ViewerAgentModel extends AgentModel {
  dispose(): void;
}

const Finished: readonly TurnStatus[] = ['done', 'stopped', 'failed'];
const isFinished = (s: TurnStatus) => Finished.includes(s);
const messageOf = (e: unknown) => (e instanceof Error ? e.message : String(e));

export function createViewerAgentModel(
  viewer: PrismViewer,
  agent: ViewerAgent,
  settings: ViewerAISettings,
  options: { openSetup?: () => void } = {},
): ViewerAgentModel {
  const { notifications, history: steps } = viewer.managers;
  const turn = observableState<AgentTurn | undefined>(undefined);
  const view = observableState<'card' | 'line' | 'none'>('none');
  const person = viewer.api.as({ kind: 'person' });
  let pending: PendingApproval[] = [];
  const decided = new Map<string, 'allow' | 'deny'>();
  let disposed = false;
  /** The card's turn, counted: what a turn's hooks and results do is dropped once another turn is on the card. */
  let generation = 0;

  /** Undo turn can restore it: the newest step of the history is this turn's (one group per turn) and safe. */
  const holdsCheckpoint = (t: AgentTurn) => {
    const top = steps.steps.value[0];
    return (
      isFinished(t.status) &&
      !t.undone &&
      !!top &&
      top.origin.kind === 'agent' &&
      top.origin.turn === t.id &&
      !top.unsafe
    );
  };
  const update = (patch: Partial<AgentTurn>, of = generation) => {
    const t = turn.value;
    if (!t || disposed || of !== generation) return;
    // A finished card does not go back to running (a hook arriving after Stop).
    const next = {
      ...t,
      ...patch,
      ...(isFinished(t.status) && patch.status && !isFinished(patch.status) ? { status: t.status } : {}),
    };
    turn.set({ ...next, checkpoint: holdsCheckpoint(next) });
  };
  const fail = (of: number) => (e: unknown) => {
    pending = [];
    decided.clear();
    update({ status: 'failed', text: messageOf(e) }, of);
  };
  // Released by the next step (the person's, another turn's), restored by an undo.
  const stepsSub = steps.steps.subscribe(() => update({}));

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
  const hooksFor = (of: number): TurnHooks => ({
    onTurn: (e) => {
      if (e.kind === 'begun') update({ id: e.turnId }, of);
    },
    onText: (text) => update({ text }, of),
    onTool: (e) => {
      const t = turn.value;
      if (!t || of !== generation) return;
      const run = runOf(e);
      const at = t.runs.findIndex((r) => r.id === run.id);
      update({ runs: at < 0 ? [...t.runs, run] : t.runs.map((r, i) => (i === at ? run : r)), status: 'working' }, of);
    },
    onStatus: (e) => {
      if (e.kind === 'model') update({ status: 'thinking' }, of);
      else if (e.kind === 'tools') update({ status: 'working' }, of);
      // `permission` opens for every call with a rule, the policy's yes too: waiting is the turn's `awaiting-approval`.
      else if (e.kind === 'retrying') update({ retries: (turn.value?.retries ?? 0) + 1 }, of);
    },
  });
  const settle = (r: TurnResult | undefined, started: number, of: number) => {
    const t = turn.value;
    if (!t || of !== generation) return;
    if (!r) {
      // Nothing was parked any more (stopped meanwhile): the agent has ended the turn.
      pending = [];
      decided.clear();
      update({ status: 'stopped', durationMs: Date.now() - started }, of);
      return;
    }
    if (r.status === 'awaiting-approval') {
      pending = r.approval ?? [];
      decided.clear();
      notifications.announce('The assistant is waiting for your approval.');
      const waiting = new Set(pending.map((p) => p.callId));
      update(
        {
          status: 'waiting',
          runs: [
            ...t.runs,
            ...pending
              .filter((p) => !t.runs.some((x) => x.id === p.callId))
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
        },
        of,
      );
      return;
    }
    pending = [];
    decided.clear();
    const status = r.status === 'completed' ? 'done' : r.status === 'stopped' ? 'stopped' : 'failed';
    const text = r.text || t.text || (r.error ?? '');
    // The finished reply, said once (the card's text streams and is not live; README §3.7).
    notifications.announce(
      status === 'done' ? `Assistant: ${text || 'done'}` : `Assistant ${status}${text ? `: ${text}` : ''}`,
    );
    update(
      {
        status,
        text,
        ...(r.usage ? { usage: { input: r.usage.input, output: r.usage.output } } : {}),
        durationMs: Date.now() - started,
      },
      of,
    );
  };

  /** Asking is refused while a turn is there: the person sees it (its card) and why. */
  const refuse = (text: string) => {
    notifications.post(text, 'info');
    view.set('card');
  };

  // ---- The history: the viewer's History rows as the assistant's items, read only while watched ----
  const records = createHistoryModel(viewer.record);
  let watchers = 0;
  let mapped: { rows: readonly HistoryRow[]; items: readonly HistoryItem[] } | undefined;
  const itemsOf = (rows: readonly HistoryRow[]) => {
    if (mapped?.rows !== rows) mapped = { rows, items: historyItemsOf(rows) };
    return mapped.items;
  };
  const history = new Observable<readonly HistoryItem[]>((subscriber) => {
    if (disposed) {
      subscriber.complete();
      return;
    }
    if (watchers++ === 0) records.setActive(true);
    const sub = records.rows.pipe(map(itemsOf)).subscribe(subscriber);
    return () => {
      sub.unsubscribe();
      if (--watchers === 0) records.setActive(false);
    };
  }) as ReadonlyState<readonly HistoryItem[]>;
  Object.defineProperty(history, 'value', { get: () => itemsOf(records.rows.value) });

  return {
    settings,
    turn: turn.asReadonly,
    view: view.asReadonly,
    history,
    ask(prompt) {
      if (disposed) return;
      const setup = settings.setup.value;
      if (!setup) {
        options.openSetup?.();
        return;
      }
      const t = turn.value;
      if (agent.parked || t?.status === 'waiting') return refuse('Answer or stop the waiting step first.');
      if (agent.session.busy || (t && !isFinished(t.status)))
        return refuse('A turn is running: stop it or wait for it.');
      const of = ++generation;
      const started = Date.now();
      pending = [];
      decided.clear();
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
      try {
        agent
          .turn(prompt, hooksFor(of), { provider: setup.provider, model: setup.model })
          .then((r) => settle(r, started, of))
          .catch(fail(of));
      } catch (e) {
        fail(of)(e);
      }
    },
    stop() {
      const t = turn.value;
      const live = !!t && !isFinished(t.status);
      if (t && live) {
        pending = [];
        decided.clear();
        update({
          status: 'stopped',
          runs: t.runs.map((r) => (r.status === 'waiting' || r.status === 'running' ? { ...r, status: 'stopped' } : r)),
          durationMs: Date.now() - t.time.getTime(),
        });
      }
      // The agent ends a parked turn too (its record entry, its history group).
      if (live || agent.session.busy || agent.parked) {
        agent.stop().catch((e) => notifications.post(`The turn could not be stopped: ${messageOf(e)}`, 'warning'));
      }
    },
    approve(runId, allow) {
      const t = turn.value;
      if (!t || t.status !== 'waiting' || decided.has(runId) || !pending.some((p) => p.callId === runId)) return;
      decided.set(runId, allow ? 'allow' : 'deny');
      // Answered at once on the card; the turn goes on once every waiting call is answered.
      update({
        runs: t.runs.map((r) =>
          r.id === runId
            ? { ...r, status: allow ? 'running' : 'denied', detail: allow ? 'Allowed.' : 'Not allowed.' }
            : r,
        ),
      });
      if (!pending.every((p) => decided.has(p.callId))) return;
      const of = generation;
      const started = t.time.getTime();
      const decisions = Object.fromEntries(decided);
      pending = [];
      update({ status: 'working' });
      try {
        agent
          .resume(decisions, hooksFor(of))
          .then((r) => settle(r, started, of))
          .catch(fail(of));
      } catch (e) {
        fail(of)(e);
      }
    },
    undo(turnId) {
      person.history.undoTurn({ turn: turnId }).then(
        (r) => {
          if (disposed) return;
          if (!r.ok) notifications.post(r.message, 'warning');
          else if (turn.value?.id === turnId) update({ undone: true });
        },
        (e) => notifications.post(`The turn could not be undone: ${messageOf(e)}`, 'warning'),
      );
    },
    dismiss: () => view.set('line'),
    hide: () => view.set('none'),
    reopen: () => view.set('card'),
    dispose() {
      if (disposed) return;
      const t = turn.value;
      // A turn running or parked ends with the model: nothing is left waiting for a card that is gone.
      if (agent.session.busy || agent.parked || (t && !isFinished(t.status))) {
        agent.stop().catch((e) => console.error('[prism-assistant] stop:', e));
      }
      disposed = true;
      stepsSub.unsubscribe();
      records.dispose();
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

const runItemOf = (r: HistoryRunRow): CommandRun => ({
  id: r.key,
  command: r.operation,
  title: r.failed ? r.title : r.text,
  source:
    r.origin === 'agent' ? { kind: 'assistant' } : r.origin === 'person' ? { kind: 'person' } : { kind: 'script' },
  status: r.failed ? 'failed' : 'done',
  ...(r.failed ? { detail: r.text } : r.undone ? { detail: 'Undone' } : {}),
  time: new Date(r.at),
});

/**
 * History's rows (`rowsOf`, newest first) as the assistant's items (newest last): the person's runs one by one, an
 * agent's turn with its runs, undone as History marks it. Session events are History's separators; the items have none.
 */
export function historyItemsOf(rows: readonly HistoryRow[]): HistoryItem[] {
  const out: HistoryItem[] = [];
  for (let i = rows.length - 1; i >= 0; i--) {
    const row = rows[i];
    if (row.kind === 'run') out.push({ kind: 'run', run: runItemOf(row) });
    else if (row.kind === 'turn') {
      out.push({
        kind: 'turn',
        turn: {
          id: row.turn,
          prompt: row.prompt,
          source: { kind: 'assistant' },
          status: row.open
            ? 'working'
            : row.outcome === 'done'
              ? 'done'
              : row.outcome === 'stopped'
                ? 'stopped'
                : 'failed',
          text: row.reply ?? '',
          runs: row.runs.map(runItemOf),
          checkpoint: false,
          ...(row.undone ? { undone: true } : {}),
          time: new Date(row.at),
        },
      });
    }
  }
  return out;
}
