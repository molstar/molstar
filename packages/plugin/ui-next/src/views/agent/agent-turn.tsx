/** @jsxImportSource @molstar/plugin-ui-next/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { distinctUntilChanged, filter, map, type Observable } from 'rxjs';
import { Icon, type IconName } from '../../components/icon/icon.js';
import { ScrollArea } from '../../components/scroll-area/scroll-area.js';
import { For, Show } from '../../dom/control-flow.js';
import type { AgentModel, AgentTurn, CommandRun, RunStatus, TurnStatus } from '../../model/agent.js';
import { classes } from './_generated/agent-turn.module.css.js';

const StatusIcons: Record<RunStatus, IconName> = {
  running: 'loader-circle',
  waiting: 'hourglass',
  done: 'check',
  failed: 'x',
  stopped: 'ban',
  denied: 'ban',
};

const TurnLabels: Record<TurnStatus, string> = {
  thinking: 'Thinking…',
  working: 'Working…',
  waiting: 'Waiting for you',
  done: 'Done',
  stopped: 'Stopped',
  failed: 'Failed',
};

export function formatDuration(ms: number | undefined) {
  if (ms === undefined) return '';
  return ms < 1000 ? `${Math.max(0.1, ms / 1000).toFixed(1)} s` : `${(ms / 1000).toFixed(1)} s`;
}

const tokens = (n: number) => (n >= 1000 ? `${(n / 1000).toFixed(1)}k` : `${n}`);

/** What a finished turn did, in one line: the commands it ran, how long, the tokens. */
export function turnSummary(t: AgentTurn) {
  const changes = t.runs.filter((r) => !r.read);
  const failed = t.runs.filter((r) => r.status === 'failed').length;
  const parts = [`${changes.length} command${changes.length === 1 ? '' : 's'}`];
  if (failed) parts.push(`${failed} failed`);
  if (t.durationMs !== undefined) parts.push(formatDuration(t.durationMs));
  if (t.usage) parts.push(`${tokens(t.usage.input)} in · ${tokens(t.usage.output)} out`);
  return parts.join(' · ');
}

export function isRunning(t: AgentTurn) {
  return t.status === 'thinking' || t.status === 'working' || t.status === 'waiting';
}

/**
 * One command run, following its record by id: its state (a spinner while running, waiting for approval, done,
 * failed, stopped), what it does in words, its outcome or error under it, and its progress while it reports one.
 * A read is quieter.
 */
export function RunRow(props: { run: Observable<CommandRun> }) {
  const run = props.run;
  const field = <K extends keyof CommandRun>(k: K) =>
    run.pipe(
      map((r) => r[k]),
      distinctUntilChanged(),
    );
  const status = field('status');
  return (
    <div class={classes.run} data-status={status} data-read={field('read').pipe(map((r) => (r ? '' : undefined)))}>
      <Show when={status}>{(s) => <Icon name={StatusIcons[s]} size={14} class={classes.mark} />}</Show>
      <span class={classes.runTitle} data-tip-overflow="">
        {field('title')}
      </span>
      <span class={classes.runTime}>{field('durationMs').pipe(map(formatDuration))}</span>
      <Show when={field('detail')}>{(d) => <span class={classes.runDetail}>{d}</span>}</Show>
      <Show when={field('progress').pipe(map((p) => p !== undefined))}>
        {() => (
          <span class={classes.bar}>
            <span
              class={classes.barFill}
              style={{ width: field('progress').pipe(map((p) => `${Math.round((p ?? 0) * 100)}%`)) }}
            />
          </span>
        )}
      </Show>
    </div>
  );
}

/** A turn's runs, each row following its record. */
export function RunList(props: { turn: Observable<AgentTurn>; class?: string }) {
  const runs = props.turn.pipe(map((t) => t.runs));
  return (
    <div class={props.class}>
      <For each={runs} key={(r) => r.id}>
        {(r) => (
          <RunRow
            run={runs.pipe(
              map((rs) => rs.find((x) => x.id === r.id)),
              filter((x): x is CommandRun => !!x),
              distinctUntilChanged(),
            )}
          />
        )}
      </For>
    </div>
  );
}

/**
 * The agent at work, over the canvas at the bottom left (design/agent.md §5): the request, the model's words as they
 * stream, each command as it runs, Allow / Deny when one waits for the person, Stop; at the end, what it did in one
 * line and Undo turn (back to the checkpoint before it).
 */
export function AgentTurnCard(props: { agent: AgentModel }) {
  const { agent } = props;
  return (
    <Show when={agent.turn.pipe(map((t) => t?.id))}>
      {(id) => {
        const turn = agent.turn.pipe(
          filter((t): t is AgentTurn => t?.id === id),
          distinctUntilChanged(),
        );
        const status = turn.pipe(
          map((t) => t.status),
          distinctUntilChanged(),
        );
        const running = turn.pipe(map(isRunning), distinctUntilChanged());
        const waiting = turn.pipe(
          map((t) => t.runs.find((r) => r.status === 'waiting')),
          distinctUntilChanged(),
        );
        return (
          <section class={classes.card} data-status={status} aria-label="Assistant" aria-live="polite">
            <header class={classes.header}>
              <Icon name="sparkles" size={14} class={classes.spark} />
              <span class={classes.prompt} data-tip-overflow="">
                {turn.pipe(map((t) => t.prompt))}
              </span>
              <span class={classes.state}>{status.pipe(map((s) => TurnLabels[s]))}</span>
              <Show
                when={running}
                fallback={() => (
                  <button
                    type="button"
                    class={classes.icon}
                    aria-label="Close"
                    data-tip="Close"
                    onClick={() => agent.dismiss()}
                  >
                    <Icon name="x" size={14} />
                  </button>
                )}
              >
                {() => (
                  <button type="button" class={classes.stop} onClick={() => agent.stop()}>
                    <Icon name="square" size={10} />
                    Stop
                  </button>
                )}
              </Show>
            </header>
            <ScrollArea class={classes.body} contentClass={classes.bodyContent}>
              <Show when={turn.pipe(map((t) => t.text))}>
                {(text) => (
                  <p class={classes.text} data-streaming={running.pipe(map((r) => (r ? '' : undefined)))}>
                    {text}
                  </p>
                )}
              </Show>
              <Show when={turn.pipe(map((t) => t.retries ?? 0))}>
                {(n) => (
                  <p class={classes.note}>
                    Retried the model {n === 1 ? 'once' : `${n} times`} (the provider was busy)
                  </p>
                )}
              </Show>
              <RunList turn={turn} class={classes.runs} />
            </ScrollArea>
            {/* A run waiting for the person: what it does, why it asks, and the answer. */}
            <Show when={waiting}>
              {(run) => (
                <div class={classes.approval} role="group" aria-label="Approval">
                  <span class={classes.approvalText}>
                    Allow “{run.title}”? <span class={classes.approvalWhy}>{run.detail}</span>
                  </span>
                  <button type="button" class={classes.deny} onClick={() => agent.approve(run.id, false)}>
                    Deny
                  </button>
                  <button type="button" class={classes.allow} onClick={() => agent.approve(run.id, true)}>
                    Allow
                  </button>
                </div>
              )}
            </Show>
            <Show when={running.pipe(map((r) => !r))}>
              {() => (
                <footer class={classes.footer}>
                  <span class={classes.summary}>{turn.pipe(map(turnSummary))}</span>
                  <Show
                    when={turn.pipe(map((t) => t.checkpoint && !t.undone))}
                    fallback={() => (
                      <Show when={turn.pipe(map((t) => !!t.undone))}>
                        {() => <span class={classes.summary}>Undone</span>}
                      </Show>
                    )}
                  >
                    {() => (
                      <button
                        type="button"
                        class={classes.undo}
                        data-tip="Back to how it was before this request"
                        onClick={() => agent.undo(id)}
                      >
                        <Icon name="undo-2" size={14} />
                        Undo turn
                      </button>
                    )}
                  </Show>
                </footer>
              )}
            </Show>
          </section>
        );
      }}
    </Show>
  );
}
