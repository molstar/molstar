/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ownsKey } from '@molstar/prism-ui/dom/events';
import { Icon, type IconName, ScrollArea, For, Show, onCleanup, StatusPill } from '@molstar/prism-ui';
import { distinctUntilChanged, filter, map, type Observable } from 'rxjs';
import type { AgentModel, AgentTurn, CommandRun, RunStatus, TurnStatus } from '../model.js';
import { classes } from './_generated/agent-turn.module.css.js';

import { derive } from '@molstar/prism-ui/reactive/state';
import { bind } from '@molstar/prism-ui/dom/scope';
import { keyedState } from '@molstar/prism-ui/reactive/keyed-state';
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
  // By id: a run's record is replaced on every change, and a row made anew would restart its spinner.
  const ids = runs.pipe(
    map((rs) => rs.map((r) => r.id)),
    distinctUntilChanged((a, b) => a.length === b.length && a.every((x, i) => x === b[i])),
  );
  // The runs by id, one map a change; a row is told only when its own run's record changes (a `find` a row, and a field,
  // on every streamed token was runs² a token).
  const byId = keyedState<string, CommandRun>();
  bind(runs, (rs) => {
    const m = new Map<string, CommandRun>();
    for (const r of rs) m.set(r.id, r);
    byId.replace(m);
  });
  return (
    <div class={props.class}>
      <For each={ids}>{(id) => <RunRow run={byId.key(id).pipe(filter((x): x is CommandRun => !!x))} />}</For>
    </div>
  );
}

/** Undo turn: back to the checkpoint before it (held only until the next action). */
export function UndoTurn(props: { agent: AgentModel; id: string }) {
  return (
    <button
      type="button"
      class={classes.undo}
      data-tip="Back to how it was before this request"
      onClick={() => props.agent.undo(props.id)}
    >
      <Icon name="undo-2" size={14} />
      Undo turn
    </button>
  );
}

/**
 * The agent at work, over the canvas at the bottom left (_spec/README.md, "The turn card"): the request, the model's
 * words as they stream, each command as it runs, Allow / Deny when one waits for the person, Stop; at the end, what it
 * did in one line, its trace (with the Inspector: `onTrace` opens it at the turn) and Undo turn (back to the checkpoint
 * before it).
 */
export function AgentTurnCard(props: { agent: AgentModel; onTrace?: (turnId: string) => void }) {
  const { agent } = props;
  // Where the card is: its viewer's keys only (several viewers on a page).
  const anchor = (<span hidden />) as HTMLElement;
  // Esc closes the card (to its pill): before the panels take it, but not from a field, a list, a menu or a dialog.
  const onKey = (e: KeyboardEvent) => {
    if (e.key !== 'Escape' || e.defaultPrevented || agent.view.value !== 'card' || !agent.turn.value) return;
    if (!ownsKey(anchor, e)) return;
    const target = e.target as HTMLElement | null;
    if (
      target?.isContentEditable ||
      target?.closest?.(
        'input, textarea, select, [role=listbox], [role=menu], [role=dialog], [data-msp-page], [popover]',
      )
    )
      return;
    e.preventDefault();
    e.stopPropagation();
    agent.dismiss();
  };
  document.addEventListener('keydown', onKey, true);
  onCleanup(() => document.removeEventListener('keydown', onKey, true));
  const which = (view: 'card' | 'line') => derive([agent.turn, agent.view], (t, v) => (v === view ? t?.id : undefined));
  return (
    <>
      {anchor}
      <Show when={which('line')}>{() => <TurnLine agent={agent} />}</Show>
      <Show when={which('card')}>{(id) => <Card agent={agent} id={id} onTrace={props.onTrace} />}</Show>
    </>
  );
}

/** How long a finished turn's line stays (the ring around its ✕ counts it down). */
const FinishedLineMs = 8000;

/**
 * A closed card: a line of the status column, like a task's: the request and its state (the spark pulsing while it
 * runs); a click opens the card again, ✕ takes the line away (the turn stays in History).
 */
function TurnLine(props: { agent: AgentModel }) {
  const turn = props.agent.turn.pipe(filter((t): t is AgentTurn => !!t));
  const status = turn.pipe(
    map((t) => t.status),
    distinctUntilChanged(),
  );
  return (
    <span class={classes.line} data-status={status}>
      <StatusPill
        icon={<Icon name="sparkles" size={14} class={classes.spark} />}
        text={turn.pipe(map((t) => t.prompt))}
        after={<span class={classes.state}>{status.pipe(map((s) => TurnLabels[s]))}</span>}
        label="Show the assistant’s turn"
        onClick={() => props.agent.reopen()}
        dismissLabel="Hide"
        onDismiss={() => props.agent.hide()}
        // Finished (done, stopped): it goes after a while; a turn that waits or failed stays.
        autoDismiss={status.pipe(map((s) => (s === 'done' || s === 'stopped' ? FinishedLineMs : undefined)))}
      />
    </span>
  );
}

/** The open card of a turn. */
function Card(props: { agent: AgentModel; id: string; onTrace?: (turnId: string) => void }) {
  const { agent, id } = props;
  return (
    <Show when={id}>
      {() => {
        const turn = agent.turn.pipe(
          filter((t): t is AgentTurn => t?.id === id),
          distinctUntilChanged(),
        );
        const status = turn.pipe(
          map((t) => t.status),
          distinctUntilChanged(),
        );
        const running = turn.pipe(map(isRunning), distinctUntilChanged());
        // The answer so far, and whether there is one (states: the paragraph is made once, its text follows).
        const text = derive(agent.turn, (t) => (t?.id === id ? t.text : ''));
        const hasText = derive(text, (x) => !!x);
        const waiting = turn.pipe(
          map((t) => t.runs.find((r) => r.status === 'waiting')),
          distinctUntilChanged(),
        );
        return (
          // Not a live region: the text streams in; the finished reply is said once by the shell's region
          // (accessibility.md; agent-model.ts announces it).
          <section class={classes.card} data-status={status} aria-label="Assistant">
            <header class={classes.header}>
              <Icon name="sparkles" size={14} class={classes.spark} />
              <span class={classes.prompt} data-tip-overflow="">
                {turn.pipe(map((t) => t.prompt))}
              </span>
              <span class={classes.state}>{status.pipe(map((s) => TurnLabels[s]))}</span>
              <Show when={running}>
                {() => (
                  <button type="button" class={classes.stop} onClick={() => agent.stop()}>
                    <Icon name="square" size={10} />
                    Stop
                  </button>
                )}
              </Show>
              {/* Closes to the pill (Esc too); a running turn goes on. */}
              <button
                type="button"
                class={classes.icon}
                aria-label="Close"
                data-tip="Close (Esc)"
                onClick={() => agent.dismiss()}
              >
                <Icon name="x" size={14} />
              </button>
            </header>
            <ScrollArea class={classes.body} contentClass={classes.bodyContent}>
              <Show when={hasText}>
                {() => (
                  // Made once: a streamed token changes its text, not the paragraph (made anew with the whole text a
                  // token, a long answer was its length² of work).
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
              {(run) => {
                let remember: HTMLInputElement | undefined;
                return (
                  <div class={classes.approval} role="group" aria-label="Approval">
                    <span class={classes.approvalText}>
                      Allow “{run.title}”? <span class={classes.approvalWhy}>{run.detail}</span>
                    </span>
                    <button type="button" class={classes.deny} onClick={() => agent.approve(run.id, false)}>
                      Deny
                    </button>
                    {run.domain && agent.allowDomain ? (
                      <>
                        <button type="button" class={classes.deny} onClick={() => agent.approve(run.id, true)}>
                          Allow once
                        </button>
                        {/* Its host from now on (agent.md §4): this session, or remembered in this browser. */}
                        <label class={classes.remember}>
                          <input type="checkbox" ref={(el: HTMLInputElement) => (remember = el)} />
                          Remember
                        </label>
                        <button
                          type="button"
                          class={classes.allow}
                          title={`Allow ${run.domain} from now on (Storage lists it, with Forget)`}
                          onClick={() => agent.allowDomain!(run.id, remember?.checked ? 'browser' : 'session')}
                        >
                          Allow {run.domain}
                        </button>
                      </>
                    ) : (
                      <button type="button" class={classes.allow} onClick={() => agent.approve(run.id, true)}>
                        Allow
                      </button>
                    )}
                  </div>
                );
              }}
            </Show>
            <Show when={running.pipe(map((r) => !r))}>
              {() => (
                <footer class={classes.footer}>
                  <span class={classes.summary}>{turn.pipe(map(turnSummary))}</span>
                  {props.onTrace ? (
                    <button
                      type="button"
                      class={classes.trace}
                      data-tip="Each step of this turn, in the Inspector"
                      onClick={() => props.onTrace?.(id)}
                    >
                      Trace
                    </button>
                  ) : null}
                  <Show
                    when={turn.pipe(map((t) => t.checkpoint && !t.undone))}
                    fallback={() => (
                      <Show when={turn.pipe(map((t) => !!t.undone))}>
                        {() => <span class={classes.undone}>Undone</span>}
                      </Show>
                    )}
                  >
                    {() => <UndoTurn agent={agent} id={id} />}
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
