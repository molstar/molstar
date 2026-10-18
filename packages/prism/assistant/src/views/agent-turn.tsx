/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ownsKey } from '@molstar/prism-ui/dom/events';
import { Icon, type IconName, ScrollArea, For, Show, StatusPill, type ReadonlyState } from '@molstar/prism-ui';
import { held, present } from '../held.js';
import type { AgentModel, AgentTurn, CommandRun, RunStatus, TurnStatus } from '../model.js';
import { classes } from './_generated/agent-turn.module.css.js';

import { derive } from '@molstar/prism-ui/reactive/state';
import { effect, watch } from '@molstar/prism-ui/dom/scope';
import { on } from '@molstar/prism-ui/dom/scoped';
import { setAttr, setShown, setStyle, setText } from '@molstar/prism-ui/dom/writes';
import { keyedState } from '@molstar/prism-ui/reactive/keyed-state';

export { held };

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
 * A read is quieter. Made once; one render writes what its record says that the DOM does not.
 */
export function RunRow(props: { run: ReadonlyState<CommandRun> }) {
  const run = props.run;
  // The mark is a region (a different icon per state); the rest of the row is written by the render below.
  const status = derive(run, (r) => r.status);
  const title = document.createTextNode('');
  const time = document.createTextNode('');
  const detailText = document.createTextNode('');
  const detail = (<span class={classes.runDetail}>{detailText}</span>) as HTMLElement;
  const fill = (<span class={classes.barFill} />) as HTMLElement;
  const bar = (<span class={classes.bar}>{fill}</span>) as HTMLElement;
  const root = (
    <div class={classes.run}>
      <Show when={status}>{(s) => <Icon name={StatusIcons[s]} size={14} class={classes.mark} />}</Show>
      <span class={classes.runTitle} data-tip-overflow="">
        {title}
      </span>
      <span class={classes.runTime}>{time}</span>
      {detail}
      {bar}
    </div>
  ) as HTMLElement;
  watch([run], () => {
    const r = run.value;
    setAttr(root, 'data-status', r.status);
    setAttr(root, 'data-read', !!r.read);
    setText(title, r.title);
    setText(time, formatDuration(r.durationMs));
    setText(detailText, r.detail ?? '');
    setShown(detail, !!r.detail);
    setShown(bar, r.progress !== undefined);
    setStyle(fill, 'width', `${Math.round((r.progress ?? 0) * 100)}%`);
  });
  return root;
}

/** A turn's runs, each row following its record. */
export function RunList(props: { turn: ReadonlyState<AgentTurn>; class?: string }) {
  const runs = derive(props.turn, (t) => t.runs);
  // By id: a run's record is replaced on every change, and a row made anew would restart its spinner.
  const ids = derive(
    runs,
    (rs) => rs.map((r) => r.id),
    (a, b) => a.length === b.length && a.every((x, i) => x === b[i]),
  );
  // The runs by id, one map a change; a row is told only when its own run's record changes (a `find` a row, and a field,
  // on every streamed token was runs² a token).
  const byId = keyedState<string, CommandRun>();
  effect(runs, (rs) => {
    const m = new Map<string, CommandRun>();
    for (const r of rs) m.set(r.id, r);
    byId.replace(m);
  });
  return (
    <div class={props.class}>
      <For each={ids}>{(id) => <RunRow run={held(byId.key(id), present)} />}</For>
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
  on(document, 'keydown', onKey, { capture: true });
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
  const turn = held(props.agent.turn, present);
  const status = derive(turn, (t) => t.status);
  const prompt = derive(turn, (t) => t.prompt);
  const promptText = document.createTextNode('');
  const stateText = document.createTextNode('');
  const root = (
    <span class={classes.line}>
      <StatusPill
        icon={<Icon name="sparkles" size={14} class={classes.spark} />}
        text={promptText}
        after={<span class={classes.state}>{stateText}</span>}
        label="Show the assistant’s turn"
        onClick={() => props.agent.reopen()}
        dismissLabel="Hide"
        onDismiss={() => props.agent.hide()}
        // Finished (done, stopped): it goes after a while; a turn that waits or failed stays.
        autoDismiss={derive(status, (s) => (s === 'done' || s === 'stopped' ? FinishedLineMs : undefined))}
      />
    </span>
  ) as HTMLElement;
  watch([status, prompt], () => {
    setAttr(root, 'data-status', status.value);
    setText(promptText, prompt.value);
    setText(stateText, TurnLabels[status.value]);
  });
  return root;
}

/** The answer so far: one paragraph, made once; a streamed token changes its text, not the paragraph. */
function StreamedText(props: { text: ReadonlyState<string>; running: ReadonlyState<boolean> }) {
  const { text, running } = props;
  const node = document.createTextNode('');
  const p = (<p class={classes.text}>{node}</p>) as HTMLElement;
  watch([text, running], () => {
    setText(node, text.value);
    setAttr(p, 'data-streaming', running.value);
  });
  return p;
}

/** What a finished turn did, in one line. */
function TurnSummaryText(props: { turn: ReadonlyState<AgentTurn> }) {
  const summary = derive(props.turn, turnSummary);
  const node = document.createTextNode('');
  watch([summary], () => setText(node, summary.value));
  return <span class={classes.summary}>{node}</span>;
}

/** The open card of a turn. */
function Card(props: { agent: AgentModel; id: string; onTrace?: (turnId: string) => void }) {
  const { agent, id } = props;
  return (
    <Show when={id}>
      {() => {
        const turn = held(agent.turn, (t): t is AgentTurn => t?.id === id);
        const status = derive(turn, (t) => t.status);
        const running = derive(turn, isRunning);
        // The answer so far, and whether there is one (states: the paragraph is made once, its text follows).
        const text = derive(agent.turn, (t) => (t?.id === id ? t.text : ''));
        const hasText = derive(text, (x) => !!x);
        const waiting = derive(turn, (t) => t.runs.find((r) => r.status === 'waiting'));
        const prompt = derive(turn, (t) => t.prompt);
        const promptText = document.createTextNode('');
        const stateText = document.createTextNode('');
        // Not a live region: the text streams in; the finished reply is said once by the shell's region
        // (accessibility.md; agent-model.ts announces it).
        const section = (
          <section class={classes.card} aria-label="Assistant">
            <header class={classes.header}>
              <Icon name="sparkles" size={14} class={classes.spark} />
              <span class={classes.prompt} data-tip-overflow="">
                {promptText}
              </span>
              <span class={classes.state}>{stateText}</span>
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
              {/* Made once: a streamed token changes its text, not the paragraph (made anew with the whole text a
                  token, a long answer was its length² of work). */}
              <Show when={hasText}>{() => <StreamedText text={text} running={running} />}</Show>
              <Show when={derive(turn, (t) => t.retries ?? 0)}>
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
            <Show when={derive(running, (r) => !r)}>
              {() => (
                <footer class={classes.footer}>
                  <TurnSummaryText turn={turn} />
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
                    when={derive(turn, (t) => t.checkpoint && !t.undone)}
                    fallback={() => (
                      <Show when={derive(turn, (t) => !!t.undone)}>
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
        ) as HTMLElement;
        watch([status, prompt], () => {
          setAttr(section, 'data-status', status.value);
          setText(promptText, prompt.value);
          setText(stateText, TurnLabels[status.value]);
        });
        return section;
      }}
    </Show>
  );
}
