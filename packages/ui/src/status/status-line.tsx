/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { derive, type ReadonlyState } from '@molstar/core/reactive';
import { Icon } from '../components/icon/icon.js';
import {
  type Child,
  DomHost,
  effect,
  For,
  type MaybeReactive,
  onCleanup,
  restartable,
  Show,
  setAttr,
  setShown,
  setStyle,
  setText,
  watch,
} from '../dom/index.js';
import { LogAll, LogProblems } from '../log/log-model.js';
import { LogKindIcons } from '../log/log-view.js';
import { state } from '../state.js';
import { classes } from './_generated/status-line.module.css.js';
import type { StatusKind, StatusModel, StatusTask } from './status-model.js';

/** The status views' own text (English by default, `DefaultStatusLabels`); parametrized ones are functions. */
export interface StatusLabels {
  /** The message's tooltip when a click opens the log at its problems (an error). */
  showProblems: string;
  /** The message's tooltip when a click opens the whole log. */
  showLog: string;
  /** A pill's ✕. */
  dismiss: string;
  /** A task's Cancel, named after its task (two running tasks are told apart). */
  cancel: (task: string) => string;
  /** A chip's ✕: its accessible name, after the chip. */
  clear: (chip: string) => string;
  /** A chip's ✕: its tooltip. */
  clearTip: string;
}

export const DefaultStatusLabels: StatusLabels = {
  showProblems: 'Show the problems in the log',
  showLog: 'Show the log',
  dismiss: 'Dismiss',
  cancel: (task) => `Cancel ${task}`,
  clear: (chip) => `Clear ${chip}`,
  clearTip: 'Clear',
};

export interface StatusLineProps {
  status: StatusModel;
  /**
   * The message was clicked: the app opens its log (the log's filter is set first: the problems for an error, else
   * everything). Without it the message is not a button.
   */
  onOpenLog?(kind: StatusKind): void;
  /** Under the tasks (an agent's turn as a `StatusPill`). */
  after?: Child;
  /**
   * Which messages this line shows (status §2), when the app has several lines (one a pane): `id`, the messages of
   * that scope (`post(…, { scope })`); `rest`, the messages without one, the tasks, and the one live region (one line
   * in the app takes it). No scope: everything (one line in the app).
   */
  scope?: { readonly id?: string; readonly rest?: boolean };
  labels?: Partial<StatusLabels>;
  class?: string;
}

/**
 * The status line (status §2): the one polite live region (`spoken`), the message on show, the running tasks (a pill
 * each, made once per id), then `after`. A column growing upward from the bottom: the app places it (a stage's
 * bottom-left slot).
 */
export function StatusLine(props: StatusLineProps) {
  const { status } = props;
  const labels = { ...DefaultStatusLabels, ...props.labels };
  const openLog = props.onOpenLog
    ? (kind: StatusKind) => {
        status.log.setFilter(kind === 'error' ? LogProblems : LogAll);
        props.onOpenLog?.(kind);
      }
    : undefined;
  const scope = props.scope;
  // The line that speaks and shows the tasks: the only one, or the one taking the rest (one live region in the app).
  const home = !scope || !!scope.rest;
  const mine = scope
    ? derive(status.current, (m) =>
        m && (m.scope === undefined ? !!scope.rest : m.scope === scope.id) ? m : undefined,
      )
    : status.current;
  const spoken = document.createTextNode('');
  if (home) watch([status.spoken], () => setText(spoken, status.spoken.value));
  return (
    <div class={props.class ? `${classes.status} ${props.class}` : classes.status}>
      {/* One polite live region for results: the latest message only, not the tasks' progress. */}
      {home ? (
        <div class={classes.announce} role="status" aria-live="polite">
          {spoken}
        </div>
      ) : null}
      <Show when={mine}>
        {(m) => (
          <StatusPill
            kind={m.kind}
            icon={<Icon name={LogKindIcons[m.kind]} />}
            text={m.text}
            tip={openLog ? (m.kind === 'error' ? labels.showProblems : labels.showLog) : undefined}
            onClick={openLog ? () => openLog(m.kind) : undefined}
            // An error stays until its ✕ (the others fade, the model's timer).
            onDismiss={m.kind === 'error' ? () => status.dismiss(m.id) : undefined}
            dismissLabel={labels.dismiss}
          />
        )}
      </Show>
      {home ? <TaskPills status={status} labels={labels} /> : null}
      {props.after}
    </div>
  );
}

export interface StatusPillProps {
  icon: Child;
  text: MaybeReactive<string>;
  /** What follows the text (a state, a progress bar). */
  after?: Child;
  /** The message's kind (its color). */
  kind?: StatusKind;
  /** The clickable body's tooltip ("Show the log"). */
  tip?: string;
  /** The clickable body's accessible name, when its text is not. */
  label?: string;
  onClick?(): void;
  onDismiss?(): void;
  /** The ✕'s name and tooltip (default `DefaultStatusLabels.dismiss`). */
  dismissLabel?: string;
  /**
   * Dismissed on its own after this many ms (a finished turn), a ring around the ✕ drawing the time left; the pointer
   * over the line or the focus in it holds the time. Undefined: it stays.
   */
  autoDismiss?: MaybeReactive<number | undefined>;
}

/**
 * One line of the status column (status §2): an icon, its text, what follows, and a ✕; the body is a button when
 * `onClick` is given. With `autoDismiss`, it goes on its own (its countdown ring), held while hovered or focused.
 */
export function StatusPill(props: StatusPillProps) {
  // The time left, held while the pointer is over the line or the focus is in it. The timer is the scope's (dom §6):
  // paused with hidden UI, gone with the line.
  let left = 0;
  let since = 0;
  let held = false;
  let running = false;
  const ring = state<number | undefined>(undefined);
  const timer = restartable(() => {
    running = false;
    props.onDismiss?.();
  });
  const run = () => {
    timer.cancel();
    running = false;
    if (held || left <= 0) return;
    since = DomHost.now();
    running = true;
    timer.start(left);
  };
  const hold = (on: boolean) => {
    if (on === held) return;
    held = on;
    if (on) {
      if (running) left -= DomHost.now() - since;
      timer.cancel();
      running = false;
    } else run();
  };
  if (props.autoDismiss !== undefined && props.onDismiss) {
    effect(props.autoDismiss, (ms) => {
      left = ms ?? 0;
      ring.set(ms);
      run();
    });
    onCleanup(() => timer.cancel());
  }
  const body = (
    <>
      {props.icon}
      <span class={classes.text}>{props.text}</span>
      {props.after ?? null}
    </>
  );
  return (
    <span
      class={classes.pill}
      data-kind={props.kind}
      onPointerEnter={() => hold(true)}
      onPointerLeave={() => hold(false)}
      onFocusIn={() => hold(true)}
      onFocusOut={(e: FocusEvent) => {
        if (!(e.currentTarget as Element).contains(e.relatedTarget as Node | null)) hold(false);
      }}
    >
      {props.onClick ? (
        <button
          type="button"
          class={classes.main}
          data-tip={props.tip}
          aria-label={props.label}
          onClick={() => props.onClick?.()}
        >
          {body}
        </button>
      ) : (
        <span class={classes.main}>{body}</span>
      )}
      {props.onDismiss ? (
        <button
          type="button"
          class={classes.button}
          aria-label={props.dismissLabel ?? DefaultStatusLabels.dismiss}
          data-tip={props.dismissLabel ?? DefaultStatusLabels.dismiss}
          onClick={() => props.onDismiss?.()}
        >
          <Icon name="x" />
          <Show when={ring}>
            {(ms) => (
              // The time left, drawn down around the ✕ (held as the timer is, by the line's hover and focus).
              <svg
                class={classes.countdown}
                viewBox="0 0 20 20"
                aria-hidden="true"
                style={{ '--msui-countdown': `${ms}ms` }}
              >
                <circle cx="10" cy="10" r="9" pathLength="100" />
              </svg>
            )}
          </Show>
        </button>
      ) : null}
    </span>
  );
}

/**
 * The status line's tasks, a pill each following its task by id: a task's record is replaced on every progress event,
 * and a pill made anew would restart its spinner and replace its Cancel under a press (status §2).
 */
export function TaskPills(props: { status: StatusModel; labels?: Partial<StatusLabels> }) {
  const { status } = props;
  return <For each={status.taskIds}>{(id) => <TaskPill task={status.task(id)} labels={props.labels} />}</For>;
}

/**
 * One task: made once; one render writes its label, its progress and whether it can be cancelled (its Cancel named
 * after it: `StatusLabels.cancel`).
 */
export function TaskPill(props: { task: ReadonlyState<StatusTask | undefined>; labels?: Partial<StatusLabels> }) {
  const { task } = props;
  const cancelLabel = props.labels?.cancel ?? DefaultStatusLabels.cancel;
  const label = document.createTextNode('');
  const fill = (<span class={classes.fill} />) as HTMLElement;
  const bar = (
    <span class={classes.progress} role="progressbar" aria-valuemin="0" aria-valuemax="100">
      {fill}
    </span>
  ) as HTMLElement;
  // Hidden by its attribute, not removed: whether a task can be cancelled may change as it runs. The click reads the
  // task's latest `abort`.
  const cancel = (
    <button type="button" class={classes.button} onClick={() => task.value?.abort?.()}>
      <Icon name="x" />
    </button>
  ) as HTMLElement;
  const root = (
    <span class={classes.pill} data-task="">
      <span class={classes.main}>
        <Icon name="loader-circle" class={classes.spin} />
        <span class={classes.text}>{label}</span>
        {bar}
      </span>
      {cancel}
    </span>
  ) as HTMLElement;
  watch([task], () => {
    // Gone (the list drops its pill in the same flush): what it showed last stays.
    const t = task.value;
    if (!t) return;
    const percent = Math.round((t.progress ?? 0) * 100);
    setText(label, t.label);
    setShown(bar, t.progress !== undefined);
    setStyle(fill, 'width', `${percent}%`);
    setAttr(bar, 'aria-valuenow', t.progress === undefined ? undefined : String(percent));
    setAttr(bar, 'aria-label', t.label);
    setAttr(cancel, 'hidden', !t.abort);
    const name = cancelLabel(t.label);
    setAttr(cancel, 'aria-label', name);
    setAttr(cancel, 'data-tip', name);
  });
  return root;
}

/** What the person works on now (a selection, a focus) as chips (status §2): a click activates one, its ✕ clears it. */
export function StatusChips(props: {
  status: StatusModel;
  lead?: Child;
  labels?: Partial<StatusLabels>;
  class?: string;
}) {
  const clear = props.labels?.clear ?? DefaultStatusLabels.clear;
  const clearTip = props.labels?.clearTip ?? DefaultStatusLabels.clearTip;
  return (
    <div class={props.class ? `${classes.chips} ${props.class}` : classes.chips}>
      {props.lead ?? null}
      <For each={props.status.chips} by={(c) => c.id}>
        {(chip) => (
          <span class={classes.chip}>
            <button type="button" class={classes.chipMain} disabled={!chip.activate} onClick={() => chip.activate?.()}>
              <span class={classes.text}>{chip.label}</span>
            </button>
            {chip.dismiss ? (
              <button
                type="button"
                class={classes.button}
                aria-label={clear(chip.label)}
                data-tip={clearTip}
                onClick={() => chip.dismiss?.()}
              >
                <Icon name="x" />
              </button>
            ) : null}
          </span>
        )}
      </For>
    </div>
  );
}
