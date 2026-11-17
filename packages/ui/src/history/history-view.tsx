/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { derive } from '@molstar/core/reactive';
import { Button } from '../components/button/button.js';
import { Icon, type IconName } from '../components/icon/icon.js';
import { keepRows, VirtualList } from '../components/virtual-list/virtual-list.js';
import { Show } from '../dom/index.js';
import { ViewEmpty, ViewRow, ViewToolbar } from '../shell/view-frame.js';
import { classes } from './_generated/history-view.module.css.js';
import type { HistoryModel, HistoryStep, UndoResult } from './history-model.js';

/** An origin's icon and name in History ("You", "The assistant"); an unknown kind is shown by its name. */
export interface HistoryOriginLook {
  readonly icon: IconName;
  readonly title: string;
}

export const DefaultHistoryOrigins: Readonly<Record<string, HistoryOriginLook>> = {
  person: { icon: 'user', title: 'You' },
  agent: { icon: 'bot', title: 'The assistant' },
  script: { icon: 'file-code', title: 'A script' },
  url: { icon: 'download', title: 'A link' },
  link: { icon: 'link', title: 'A linked view' },
};

const time = new Intl.DateTimeFormat(undefined, { hour: '2-digit', minute: '2-digit', second: '2-digit' });

/** A row of History (oldest first): the start, a step done or undone, or one run of a step's group. */
export type HistoryRow =
  | { readonly kind: 'start'; readonly key: string; readonly current: boolean }
  | {
      readonly kind: 'step';
      readonly key: string;
      readonly step: HistoryStep;
      /** Undone: Redo puts it back (muted, after the position). */
      readonly undone: boolean;
      /** The position: the latest step done. */
      readonly current: boolean;
    }
  | { readonly kind: 'run'; readonly key: string; readonly label: string; readonly undone: boolean };

/**
 * History's rows (history §4): the start, the steps done oldest first (a group's runs under it), the position at the
 * latest, then the undone steps Redo would put back, the next first. Pure: a headless host lists the same.
 */
export function historyRows(steps: readonly HistoryStep[], redoable: readonly HistoryStep[]): HistoryRow[] {
  const rows: HistoryRow[] = [{ kind: 'start', key: 'start', current: steps.length === 0 }];
  const add = (step: HistoryStep, undone: boolean, current: boolean) => {
    rows.push({ kind: 'step', key: `s${step.id}`, step, undone, current });
    const runs = step.runs;
    if (runs.length > 1) {
      for (let i = 0; i < runs.length; i++) rows.push({ kind: 'run', key: `r${step.id}:${i}`, label: runs[i], undone });
    }
  };
  for (let i = steps.length - 1; i >= 0; i--) add(steps[i], false, i === 0);
  for (const step of redoable) add(step, true, false);
  return rows;
}

const sameRow = (a: HistoryRow, b: HistoryRow) => {
  if (a.kind !== b.kind) return false;
  if (a.kind === 'start') return a.current === (b as typeof a).current;
  if (a.kind === 'run') return a.label === (b as typeof a).label && a.undone === (b as typeof a).undone;
  const c = b as typeof a;
  return a.step === c.step && a.undone === c.undone && a.current === c.current;
};

/** ARIA's key syntax as a hint ("Control+Shift+Z" → "Ctrl+Shift+Z"). */
const keyHint = (keys: string | undefined) => (keys ? keys.replace('Control', 'Ctrl') : undefined);

export interface HistoryViewProps {
  history: HistoryModel<any>;
  /** The origins' icons and names (`DefaultHistoryOrigins`). */
  origins?: Readonly<Record<string, HistoryOriginLook>>;
  /** Undo's and Redo's keys, in ARIA's syntax (`Control+Z`, `Control+Shift+Z`): their tooltips say them. */
  keys?: { readonly undo?: string; readonly redo?: string };
  /** What an action gave back (the app says a failure: its status line). */
  onResult?(result: UndoResult): void;
  class?: string;
}

/**
 * History (history §4): Undo and Redo in a toolbar (their tooltips say the step and the key, or why not), the reason
 * when Undo is blocked, then the steps as a list of rows: oldest first, a group's runs under it, its origin and time as
 * the description, the position marked, "Undo to here" on a step before it and "Redo to here" on one undone (muted
 * after it), an unsafe step's warning. Every action is a model method (`undo`, `redo`, `moveTo`), so the same works
 * headless. It fills its parent (a flex column).
 */
export function HistoryView(props: HistoryViewProps) {
  const { history } = props;
  const origins = props.origins ?? DefaultHistoryOrigins;
  const undoKey = props.keys?.undo ?? 'Control+Z';
  const redoKey = props.keys?.redo ?? 'Control+Shift+Z';
  const keep = keepRows<HistoryRow>((r) => r.key, sameRow);
  const rows = derive([history.steps, history.redoable], (steps, redoable) => keep(historyRows(steps, redoable)));
  const busy = derive(history.running, (n) => n > 0);
  const undoTip = derive(
    [history.undoBlocked, history.steps],
    (reason, steps) => reason ?? `Undo ${steps[0]?.label ?? ''} (${keyHint(undoKey)})`,
  );
  const redoTip = derive([history.redoable, busy], (r, b) =>
    b ? 'A change is running' : r.length ? `Redo ${r[0].label} (${keyHint(redoKey)})` : 'Nothing to redo',
  );
  // The reason Undo cannot run, when there is something it would undo (an unsafe step, a run in flight).
  const blocked = derive([history.undoBlocked, history.steps], (reason, steps) =>
    steps.length && reason ? reason : undefined,
  );
  const count = derive(history.steps, (steps) =>
    steps.length ? `${steps.length} step${steps.length === 1 ? '' : 's'}` : '',
  );
  const act = (p: Promise<UndoResult>) => void p.then((r) => props.onResult?.(r));

  const row = (r: HistoryRow) => {
    if (r.kind === 'start') {
      return (
        <ViewRow
          class={classes.row}
          role="listitem"
          icon="history"
          label="Start"
          description="Before the first step"
          data-state={r.current ? 'current' : 'done'}
          aria-current={r.current ? 'step' : undefined}
          trailing={
            r.current ? undefined : (
              <Button
                class={classes.action}
                variant="ghost"
                size="xs"
                icon={<Icon name="undo-2" />}
                title="Undo every step"
                disabled={busy}
                onClick={() => act(history.moveTo(undefined))}
              />
            )
          }
        />
      );
    }
    if (r.kind === 'run') {
      return (
        <ViewRow
          class={classes.row}
          role="listitem"
          level={1}
          label={r.label}
          data-state={r.undone ? 'undone' : 'run'}
        />
      );
    }
    const { step } = r;
    const origin = origins[step.origin.kind];
    const trailing = (
      <>
        {step.unsafe ? (
          <span class={classes.unsafe} data-tip={step.unsafe} role="img" aria-label={step.unsafe}>
            <Icon name="triangle-alert" />
          </span>
        ) : null}
        {r.current ? null : (
          <Button
            class={classes.action}
            variant="ghost"
            size="xs"
            icon={<Icon name={r.undone ? 'redo-2' : 'undo-2'} />}
            title={r.undone ? 'Redo to here' : 'Undo to here'}
            disabled={busy}
            onClick={() => act(history.moveTo(step.id))}
          />
        )}
      </>
    );
    return (
      <ViewRow
        class={classes.row}
        role="listitem"
        icon={origin?.icon ?? 'history'}
        label={step.label}
        description={`${origin?.title ?? step.origin.kind} · ${time.format(step.time)}`}
        data-state={r.undone ? 'undone' : r.current ? 'current' : 'done'}
        data-unsafe={step.unsafe ? '' : undefined}
        aria-current={r.current ? 'step' : undefined}
        trailing={trailing}
      />
    );
  };

  return (
    <div class={props.class ? `${classes.history} ${props.class}` : classes.history}>
      <ViewToolbar label="History">
        <Button
          variant="ghost"
          size="sm"
          icon={<Icon name="undo-2" />}
          title={undoTip}
          aria-label="Undo"
          aria-description={history.undoBlocked}
          aria-keyshortcuts={undoKey}
          disabled={derive(history.canUndo, (can) => !can)}
          onClick={() => act(history.undo())}
        >
          Undo
        </Button>
        <Button
          variant="ghost"
          size="sm"
          icon={<Icon name="redo-2" />}
          title={redoTip}
          aria-label="Redo"
          aria-keyshortcuts={redoKey}
          disabled={derive([history.canRedo, busy], (can, b) => !can || b)}
          onClick={() => act(history.redo())}
        >
          Redo
        </Button>
        <span class={classes.count}>{count}</span>
      </ViewToolbar>
      <Show when={blocked}>
        {(reason) => (
          <p class={classes.blocked} role="note">
            <Icon name="info" />
            <span>{reason}</span>
          </p>
        )}
      </Show>
      <Show
        when={derive(rows, (list) => list.length > 1)}
        fallback={() => (
          <ViewEmpty>Nothing has happened yet. Each change shows here; Undo takes the latest back.</ViewEmpty>
        )}
      >
        {() => (
          <VirtualList items={rows} itemKey={(r) => r.key} class={classes.steps} edge role="list" aria-label="Steps">
            {row}
          </VirtualList>
        )}
      </Show>
    </div>
  );
}
