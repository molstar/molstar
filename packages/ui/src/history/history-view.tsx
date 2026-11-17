/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { derive, type ReadonlyState, State } from '@molstar/core/reactive';
import { ariaKeyShortcut, formatChord } from '../commands/keys.js';
import type { CommandManager } from '../commands/manager.js';
import type { KeyChord } from '../commands/types.js';
import { Button } from '../components/button/button.js';
import { Icon, type IconName } from '../components/icon/icon.js';
import { keepRows, VirtualList } from '../components/virtual-list/virtual-list.js';
import { Show } from '../dom/index.js';
import { ViewEmpty, ViewRow, ViewToolbar } from '../shell/view-frame.js';
import { classes } from './_generated/history-view.module.css.js';
import { HistoryCommandIds, HistoryShortcuts } from './history-commands.js';
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

/** History's words (history §6); the refusals are the model's (`HistoryModel.labels`). */
export interface HistoryViewLabels {
  /** The toolbar's accessible name. */
  readonly toolbar: string;
  readonly undo: string;
  readonly redo: string;
  /** Undo's tooltip: the step and its key (none when it has no key). */
  readonly undoTip: (step: string, key: string | undefined) => string;
  readonly redoTip: (step: string, key: string | undefined) => string;
  /** A row's action before the position. */
  readonly undoTo: (step: string) => string;
  /** A row's action after it (an undone step). */
  readonly redoTo: (step: string) => string;
  readonly start: string;
  readonly startDescription: string;
  /** The start's action. */
  readonly undoAll: string;
  readonly count: (steps: number) => string;
  /** The list's accessible name. */
  readonly list: string;
  readonly empty: string;
  /** A group's runs above those shown (folded, or no longer kept). */
  readonly earlierRuns: (count: number) => string;
  /** A group's every kept run shown. */
  readonly allRuns: (count: number) => string;
  /** The toggle (pressed: every kept run listed). */
  readonly showRuns: string;
}

export const DefaultHistoryViewLabels: HistoryViewLabels = {
  toolbar: 'History',
  undo: 'Undo',
  redo: 'Redo',
  undoTip: (step, key) => `Undo ${step}${key ? ` (${key})` : ''}`,
  redoTip: (step, key) => `Redo ${step}${key ? ` (${key})` : ''}`,
  undoTo: (step) => `Undo to ${step}`,
  redoTo: (step) => `Redo to ${step}`,
  start: 'Start',
  startDescription: 'Before the first step',
  undoAll: 'Undo every step',
  count: (n) => (n ? `${n} step${n === 1 ? '' : 's'}` : ''),
  list: 'Steps',
  empty: 'Nothing has happened yet. Each change shows here; Undo takes the latest back.',
  earlierRuns: (n) => `${n} earlier run${n === 1 ? '' : 's'}`,
  allRuns: (n) => `All ${n} runs`,
  showRuns: 'Show all',
};

const time = new Intl.DateTimeFormat(undefined, { hour: '2-digit', minute: '2-digit', second: '2-digit' });

/** A row of History (oldest first): the start, a step done or undone, a group's earlier runs, or one of its runs. */
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
  | {
      readonly kind: 'more';
      readonly key: string;
      readonly step: HistoryStep;
      /** Runs above those shown: folded, or past the model's `maxRuns`. */
      readonly earlier: number;
      /** Every kept run is listed. */
      readonly expanded: boolean;
      /** More runs are kept than a folded step shows: it can open (and close). */
      readonly expandable: boolean;
      readonly undone: boolean;
    }
  | { readonly kind: 'run'; readonly key: string; readonly label: string; readonly undone: boolean };

export interface HistoryRowsOptions {
  /** Runs listed under a folded step (10): its latest; a row says how many are above. */
  readonly runs?: number;
  /** The steps whose every kept run is listed (at most the model's `maxRuns`). */
  readonly expanded?: ReadonlySet<number>;
}

/**
 * History's rows (history §4): the start, the steps done oldest first (a group's latest runs under it, a row for the
 * earlier ones), the position at the latest, then the undone steps Redo would put back, the next first. Each step costs
 * at most `runs` rows unless expanded, and reads only those runs (`lastRuns`): a long turn published per run stays
 * linear while History is open. Pure: a headless host lists the same.
 */
export function historyRows(
  steps: readonly HistoryStep[],
  redoable: readonly HistoryStep[],
  options: HistoryRowsOptions = {},
): HistoryRow[] {
  const cap = Math.max(1, options.runs ?? 10);
  const rows: HistoryRow[] = [{ kind: 'start', key: 'start', current: steps.length === 0 }];
  const add = (step: HistoryStep, undone: boolean, current: boolean) => {
    rows.push({ kind: 'step', key: `s${step.id}`, step, undone, current });
    const count = step.runCount;
    if (count <= 1) return;
    const kept = count - step.droppedRuns;
    const expanded = kept > cap && !!options.expanded?.has(step.id);
    const shown = step.lastRuns(expanded ? kept : cap);
    const first = count - shown.length;
    if (first > 0 || expanded) {
      rows.push({ kind: 'more', key: `m${step.id}`, step, earlier: first, expanded, expandable: kept > cap, undone });
    }
    for (let i = 0; i < shown.length; i++) {
      rows.push({ kind: 'run', key: `r${step.id}:${first + i}`, label: shown[i], undone });
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
  if (a.kind === 'more') {
    const c = b as typeof a;
    return (
      a.step.id === c.step.id &&
      a.earlier === c.earlier &&
      a.expanded === c.expanded &&
      a.expandable === c.expandable &&
      a.undone === c.undone
    );
  }
  const c = b as typeof a;
  return a.step === c.step && a.undone === c.undone && a.current === c.current;
};

/** ARIA's key syntax as a hint ("Control+Shift+Z" → "Ctrl+Shift+Z"). */
const keyHint = (keys: string) => keys.replace('Control', 'Ctrl');

interface KeyLook {
  /** `aria-keyshortcuts`. */
  readonly aria?: string;
  /** As a tooltip says it (⌘Z, Ctrl+Z). */
  readonly hint?: string;
}

const lookOf = (chord: KeyChord | undefined): KeyLook =>
  chord ? { aria: ariaKeyShortcut(chord), hint: formatChord(chord) } : {};

/**
 * Undo's or Redo's key (history §4): given in ARIA's syntax; else the command's key now from the app's command
 * manager (after its keymap: a remap shows); else its own key for this platform (`KeyPlatform`).
 */
function keyOf(
  which: 'undo' | 'redo',
  given: string | undefined,
  commands: CommandManager<any> | undefined,
): ReadonlyState<KeyLook> {
  if (given) return new State<KeyLook>({ aria: given, hint: keyHint(given) });
  if (!commands) return new State(lookOf(HistoryShortcuts[which]));
  return derive([commands.keymap.layers, commands.commands], () => {
    const c = commands.registered(HistoryCommandIds[which]);
    return lookOf(c ? commands.shortcut(c) : HistoryShortcuts[which]);
  });
}

export interface HistoryViewProps {
  history: HistoryModel<any>;
  /** The origins' icons and names (`DefaultHistoryOrigins`). */
  origins?: Readonly<Record<string, HistoryOriginLook>>;
  /** The command manager `historyCommands` are registered with: Undo's and Redo's keys are read from it. */
  commands?: CommandManager<any>;
  /** Undo's and Redo's keys, in ARIA's syntax (`Meta+Z`), over the commands'; else ⌘Z / Ctrl+Z by platform. */
  keys?: { readonly undo?: string; readonly redo?: string };
  /** Runs listed under a group's step before "Show all" (10). */
  runs?: number;
  labels?: Partial<HistoryViewLabels>;
  /** What an action gave back (the app says a failure: its status line). */
  onResult?(result: UndoResult): void;
  class?: string;
}

/**
 * History (history §4): Undo and Redo in a toolbar (their tooltips say the step and the key, or why not), the reason
 * when Undo is blocked (not while a change runs: the buttons are greyed), then the steps as a list of rows: oldest
 * first, a group's latest runs under it, its origin and time as the description, the position marked, "Undo to <step>"
 * on a step before it and "Redo to <step>" on one undone (muted after it), an unsafe step's warning. Every action is a
 * model method (`undo`, `redo`, `moveTo`), so the same works headless. It fills its parent (a flex column).
 */
export function HistoryView(props: HistoryViewProps) {
  const { history } = props;
  const l = { ...DefaultHistoryViewLabels, ...props.labels };
  const origins = props.origins ?? DefaultHistoryOrigins;
  const undoKey = keyOf('undo', props.keys?.undo, props.commands);
  const redoKey = keyOf('redo', props.keys?.redo, props.commands);
  const expanded = new State<ReadonlySet<number>>(new Set());
  const toggle = (id: number) => {
    const next = new Set(expanded.value);
    if (!next.delete(id)) next.add(id);
    expanded.set(next);
  };
  const keep = keepRows<HistoryRow>((r) => r.key, sameRow);
  const rows = derive([history.steps, history.redoable, expanded], (steps, redoable, open) =>
    keep(historyRows(steps, redoable, { runs: props.runs, expanded: open })),
  );
  const busy = derive(history.running, (n) => n > 0);
  const undoTip = derive(
    [history.undoBlocked, history.steps, undoKey],
    (reason, steps, key) => reason ?? l.undoTip(steps[0]?.label ?? '', key.hint),
  );
  const redoTip = derive([history.redoable, busy, redoKey], (r, b, key) =>
    b ? history.labels.running : r.length ? l.redoTip(r[0].label, key.hint) : history.labels.nothingToRedo,
  );
  // The reason Undo cannot run, when there is something it would undo (an unsafe step), but not a run in flight: that
  // is transient (the buttons are greyed meanwhile), and a note on every run would flash.
  const blocked = derive([history.undoBlocked, history.undoBlockedCode, history.steps], (reason, code, steps) =>
    steps.length && reason && code !== 'running' ? reason : undefined,
  );
  const count = derive(history.steps, (steps) => l.count(steps.length));
  const act = (p: Promise<UndoResult>) => void p.then((r) => props.onResult?.(r));

  const row = (r: HistoryRow) => {
    if (r.kind === 'start') {
      return (
        <ViewRow
          class={classes.row}
          role="listitem"
          icon="history"
          label={l.start}
          description={l.startDescription}
          data-state={r.current ? 'current' : 'done'}
          aria-current={r.current ? 'step' : undefined}
          trailing={
            r.current ? undefined : (
              <Button
                class={classes.action}
                variant="ghost"
                size="xs"
                icon={<Icon name="undo-2" />}
                title={l.undoAll}
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
    if (r.kind === 'more') {
      return (
        <ViewRow
          class={classes.row}
          role="listitem"
          level={1}
          label={r.earlier ? l.earlierRuns(r.earlier) : l.allRuns(r.step.runCount)}
          data-state={r.undone ? 'undone' : 'run'}
          data-more=""
          trailing={
            r.expandable ? (
              <Button variant="ghost" size="xs" pressed={r.expanded} onClick={() => toggle(r.step.id)}>
                {l.showRuns}
              </Button>
            ) : undefined
          }
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
            title={r.undone ? l.redoTo(step.label) : l.undoTo(step.label)}
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
      <ViewToolbar label={l.toolbar}>
        <Button
          variant="ghost"
          size="sm"
          icon={<Icon name="undo-2" />}
          title={undoTip}
          aria-label={l.undo}
          aria-description={history.undoBlocked}
          aria-keyshortcuts={derive(undoKey, (k) => k.aria)}
          disabled={derive(history.canUndo, (can) => !can)}
          onClick={() => act(history.undo())}
        >
          {l.undo}
        </Button>
        <Button
          variant="ghost"
          size="sm"
          icon={<Icon name="redo-2" />}
          title={redoTip}
          aria-label={l.redo}
          aria-keyshortcuts={derive(redoKey, (k) => k.aria)}
          disabled={derive([history.canRedo, busy], (can, b) => !can || b)}
          onClick={() => act(history.redo())}
        >
          {l.redo}
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
      <Show when={derive(rows, (list) => list.length > 1)} fallback={() => <ViewEmpty>{l.empty}</ViewEmpty>}>
        {() => (
          <VirtualList items={rows} itemKey={(r) => r.key} class={classes.steps} edge role="list" aria-label={l.list}>
            {row}
          </VirtualList>
        )}
      </Show>
    </div>
  );
}
