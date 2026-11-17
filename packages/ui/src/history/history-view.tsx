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
import { VirtualList } from '../components/virtual-list/virtual-list.js';
import { Show } from '../dom/index.js';
import { ViewEmpty, ViewRow, ViewToolbar } from '../shell/view-frame.js';
import { classes } from './_generated/history-view.module.css.js';
import { HistoryCommandIds, HistoryShortcuts } from './history-commands.js';
import type { HistoryModel, HistoryStep, UndoResult } from './history-model.js';

/** An origin's icon in History, by its kind (an unknown kind: a history icon). Its name is the labels' (`origin`). */
export const DefaultHistoryOriginIcons: Readonly<Record<string, IconName>> = {
  person: 'user',
  agent: 'bot',
  script: 'file-code',
  url: 'download',
  link: 'link',
};

const DefaultOriginNames: Readonly<Record<string, string>> = {
  person: 'You',
  agent: 'The assistant',
  script: 'A script',
  url: 'A link',
  link: 'A linked view',
};

const time = new Intl.DateTimeFormat(undefined, { hour: '2-digit', minute: '2-digit', second: '2-digit' });

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
  /** The toggle's text (pressed: every kept run listed). */
  readonly showRuns: string;
  /** The toggle's accessible name, by its step (one per group: "Show all" alone told them apart by nothing). */
  readonly showRunsName: (step: string) => string;
  /** An origin's name by its kind ("You", "The assistant"; an unknown kind: itself). */
  readonly origin: (kind: string) => string;
  /** A step's description: its origin's name and its time (ms, the model's `now`). */
  readonly description: (origin: string, time: number) => string;
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
  showRunsName: (step) => `Show all runs of ${step}`,
  origin: (kind) => DefaultOriginNames[kind] ?? kind,
  description: (origin, at) => `${origin} · ${time.format(at)}`,
};

/** A row of History (oldest first): the start, a step done or undone, a group's earlier runs, or one of its runs. */
export type HistoryRow =
  | {
      readonly kind: 'start';
      readonly key: string;
      readonly current: boolean;
      /** "Undo every step" can run: no step is unsafe. */
      readonly reachable: boolean;
    }
  | {
      readonly kind: 'step';
      readonly key: string;
      readonly step: HistoryStep;
      /** Undone: Redo puts it back (muted, after the position). */
      readonly undone: boolean;
      /** The position: the latest step done. */
      readonly current: boolean;
      /** Its move can run: no later step on the way is unsafe (`undoTo` is all or none, history §2). */
      readonly reachable: boolean;
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

/** History's rows now: the steps (newest first), the redoable ones (the next first), the steps whose runs all show. */
export type HistoryRows = (
  steps: readonly HistoryStep[],
  redoable: readonly HistoryStep[],
  expanded?: ReadonlySet<number>,
) => readonly HistoryRow[];

const NoneExpanded: ReadonlySet<number> = new Set();

interface CachedStep {
  readonly step: HistoryStep;
  readonly undone: boolean;
  readonly current: boolean;
  readonly reachable: boolean;
  readonly expanded: boolean;
  readonly rows: readonly HistoryRow[];
}

/**
 * History's rows (history §4), kept between calls: the start, the steps done oldest first (a group's latest `runs` (10)
 * runs under it, a row for the earlier ones), the position at the latest, then the undone steps Redo would put back, the
 * next first. A step whose object and place are unchanged keeps its rows as they were, and a changed one keeps those of
 * its rows that are the same: a long turn published per run rebuilds its own step's rows and the list's order, not
 * every row (a row is `lastRuns`' at most `runs`, never the whole list). DOM-free: a headless host lists the same.
 */
export function createHistoryRows(options: { readonly runs?: number } = {}): HistoryRows {
  const cap = Math.max(1, options.runs ?? 10);
  let cache = new Map<number, CachedStep>();
  let start: Extract<HistoryRow, { kind: 'start' }> | undefined;
  return (steps, redoable, expanded = NoneExpanded) => {
    // A move to a step undoes every later one, all or none: reachable while none of those is unsafe (history §2).
    let firstUnsafe = steps.length;
    for (let i = 0; i < steps.length; i++) {
      if (steps[i].unsafe) {
        firstUnsafe = i;
        break;
      }
    }
    const current = steps.length === 0;
    const all = firstUnsafe === steps.length;
    if (!start || start.current !== current || start.reachable !== all) {
      start = { kind: 'start', key: 'start', current, reachable: all };
    }
    const rows: HistoryRow[] = [start];
    const next = new Map<number, CachedStep>();
    const add = (step: HistoryStep, undone: boolean, current: boolean, reachable: boolean) => {
      const open = expanded.has(step.id);
      let c = cache.get(step.id);
      if (
        !c ||
        c.step !== step ||
        c.undone !== undone ||
        c.current !== current ||
        c.reachable !== reachable ||
        c.expanded !== open
      ) {
        c = {
          step,
          undone,
          current,
          reachable,
          expanded: open,
          rows: stepRows(step, undone, current, reachable, open, cap, c?.rows),
        };
      }
      next.set(step.id, c);
      for (const r of c.rows) rows.push(r);
    };
    for (let i = steps.length - 1; i >= 0; i--) add(steps[i], false, i === 0, i <= firstUnsafe);
    for (const step of redoable) add(step, true, false, true);
    cache = next;
    return rows;
  };
}

/** A step's rows, those the same as before (`previous`, by key) kept as they were. */
function stepRows(
  step: HistoryStep,
  undone: boolean,
  current: boolean,
  reachable: boolean,
  open: boolean,
  cap: number,
  previous: readonly HistoryRow[] | undefined,
): HistoryRow[] {
  const made: HistoryRow[] = [{ kind: 'step', key: `s${step.id}`, step, undone, current, reachable }];
  const count = step.runCount;
  if (count > 1) {
    const kept = count - step.droppedRuns;
    const expanded = kept > cap && open;
    const shown = step.lastRuns(expanded ? kept : cap);
    const first = count - shown.length;
    if (first > 0 || expanded) {
      made.push({ kind: 'more', key: `m${step.id}`, step, earlier: first, expanded, expandable: kept > cap, undone });
    }
    for (let i = 0; i < shown.length; i++) {
      made.push({ kind: 'run', key: `r${step.id}:${first + i}`, label: shown[i], undone });
    }
  }
  if (!previous?.length) return made;
  const old = new Map<string, HistoryRow>();
  for (const r of previous) old.set(r.key, r);
  for (let i = 0; i < made.length; i++) {
    const was = old.get(made[i].key);
    if (was && sameRow(was, made[i])) made[i] = was;
  }
  return made;
}

function sameRow(a: HistoryRow, b: HistoryRow): boolean {
  if (a.kind !== b.kind) return false;
  if (a.kind === 'start') return false;
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
  return a.step === c.step && a.undone === c.undone && a.current === c.current && a.reachable === c.reachable;
}

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
  /** The origins' icons by kind (`DefaultHistoryOriginIcons`); their names are the labels' (`origin`). */
  originIcons?: Readonly<Record<string, IconName>>;
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
 * on a step before it and "Redo to <step>" on one undone (muted after it), greyed when it cannot be reached (past an
 * unsafe step), an unsafe step's warning. Every action is a
 * model method (`undo`, `redo`, `moveTo`), so the same works headless. It fills its parent (a flex column).
 */
export function HistoryView(props: HistoryViewProps) {
  const { history } = props;
  const l = { ...DefaultHistoryViewLabels, ...props.labels };
  const icons = props.originIcons ?? DefaultHistoryOriginIcons;
  const undoKey = keyOf('undo', props.keys?.undo, props.commands);
  const redoKey = keyOf('redo', props.keys?.redo, props.commands);
  const expanded = new State<ReadonlySet<number>>(new Set());
  const toggle = (id: number) => {
    const next = new Set(expanded.value);
    if (!next.delete(id)) next.add(id);
    expanded.set(next);
  };
  const list = createHistoryRows({ runs: props.runs });
  const rows = derive([history.steps, history.redoable, expanded], list);
  const busy = derive(history.running, (n) => n > 0);
  const undoTip = derive(
    [history.undoBlocked, history.steps, undoKey],
    (blocked, steps, key) => blocked?.reason ?? l.undoTip(steps[0]?.label ?? '', key.hint),
  );
  const redoTip = derive([history.redoable, busy, redoKey], (r, b, key) =>
    b ? history.labels.running : r.length ? l.redoTip(r[0].label, key.hint) : history.labels.nothingToRedo,
  );
  // The reason Undo cannot run, when there is something it would undo (an unsafe step), but not a run in flight: that
  // is transient (the buttons are greyed meanwhile), and a note on every run would flash.
  const blocked = derive([history.undoBlocked, history.steps], (b, steps) =>
    steps.length && b && b.code !== 'running' ? b.reason : undefined,
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
                disabled={r.reachable ? busy : true}
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
              <Button
                variant="ghost"
                size="xs"
                pressed={r.expanded}
                aria-label={l.showRunsName(r.step.label)}
                onClick={() => toggle(r.step.id)}
              >
                {l.showRuns}
              </Button>
            ) : undefined
          }
        />
      );
    }
    const { step } = r;

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
            disabled={r.reachable ? busy : true}
            onClick={() => act(history.moveTo(step.id))}
          />
        )}
      </>
    );
    return (
      <ViewRow
        class={classes.row}
        role="listitem"
        icon={icons[step.origin.kind] ?? 'history'}
        label={step.label}
        description={l.description(l.origin(step.origin.kind), step.time)}
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
          aria-description={derive(history.undoBlocked, (b) => b?.reason)}
          aria-keyshortcuts={derive(undoKey, (k) => k.aria)}
          disabled={derive(history.undoBlocked, (b) => !!b)}
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
