/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import {
  Button,
  Choices,
  domView,
  LogView,
  onActive,
  onCleanup,
  ScrollArea,
  type ShellTab,
  type ShellViewContext,
} from '@molstar/prism-ui';
import type { Child } from '@molstar/prism-ui/dom/jsx-types';
import { observableState } from '@molstar/prism-ui/core';
import { Icon, type IconName } from '@molstar/prism-ui/components/icon/icon';
import { For, Show } from '@molstar/prism-ui/dom/control-flow';
import { distinctUntilChanged, map, type Observable, of, shareReplay } from 'rxjs';
import type { ViewerApi } from '../api/area.js';
import {
  createHistoryModel,
  type HistoryRow,
  type HistoryRunRow,
  type HistorySessionRow,
  type HistoryTurnRow,
} from '../journal/history-model.js';
import type { OpenedRecord } from '../journal/kept-records.js';
import { createRecordIndex, type KeptRecord } from '../journal/record.js';
import { recordSpan } from '../views/home/kept-records.js';
import type { PrismViewer } from '../viewer.js';
import { classes } from './_generated/history.module.css.js';

// The History tab (README §3.3, journal.md §3): the session's record, newest first. A person's run is a row (who made it,
// what it did, when), an agent's turn a group (its prompt, its reply, its runs under it), an undone step is struck, a
// session event a small separator. Undo takes back the latest step, "Undo turn" the agent's latest turn. The model
// (`journal/history-model.ts`) reads the record while the tab shows and not otherwise.

const Origins: Readonly<Record<string, { icon: IconName; title: string }>> = {
  person: { icon: 'user', title: 'You' },
  agent: { icon: 'bot', title: 'The assistant' },
  script: { icon: 'file-code', title: 'A script' },
  url: { icon: 'download', title: 'A link to this page' },
  link: { icon: 'link', title: 'The sequence' },
};

const time = new Intl.DateTimeFormat(undefined, { hour: '2-digit', minute: '2-digit', second: '2-digit' });

/** The name IndexedDB keeps the record in (journal.md §1). */
const RecordsDatabase = 'prism-viewer-records';

/** What kept records are listed under: the loaded structures' names, else "Session". */
export function recordTitle(structures: readonly { label: string }[]): string {
  const names = [...new Set(structures.map((s) => s.label).filter(Boolean))];
  if (!names.length) return 'Session';
  return names.length > 3 ? `${names.slice(0, 3).join(', ')} and ${names.length - 3} more` : names.join(', ');
}

/** A keep in flight, per record: History's button and the command share it (one copy, one message). */
const keeping = new WeakMap<object, Promise<boolean>>();

/**
 * "Keep this session's record" (journal.md §1): copies the record into the spec's `persist.rows`, else the browser's
 * IndexedDB, and goes on writing there. Only the record: no state, no files. Says so when it cannot (no IndexedDB), and
 * returns whether it is kept.
 */
export function keepSessionRecord(viewer: Pick<PrismViewer, 'record' | 'sources' | 'managers'>): Promise<boolean> {
  const { record } = viewer;
  if (record.kept) return Promise.resolve(true);
  let p = keeping.get(record);
  if (!p) {
    p = keep(viewer).finally(() => keeping.delete(record));
    keeping.set(record, p);
  }
  return p;
}

async function keep(viewer: Pick<PrismViewer, 'record' | 'sources' | 'managers'>): Promise<boolean> {
  const { record, sources, managers } = viewer;
  try {
    // Loaded when asked: the viewer without a kept record never pays for it.
    // The spec's rows are the host's; the ones made here are the record's to close (when it is forgotten or ends).
    const made = record.keepIn
      ? undefined
      : (await import('../journal/indexeddb-rows.js')).createIndexedDbRows<import('@molstar/ai').FiberTable>(
          RecordsDatabase,
        );
    const rows = record.keepIn ?? made!;
    let handed = false;
    try {
      const index = await createRecordIndex(rows);
      const title = recordTitle(sources.structures.value);
      handed = true;
      await record.keep(rows, index, title, made && (() => made.close()));
    } finally {
      if (!handed) made?.close();
    }
    managers.notifications.post("This session's record is kept in this browser.", 'info');
    return true;
  } catch (e) {
    managers.notifications.post(
      `The record could not be kept: ${e instanceof Error ? e.message : String(e)}`,
      'warning',
    );
    return false;
  }
}

/**
 * History and the log in one panel (what happened: the steps, and what was said): on a desktop the two halves of it,
 * each scrolling on its own; compact, two tabs, the log's chosen when shown at `log` (the status line's log link).
 */
export function historyTab(
  viewer: PrismViewer,
  api: ViewerApi,
  status: Parameters<typeof LogView>[0]['status'],
): ShellTab {
  return {
    id: 'history',
    title: 'History',
    icon: 'history',
    order: 2.55,
    place: 'foot',
    render: domView<ShellViewContext>(({ ctx }) => {
      const compact = ctx.shell.mode.pipe(
        map((m) => m === 'compact'),
        distinctUntilChanged(),
      );
      // Split, both halves are always in view (the log at `log` too).
      const half = (title: string, body: () => Child) => (
        <section class={classes.half} aria-label={title}>
          <h3 class={classes.halfTitle}>{title}</h3>
          <ScrollArea class={classes.halfBody}>{body()}</ScrollArea>
        </section>
      );
      // A kept record opened from Home shows in its place, read-only, until "Back to this session".
      const steps = () => (
        <Show when={viewer.keptRecords.opened} fallback={() => <HistoryView viewer={viewer} api={api} />}>
          {(opened) => <KeptRecordView viewer={viewer} opened={opened} />}
        </Show>
      );
      const log = () => <LogView status={status} />;
      const shown = observableState<'steps' | 'log'>(ctx.focus === 'log' ? 'log' : 'steps');
      return (
        <Show
          when={compact}
          fallback={() => (
            <div class={classes.split}>
              {half('Steps', steps)}
              {half('Log', log)}
            </div>
          )}
        >
          {() => (
            <div class={classes.tabs}>
              <Choices
                choices={[
                  { id: 'steps', label: 'Steps', icon: 'history' },
                  { id: 'log', label: 'Log', icon: 'scroll-text' },
                ]}
                active={shown.asReadonly}
                onPick={(id) => shown.set(id as 'steps' | 'log')}
              />
              <Show when={shown.asReadonly.pipe(map((v) => v === 'log'))} fallback={steps}>
                {log}
              </Show>
            </div>
          )}
        </Show>
      );
    }),
  };
}

function HistoryView(props: { viewer: PrismViewer; api: ViewerApi }) {
  const { viewer, api } = props;
  const { history, notifications } = viewer.managers;
  // The record is read while the tab shows; a hidden tab reads nothing and holds no subscription.
  const model = createHistoryModel(viewer.record);
  onActive((active) => model.setActive(active));
  onCleanup(() => model.dispose());
  // Turns the person opened stay open when their row is made again (a reply arriving).
  const opened = new Set<string>();
  const keepBusy = observableState(false);
  // Keeping one more than the bound: the oldest kept record, to forget first (asked).
  const over = observableState<KeptRecord | undefined>(undefined);
  // Undo turn with later steps above it asks first (the history's `laterSteps`): the row says how many, and offers it.
  const confirm = observableState<{ turn: string; text: string } | undefined>(undefined);

  // Undo says what it undoes, or why it cannot.
  const title = history.undoBlocked.pipe(map((reason) => reason ?? 'Undo the latest step'));
  const undoLabel = history.steps.pipe(map((steps) => (steps.length ? `Undo ${steps[0].label}` : 'Undo')));
  const redoLabel = history.redoable.pipe(map((r) => (r.length ? `Redo ${r[0].label}` : 'Nothing to redo')));
  const count = history.steps.pipe(
    map((steps) => (steps.length ? `${steps.length} step${steps.length === 1 ? '' : 's'}` : '')),
  );
  // "Undo turn" is for the agent's latest turn: the newest step that is its.
  // Shared: computed once per change, not once per turn row that reads it.
  const latestTurn = history.steps.pipe(
    map((steps) => steps.find((s) => s.origin.kind === 'agent')?.origin.turn),
    distinctUntilChanged(),
    shareReplay({ bufferSize: 1, refCount: true }),
  );
  const say = (r: Awaited<ReturnType<ViewerApi['history']['undo']>>) => {
    if (!r.ok) notifications.post(r.message, 'warning');
  };
  return (
    <div class={classes.history} aria-label="History">
      <div class={classes.bar}>
        <Button
          size="sm"
          icon={<Icon name="undo-2" size={14} />}
          title={title}
          disabled={history.canUndo.pipe(map((can) => !can))}
          onClick={() => void api.history.undo().then(say)}
        >
          {undoLabel}
        </Button>
        <Button
          size="sm"
          icon={<Icon name="redo-2" size={14} />}
          title={redoLabel}
          aria-label="Redo"
          disabled={history.redoable.pipe(map((r) => r.length === 0))}
          onClick={() => void api.history.redo().then(say)}
        />
        <span class={classes.count}>{count}</span>
      </div>
      <div class={classes.keep}>
        <Show
          when={model.kept}
          fallback={() => (
            <Show
              when={over.asReadonly}
              fallback={() => (
                <Button
                  size="xs"
                  variant="ghost"
                  icon={<Icon name="save" size={14} />}
                  title="Keep this session's record in this browser (the record only: no scene, no files)"
                  disabled={keepBusy.asReadonly}
                  onClick={() => {
                    keepBusy.set(true);
                    void viewer.keptRecords
                      .overBound()
                      .then((oldest) => (oldest ? void over.set(oldest) : keepSessionRecord(viewer).then(() => {})))
                      .finally(() => keepBusy.set(false));
                  }}
                >
                  Keep this record
                </Button>
              )}
            >
              {(oldest) => (
                // At the bound (journal.md §5): keeping one more forgets the oldest, asked first.
                <div class={classes.ask}>
                  <span class={classes.caption}>
                    {`${viewer.keptRecords.max} records are kept: forget the oldest, “${oldest.title}”?`}
                  </span>
                  <Button
                    size="xs"
                    disabled={keepBusy.asReadonly}
                    onClick={() => {
                      keepBusy.set(true);
                      void viewer.keptRecords
                        .forget(oldest.id)
                        .then(() => keepSessionRecord(viewer))
                        .finally(() => {
                          keepBusy.set(false);
                          over.set(undefined);
                        });
                    }}
                  >
                    Forget it and keep
                  </Button>
                  <Button size="xs" variant="ghost" onClick={() => over.set(undefined)}>
                    Cancel
                  </Button>
                </div>
              )}
            </Show>
          )}
        >
          {() => (
            <>
              <span class={classes.caption}>Kept in this browser</span>
              <Button
                size="xs"
                variant="ghost"
                title="Forget the kept copy: it goes on in memory, and can be kept again"
                onClick={() => void viewer.keptRecords.forget(viewer.record.id)}
              >
                Forget
              </Button>
            </>
          )}
        </Show>
      </div>
      {/* The rows in a contained box of their own: the bar's changes do not lay them out, nor theirs the panel. */}
      <div class={classes.list}>
        <For
          each={model.rows}
          key={(r) => r.key}
          fallback={() => <p class={classes.empty}>Nothing has happened yet.</p>}
        >
          {(row) => (
            <Row
              row={row}
              opened={opened}
              latestTurn={latestTurn}
              confirm={confirm.asReadonly}
              onUndoTurn={(t, discardLater) => {
                confirm.set(undefined);
                void api.history.undoTurn({ turn: t, discardLater }).then((r) => {
                  if (!r.ok && r.code === 'confirm') confirm.set({ turn: t, text: r.message });
                  else say(r);
                });
              }}
            />
          )}
        </For>
      </div>
      <Show when={model.more}>
        {() => (
          <div class={classes.more}>
            <Button size="xs" variant="ghost" onClick={() => void model.showEarlier()}>
              Show earlier
            </Button>
          </div>
        )}
      </Show>
    </div>
  );
}

/**
 * A record kept in an earlier session, read-only (journal.md §3, "After a reload"): beside whatever is loaded now,
 * nothing in it can be undone or replayed. Its title and when, "Back to this session", Forget; its rows as History's.
 */
function KeptRecordView(props: { viewer: Pick<PrismViewer, 'keptRecords'>; opened: OpenedRecord }) {
  const { viewer, opened } = props;
  const model = createHistoryModel(opened.reader);
  onActive((active) => model.setActive(active));
  onCleanup(() => model.dispose());
  const opens = new Set<string>();
  const none = of(undefined);
  const asking = observableState(false);
  return (
    <div class={classes.history} aria-label="Kept record">
      <div class={classes.kept}>
        <div class={classes.keptText}>
          <span class={classes.keptTitle}>{opened.record.title}</span>
          <span class={classes.caption}>
            {opened.fromFile ? 'From the session file' : 'Kept record'}, read-only · {recordSpan(opened.record)}
          </span>
        </div>
        {/* A file's record is not kept: nothing to forget. */}
        {opened.fromFile ? null : (
          <Show
            when={asking.asReadonly}
            fallback={() => (
              <Button
                size="xs"
                variant="ghost"
                icon={<Icon name="trash-2" />}
                title="Forget this record"
                onClick={() => asking.set(true)}
              />
            )}
          >
            {() => (
              <Button
                size="xs"
                variant="ghost"
                title="Forget it: it is removed from this browser"
                onClick={() => void viewer.keptRecords.forget(opened.record.id)}
              >
                Forget?
              </Button>
            )}
          </Show>
        )}
        <Button
          size="xs"
          icon={<Icon name="chevron-left" />}
          title="Back to this session's history"
          onClick={() => viewer.keptRecords.close()}
        >
          This session
        </Button>
      </div>
      <div class={classes.list}>
        <For each={model.rows} key={(r) => r.key} fallback={() => <p class={classes.empty}>Nothing in it.</p>}>
          {(row) => <Row row={row} opened={opens} latestTurn={none} confirm={none} onUndoTurn={() => {}} />}
        </For>
      </div>
      <Show when={model.more}>
        {() => (
          <div class={classes.more}>
            <Button size="xs" variant="ghost" onClick={() => void model.showEarlier()}>
              Show earlier
            </Button>
          </div>
        )}
      </Show>
    </div>
  );
}

interface RowProps {
  row: HistoryRow;
  opened: Set<string>;
  latestTurn: Observable<string | undefined>;
  /** Undo turn asked first: the turn, and how many later steps would go with it. */
  confirm: Observable<{ turn: string; text: string } | undefined>;
  onUndoTurn(turn: string, discardLater?: boolean): void;
}

function Row(props: RowProps) {
  const { row } = props;
  if (row.kind === 'turn')
    return (
      <Turn
        row={row}
        opened={props.opened}
        latestTurn={props.latestTurn}
        confirm={props.confirm}
        onUndoTurn={props.onUndoTurn}
      />
    );
  if (row.kind === 'run') return <Run row={row} />;
  return <Separator row={row} />;
}

function Run(props: { row: HistoryRunRow; nested?: boolean }) {
  const { row } = props;
  const origin = Origins[row.origin] ?? { icon: 'user' as IconName, title: row.origin };
  return (
    <div
      class={classes.step}
      data-origin={row.origin}
      data-nested={props.nested ? '' : undefined}
      data-undone={row.undone ? '' : undefined}
      data-failed={row.failed ? '' : undefined}
    >
      <Icon name={row.failed ? 'circle-alert' : origin.icon} size={14} class={classes.origin} label={origin.title} />
      <span class={classes.label}>{row.text}</span>
      <time class={classes.time} dateTime={new Date(row.at).toISOString()}>
        {time.format(row.at)}
      </time>
      {row.undone || row.failed ? (
        <div class={classes.detail}>
          <span>{row.failed ? 'Failed' : 'Undone'}</span>
        </div>
      ) : null}
    </div>
  );
}

function Separator(props: { row: HistorySessionRow }) {
  const { row } = props;
  return (
    <div class={classes.separator}>
      <span class={classes.separatorText}>{row.text}</span>
      <time class={classes.time} dateTime={new Date(row.at).toISOString()}>
        {time.format(row.at)}
      </time>
    </div>
  );
}

const firstLine = (text: string) => text.trimStart().split('\n', 1)[0];

function Turn(props: Omit<RowProps, 'row'> & { row: HistoryTurnRow }) {
  const { row, opened } = props;
  const expanded = observableState(opened.has(row.key));
  const toggle = () => {
    const next = !expanded.value;
    if (next) opened.add(row.key);
    else opened.delete(row.key);
    expanded.set(next);
  };
  const reply = row.reply?.trim();
  const caption = row.open
    ? 'Under way'
    : row.undone
      ? 'Undone'
      : row.outcome === 'stopped'
        ? 'Stopped'
        : row.outcome === 'failed'
          ? 'Failed'
          : '';
  return (
    <div
      class={classes.turn}
      data-origin="agent"
      data-undone={row.undone ? '' : undefined}
      data-open={row.open ? '' : undefined}
    >
      <button
        type="button"
        class={classes.turnHead}
        aria-expanded={expanded.asReadonly.pipe(map(String))}
        onClick={toggle}
      >
        <Icon name="bot" size={14} class={classes.origin} label="The assistant" />
        <span class={classes.label}>{firstLine(row.prompt)}</span>
        <time class={classes.time} dateTime={new Date(row.at).toISOString()}>
          {time.format(row.at)}
        </time>
      </button>
      <Show
        when={expanded.asReadonly}
        fallback={() => (reply ? <p class={classes.reply}>{firstLine(reply)}</p> : null)}
      >
        {() => (
          <div class={classes.turnBody}>
            {firstLine(row.prompt) !== row.prompt.trim() ? <p class={classes.prompt}>{row.prompt}</p> : null}
            {reply ? <p class={classes.replyFull}>{reply}</p> : null}
            <div class={classes.runs}>
              {row.runs.map((r) => (
                <Run row={r} nested />
              ))}
            </div>
          </div>
        )}
      </Show>
      <div class={classes.detail}>
        {caption ? <span>{caption}</span> : null}
        {row.runs.length ? (
          <span>
            {row.runs.length} change{row.runs.length === 1 ? '' : 's'}
          </span>
        ) : null}
        <Show when={props.latestTurn.pipe(map((t) => t === row.turn && !row.undone && !row.open))}>
          {() => (
            <Button
              size="xs"
              variant="ghost"
              title="Undo everything this turn changed"
              onClick={() => props.onUndoTurn(row.turn)}
            >
              Undo turn
            </Button>
          )}
        </Show>
      </div>
      <Show when={props.confirm.pipe(map((c) => (c?.turn === row.turn ? c.text : undefined)))}>
        {(text) => (
          <div class={classes.detail}>
            <span class={classes.unsafe}>{text}.</span>
            <Button
              size="xs"
              variant="ghost"
              title="Undo this turn and every step made after it"
              onClick={() => props.onUndoTurn(row.turn, true)}
            >
              Undo them too
            </Button>
          </div>
        )}
      </Show>
    </div>
  );
}
