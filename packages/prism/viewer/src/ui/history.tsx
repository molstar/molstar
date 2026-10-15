/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Button, domView, onActive, onCleanup, type ShellTab } from '@molstar/prism-ui';
import { observableState } from '@molstar/prism-ui/core';
import { Icon, type IconName } from '@molstar/prism-ui/components/icon/icon';
import { For, Show } from '@molstar/prism-ui/dom/control-flow';
import { distinctUntilChanged, map, type Observable } from 'rxjs';
import type { ViewerApi } from '../api/area.js';
import {
  createHistoryModel,
  type HistoryRow,
  type HistoryRunRow,
  type HistorySessionRow,
  type HistoryTurnRow,
} from '../journal/history-model.js';
import { createRecordIndex } from '../journal/record.js';
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

/**
 * "Keep this session's record" (journal.md §1): copies the record into IndexedDB and goes on writing there. Only the
 * record: no state, no files. Says so when it cannot (no IndexedDB), and returns whether it is kept.
 */
export async function keepSessionRecord(
  viewer: Pick<PrismViewer, 'record' | 'sources' | 'managers'>,
): Promise<boolean> {
  const { record, sources, managers } = viewer;
  if (record.kept) return true;
  try {
    // Loaded when asked: the viewer without a kept record never pays for it.
    const { createIndexedDbRows } = await import('@molstar/ai-components');
    const rows = createIndexedDbRows<import('@molstar/ai-fibers').FiberTable>(RecordsDatabase);
    await record.keep(rows, await createRecordIndex(rows), recordTitle(sources.structures.value));
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

export function historyTab(viewer: PrismViewer, api: ViewerApi): ShellTab {
  return {
    id: 'history',
    title: 'History',
    icon: 'history',
    order: 2.55,
    place: 'foot',
    render: domView(() => <HistoryView viewer={viewer} api={api} />),
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

  // Undo says what it undoes, or why it cannot.
  const title = history.undoBlocked.pipe(map((reason) => reason ?? 'Undo the latest step'));
  const undoLabel = history.steps.pipe(map((steps) => (steps.length ? `Undo ${steps[0].label}` : 'Undo')));
  const count = history.steps.pipe(
    map((steps) => (steps.length ? `${steps.length} step${steps.length === 1 ? '' : 's'}` : '')),
  );
  // "Undo turn" is for the agent's latest turn: the newest step that is its.
  const latestTurn = history.steps.pipe(
    map((steps) => steps.find((s) => s.origin.kind === 'agent')?.origin.turn),
    distinctUntilChanged(),
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
        <span class={classes.count}>{count}</span>
      </div>
      <div class={classes.keep}>
        <Show
          when={model.kept}
          fallback={() => (
            <Button
              size="xs"
              variant="ghost"
              icon={<Icon name="save" size={14} />}
              title="Keep this session's record in this browser (the record only: no scene, no files)"
              onClick={() => void keepSessionRecord(viewer)}
            >
              Keep this record
            </Button>
          )}
        >
          {() => <span class={classes.caption}>Kept in this browser</span>}
        </Show>
      </div>
      <For each={model.rows} key={(r) => r.key} fallback={() => <p class={classes.empty}>Nothing has happened yet.</p>}>
        {(row) => (
          <Row
            row={row}
            opened={opened}
            latestTurn={latestTurn}
            onUndoTurn={(t) => void api.history.undoTurn({ turn: t }).then(say)}
          />
        )}
      </For>
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
  onUndoTurn(turn: string): void;
}

function Row(props: RowProps) {
  const { row } = props;
  if (row.kind === 'turn')
    return <Turn row={row} opened={props.opened} latestTurn={props.latestTurn} onUndoTurn={props.onUndoTurn} />;
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
    </div>
  );
}
