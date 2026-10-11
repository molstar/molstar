/** @jsxImportSource @molstar/plugin-prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { combineLatest, distinctUntilChanged, filter, map, type Observable } from 'rxjs';
import { Icon } from '../../components/icon/icon.js';
import { ScrollArea } from '../../components/scroll-area/scroll-area.js';
import { For, Show } from '../../dom/control-flow.js';
import type { AgentModel, AgentTurn, CommandRun, HistoryItem } from '../../model/agent.js';
import { observableState } from '../../reactive/observable-state.js';
import type { StatusModel } from '../../shell/status.js';
import { LogView } from '../shell/log-view.js';
import { isRunning, RunList, RunRow, turnSummary, UndoTurn } from './agent-turn.js';
import { classes } from './_generated/history-tab.module.css.js';

const keyOf = (h: HistoryItem) => (h.kind === 'run' ? h.run.id : h.turn.id);

function ago(time: Date) {
  const s = Math.max(0, (Date.now() - time.getTime()) / 1000);
  if (s < 45) return 'just now';
  if (s < 3600) return `${Math.round(s / 60)} min ago`;
  return time.toLocaleTimeString(undefined, { hour: '2-digit', minute: '2-digit' });
}

/**
 * The History tab (the Log grown up, design/agent.md §5): every command, newest first, the person's one by one and an
 * agent's grouped under its request, with the checkpoint before it (Undo turn); the status messages under Messages.
 */
export function HistoryTab(props: { agent: AgentModel; status: StatusModel }) {
  const { agent } = props;
  const view = observableState<'commands' | 'messages'>('commands');
  // By id, newest first: an item's record is replaced as its turn runs, and an item made anew would lose its fold.
  const ids = agent.history.pipe(
    map((h) => h.map(keyOf).reverse()),
    distinctUntilChanged((a, b) => a.length === b.length && a.every((x, i) => x === b[i])),
  );
  const live = (id: string) =>
    agent.history.pipe(
      map((h) => h.find((x) => keyOf(x) === id)),
      filter((x): x is HistoryItem => !!x),
      distinctUntilChanged(),
    );
  return (
    <div class={classes.tab}>
      <div class={classes.switch} role="tablist" aria-label="History">
        {(['commands', 'messages'] as const).map((v) => (
          <button
            type="button"
            role="tab"
            class={classes.switchButton}
            aria-selected={view.asReadonly.pipe(map((x) => x === v))}
            onClick={() => view.set(v)}
          >
            {v === 'commands' ? 'Commands' : 'Messages'}
          </button>
        ))}
      </div>
      <Show
        when={view.asReadonly.pipe(map((v) => v === 'commands'))}
        fallback={() => <LogView status={props.status} />}
      >
        {() => (
          <ScrollArea class={classes.list} contentClass={classes.listContent} role="log" aria-label="Commands">
            <For
              each={ids}
              fallback={() => (
                <p class={classes.empty}>Nothing run yet. Commands and the assistant’s work show here.</p>
              )}
            >
              {(id) => {
                const item = agent.history.value.find((x) => keyOf(x) === id)!;
                return item.kind === 'run' ? (
                  <div class={classes.item}>
                    <RunRow run={live(id).pipe(map((x) => (x as { kind: 'run'; run: CommandRun }).run))} />
                    <span class={classes.when}>{ago(item.run.time)}</span>
                  </div>
                ) : (
                  <TurnItem
                    agent={agent}
                    turn={live(id).pipe(map((x) => (x as { kind: 'turn'; turn: AgentTurn }).turn))}
                    id={id}
                    time={item.turn.time}
                  />
                );
              }}
            </For>
          </ScrollArea>
        )}
      </Show>
    </div>
  );
}

/** An agent's turn: who and what was asked, its answer, its commands (folded once done), and its checkpoint. */
function TurnItem(props: { agent: AgentModel; turn: Observable<AgentTurn>; id: string; time: Date }) {
  const { agent, turn, id } = props;
  const open = observableState(false);
  const running = turn.pipe(map(isRunning), distinctUntilChanged());
  const undone = turn.pipe(
    map((t) => !!t.undone),
    distinctUntilChanged(),
  );
  return (
    <section
      class={classes.turn}
      data-undone={undone.pipe(map((u) => (u ? '' : undefined)))}
      aria-label="Assistant turn"
    >
      <div class={classes.turnHead}>
        <Icon name="sparkles" size={14} class={classes.spark} />
        <span class={classes.who}>
          {turn.pipe(map((t) => (t.source.kind === 'agent' ? t.source.name : 'Assistant')))}
        </span>
        <span class={classes.when}>{ago(props.time)}</span>
      </div>
      <p class={classes.prompt}>“{turn.pipe(map((t) => t.prompt))}”</p>
      <Show when={turn.pipe(map((t) => (isRunning(t) ? '' : t.text)))}>
        {(text) => <p class={classes.answer}>{text}</p>}
      </Show>
      <button type="button" class={classes.fold} aria-expanded={open.asReadonly} onClick={() => open.set(!open.value)}>
        <span>{turn.pipe(map(turnSummary))}</span>
        <Icon name="chevron-down" size={14} class={classes.chevron} />
      </button>
      {/* Its commands: shown while it runs, folded once it is done. */}
      <Show when={combineLatest([open.asReadonly, running]).pipe(map(([o, r]) => o || r))} reveal>
        {() => <RunList turn={turn} class={classes.runs} />}
      </Show>
      {/* Its checkpoint, held until the next action: Undo turn goes back to it; after that, it says why it cannot. */}
      <Show when={turn.pipe(map((t) => !isRunning(t)))}>
        {() => (
          <div class={classes.checkpoint}>
            <Icon name="undo-2" size={12} class={classes.checkpointMark} />
            <span class={classes.checkpointText}>
              {turn.pipe(
                map((t) =>
                  t.undone
                    ? 'Undone: back to how it was before'
                    : t.checkpoint
                      ? 'Checkpoint before it'
                      : t.runs.some((r) => !r.read)
                        ? 'No longer undoable: something came after it'
                        : 'Nothing to undo',
                ),
              )}
            </span>
            <Show when={turn.pipe(map((t) => t.checkpoint && !t.undone))}>
              {() => <UndoTurn agent={agent} id={id} />}
            </Show>
          </div>
        )}
      </Show>
    </section>
  );
}
