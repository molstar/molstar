/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import {
  Icon,
  ScrollArea,
  For,
  Show,
  observableState,
  type NotificationsManager,
  LogView,
  type ReadonlyState,
} from '@molstar/prism-ui';
import type { AgentModel, AgentTurn, HistoryItem } from '../model.js';
import { held, isRunning, RunList, RunRow, turnSummary, UndoTurn } from './agent-turn.js';
import { classes } from './_generated/history-tab.module.css.js';

import { derive } from '@molstar/prism-ui/reactive/state';
import { effect, watch } from '@molstar/prism-ui/dom/scope';
import { setAttr, setText } from '@molstar/prism-ui/dom/writes';
import { keyedState } from '@molstar/prism-ui/reactive/keyed-state';
const keyOf = (h: HistoryItem) => (h.kind === 'run' ? h.run.id : h.turn.id);

const clock = new Intl.DateTimeFormat(undefined, { hour: '2-digit', minute: '2-digit' });

/**
 * When an item ran, as the time of day (as the Log says it): a relative "just now" / "N min ago" was written once, when
 * the row was made, and never again, so it went stale while newer rows said later times (perf audit 2, 11/B4). A clock
 * ticking every row would keep a loaded, idle app busy.
 */
function When(props: { time: Date }) {
  return (
    <time class={classes.when} dateTime={props.time.toISOString()} title={props.time.toLocaleString()}>
      {clock.format(props.time)}
    </time>
  );
}

/**
 * The History tab (the Log grown up, _spec/README.md, "History"): every command, newest first, the person's one by one
 * and an agent's grouped under its request, with the checkpoint before it (Undo turn); the status messages under
 * Messages.
 */
export function HistoryTab(props: { agent: AgentModel; status: NotificationsManager }) {
  const { agent } = props;
  const view = observableState<'commands' | 'messages'>('commands');
  // By id, newest first: an item's record is replaced as its turn runs, and an item made anew would lose its fold.
  const ids = derive(
    agent.history,
    (h) => h.map(keyOf).reverse(),
    (a, b) => a.length === b.length && a.every((x, i) => x === b[i]),
  );
  // The items by key, one map a change; an item is told only when its own record changes (a `find` an item, for each of
  // its fields, on every publish was items² a frame while a turn streams).
  const byKey = keyedState<string, HistoryItem>();
  effect(agent.history, (h) => {
    const m = new Map<string, HistoryItem>();
    for (const x of h) m.set(keyOf(x), x);
    byKey.replace(m);
  });
  const live = (id: string) => held(byKey.key(id), (x): x is HistoryItem => !!x);
  const tabs = (['commands', 'messages'] as const).map((v) => ({
    v,
    button: (
      <button type="button" role="tab" class={classes.switchButton} onClick={() => view.set(v)}>
        {v === 'commands' ? 'Commands' : 'Messages'}
      </button>
    ) as HTMLElement,
  }));
  watch([view], () => {
    for (const t of tabs) setAttr(t.button, 'aria-selected', String(t.v === view.value));
  });
  return (
    <div class={classes.tab}>
      <div class={classes.switch} role="tablist" aria-label="History">
        {tabs.map((t) => t.button)}
      </div>
      <Show when={derive(view.asReadonly, (v) => v === 'commands')} fallback={() => <LogView status={props.status} />}>
        {() => (
          <ScrollArea
            class={classes.list}
            contentClass={classes.listContent}
            role="log"
            aria-live="off"
            aria-label="Commands"
          >
            <For
              each={ids}
              fallback={() => (
                <p class={classes.empty}>Nothing run yet. Commands and the assistant’s work show here.</p>
              )}
            >
              {(id) => {
                const item = byKey.get(id) ?? agent.history.value.find((x) => keyOf(x) === id)!;
                return item.kind === 'run' ? (
                  <div class={classes.item}>
                    <RunRow run={derive(live(id), (x) => (x as Extract<HistoryItem, { kind: 'run' }>).run)} />
                    <When time={item.run.time} />
                  </div>
                ) : (
                  <TurnItem
                    agent={agent}
                    turn={derive(live(id), (x) => (x as Extract<HistoryItem, { kind: 'turn' }>).turn)}
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
function TurnItem(props: { agent: AgentModel; turn: ReadonlyState<AgentTurn>; id: string; time: Date }) {
  const { agent, turn, id } = props;
  const open = observableState(false);
  const running = derive(turn, isRunning);
  // What the item's header and fold say, each from the one field of the record it is made of.
  const undone = derive(turn, (t) => !!t.undone);
  const who = derive(turn, (t) => (t.source.kind === 'agent' ? t.source.name : 'Assistant'));
  const prompt = derive(turn, (t) => t.prompt);
  const summary = derive(turn, turnSummary);
  const whoText = document.createTextNode('');
  const promptText = document.createTextNode('');
  const summaryText = document.createTextNode('');
  let fold: HTMLElement | undefined;
  const section = (
    <section class={classes.turn} aria-label="Assistant turn">
      <div class={classes.turnHead}>
        <Icon name="sparkles" size={14} class={classes.spark} />
        <span class={classes.who}>{whoText}</span>
        <When time={props.time} />
      </div>
      <p class={classes.prompt}>“{promptText}”</p>
      <Show when={derive(turn, (t) => (isRunning(t) ? '' : t.text))}>
        {(text) => <p class={classes.answer}>{text}</p>}
      </Show>
      <button
        type="button"
        class={classes.fold}
        ref={(el: HTMLElement) => (fold = el)}
        onClick={() => open.set(!open.value)}
      >
        <span>{summaryText}</span>
        <Icon name="chevron-down" size={14} class={classes.chevron} />
      </button>
      {/* Its commands: shown while it runs, folded once it is done. */}
      <Show when={derive([open.asReadonly, running], (o, r) => o || r)} reveal>
        {() => <RunList turn={turn} class={classes.runs} />}
      </Show>
      {/* Its checkpoint, held until the next action: Undo turn goes back to it; after that, it says why it cannot. */}
      <Show when={derive(turn, (t) => !isRunning(t))}>{() => <Checkpoint agent={agent} turn={turn} id={id} />}</Show>
    </section>
  ) as HTMLElement;
  watch([undone, who, prompt, summary, open], () => {
    setAttr(section, 'data-undone', undone.value);
    setText(whoText, who.value);
    setText(promptText, prompt.value);
    setText(summaryText, summary.value);
    if (fold) setAttr(fold, 'aria-expanded', String(open.value));
  });
  return section;
}

/** A finished turn's checkpoint: what can be done about it, and Undo turn while it can. */
function Checkpoint(props: { agent: AgentModel; turn: ReadonlyState<AgentTurn>; id: string }) {
  const { agent, turn, id } = props;
  const note = derive(turn, (t) =>
    t.undone
      ? 'Undone: back to how it was before'
      : t.checkpoint
        ? 'Checkpoint before it'
        : t.runs.some((r) => !r.read)
          ? 'No longer undoable: something came after it'
          : 'Nothing to undo',
  );
  const text = document.createTextNode('');
  watch([note], () => setText(text, note.value));
  return (
    <div class={classes.checkpoint}>
      <Icon name="undo-2" size={12} class={classes.checkpointMark} />
      <span class={classes.checkpointText}>{text}</span>
      <Show when={derive(turn, (t) => t.checkpoint && !t.undone)}>{() => <UndoTurn agent={agent} id={id} />}</Show>
    </div>
  );
}
