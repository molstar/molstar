/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Button, domView, type ShellTab } from '@molstar/prism-ui';
import type { HistoryStep } from '@molstar/prism-ui/core';
import { Icon, type IconName } from '@molstar/prism-ui/components/icon/icon';
import { For, Show } from '@molstar/prism-ui/dom/control-flow';
import { map } from 'rxjs';
import type { ViewerApi } from '../api/area.js';
import type { ViewerManagers } from '../api/managers.js';
import { classes } from './_generated/history.module.css.js';

// The History tab (README §3.3, sequence.md §7.3): the undoable steps, newest first, each with who made it and when;
// Undo for the latest, "Undo turn" for an agent's turn. Steps are keyed by id, so a row is made once and again only when
// its step changes (an agent's turn takes more runs, an outside change makes it unsafe). The journal (phase 4) adds the
// record's rows here: the calls inside a step, their results.

const Origins: Readonly<Record<string, { icon: IconName; title: string }>> = {
  person: { icon: 'user', title: 'You' },
  agent: { icon: 'bot', title: 'The assistant' },
  script: { icon: 'file-code', title: 'A script' },
  url: { icon: 'download', title: 'A link to this page' },
  link: { icon: 'link', title: 'The sequence' },
};

const time = new Intl.DateTimeFormat(undefined, { hour: '2-digit', minute: '2-digit', second: '2-digit' });

export function historyTab(managers: ViewerManagers, api: ViewerApi): ShellTab {
  return {
    id: 'history',
    title: 'History',
    icon: 'history',
    order: 2.55,
    place: 'foot',
    render: domView(() => <HistoryView managers={managers} api={api} />),
  };
}

function HistoryView(props: { managers: ViewerManagers; api: ViewerApi }) {
  const { managers, api } = props;
  const { history, notifications } = managers;
  // Undo says what it undoes, or why it cannot.
  const title = history.undoBlocked.pipe(map((reason) => reason ?? 'Undo the latest step'));
  const undoLabel = history.steps.pipe(map((steps) => (steps.length ? `Undo ${steps[0].label}` : 'Undo')));
  const count = history.steps.pipe(
    map((steps) => (steps.length ? `${steps.length} step${steps.length === 1 ? '' : 's'}` : '')),
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
      <For each={history.steps} key={(s) => s.id} fallback={() => <p class={classes.empty}>Nothing to undo yet.</p>}>
        {(step) => <Step step={step} onUndoTurn={(turn) => void api.history.undoTurn({ turn }).then(say)} />}
      </For>
    </div>
  );
}

function Step(props: { step: HistoryStep; onUndoTurn(turn: string): void }) {
  const { step } = props;
  const origin = Origins[step.origin.kind] ?? { icon: 'user' as IconName, title: step.origin.kind };
  const turn = step.origin.kind === 'agent' ? step.origin.turn : undefined;
  return (
    <div class={classes.step} data-origin={step.origin.kind} data-unsafe={step.unsafe ? '' : undefined}>
      <Icon name={origin.icon} size={14} class={classes.origin} label={origin.title} />
      <span class={classes.label}>{step.label}</span>
      <time class={classes.time} dateTime={new Date(step.time).toISOString()}>
        {time.format(step.time)}
      </time>
      <Show when={step.unsafe || turn || step.runs.length > 1}>
        {() => (
          <div class={classes.detail}>
            {step.runs.length > 1 ? <span>{step.runs.length} changes</span> : null}
            {step.unsafe ? <span class={classes.unsafe}>{step.unsafe}</span> : null}
            {turn ? (
              <Button
                size="xs"
                variant="ghost"
                title="Undo everything this turn changed"
                onClick={() => props.onUndoTurn(turn)}
              >
                Undo turn
              </Button>
            ) : null}
          </div>
        )}
      </Show>
    </div>
  );
}
