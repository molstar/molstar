/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Button, onActive, onCleanup } from '@molstar/prism-ui';
import { Icon } from '@molstar/prism-ui/components/icon/icon';
import { For, Show } from '@molstar/prism-ui/dom/control-flow';
import { observableState } from '@molstar/prism-ui/reactive/observable-state';
import type { KeptRecord } from '../../journal/record.js';
import type { KeptRecordsHomeModel } from '../../model/home.js';
import { classes } from './_generated/kept-records.module.css.js';

import { derive } from '@molstar/prism-ui/reactive/state';
// Home's kept records (journal.md §3, "After a reload"): the records kept in earlier sessions, newest first, each its
// title and when it was written; a click opens it in History, read-only, beside whatever is loaded; Forget asks once on
// its row. Listed only while Home shows; with none, nothing.

const day = new Intl.DateTimeFormat(undefined, { month: 'short', day: 'numeric' });
const time = new Intl.DateTimeFormat(undefined, { hour: '2-digit', minute: '2-digit' });

/** "Oct 14, 10:02–11:30", or across days "Oct 12 – Oct 14". */
export function recordSpan(r: Pick<KeptRecord, 'first' | 'last'>): string {
  const a = new Date(r.first);
  const b = new Date(r.last);
  if (a.toDateString() !== b.toDateString()) return `${day.format(a)} – ${day.format(b)}`;
  const from = time.format(a);
  const to = time.format(b);
  return from === to ? `${day.format(a)}, ${from}` : `${day.format(a)}, ${from}–${to}`;
}

export function KeptRecordsList(props: { model: KeptRecordsHomeModel }) {
  const { model } = props;
  let stop: (() => void) | undefined;
  onActive((active) => {
    if (active) stop ??= model.watch();
    else {
      stop?.();
      stop = undefined;
    }
  });
  onCleanup(() => stop?.());
  // The row asking "Forget?" (one at a time).
  const asking = observableState<string | undefined>(undefined);
  return (
    <Show when={derive(model.list, (l) => l.length > 0)}>
      {() => (
        <section class={classes.group} aria-label="Kept records">
          <h3 class={classes.caption}>Kept records</h3>
          <div class={classes.listContent} role="list">
            <For each={model.list} key={(r) => r.id}>
              {(r) => {
                const isAsking = derive(asking.asReadonly, (id) => id === r.id);
                return (
                  <div class={classes.record} role="listitem">
                    <button
                      type="button"
                      class={classes.open}
                      title="Open in History (read-only)"
                      onClick={() => model.open(r.id)}
                    >
                      <span class={classes.title} data-tip-overflow="">
                        {r.title}
                      </span>
                      <span class={classes.when}>{recordSpan(r)}</span>
                    </button>
                    <Show
                      when={isAsking}
                      fallback={() => (
                        <Button
                          size="xs"
                          variant="ghost"
                          class={classes.forget}
                          icon={<Icon name="trash-2" />}
                          title="Forget this record"
                          onClick={() => asking.set(r.id)}
                        />
                      )}
                    >
                      {() => (
                        <Button
                          size="xs"
                          variant="ghost"
                          class={`${classes.forget} ${classes.asking}`}
                          title="Forget it: it is removed from this browser"
                          onClick={() => {
                            asking.set(undefined);
                            void model.forget(r.id);
                          }}
                        >
                          Forget?
                        </Button>
                      )}
                    </Show>
                  </div>
                );
              }}
            </For>
          </div>
        </section>
      )}
    </Show>
  );
}
