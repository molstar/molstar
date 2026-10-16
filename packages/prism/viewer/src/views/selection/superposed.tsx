/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ScrollArea } from '@molstar/prism-ui';
import { Icon } from '@molstar/prism-ui/components/icon/icon';
import { For, Show } from '@molstar/prism-ui/dom/control-flow';
import type { ReadonlyState } from '@molstar/prism-ui/reactive/readonly-state';
import { distinctUntilChanged, map } from 'rxjs';
import type { SuperposedEntry } from '../../model/superposition.js';
import { classes } from './_generated/selection-lists.module.css.js';

/**
 * The superpositions made (under the Superpose section's footer, as the Measure section's Added list): a row each, the
 * moved structure and how it was superposed, ✕ puts it back where it was loaded; Clear all with more than one.
 */
export function Superposed(props: {
  made: ReadonlyState<readonly SuperposedEntry[]>;
  remove(id: string): void;
  clear?(): void;
}) {
  const count = props.made.pipe(
    map((l) => l.length),
    distinctUntilChanged(),
  );
  return (
    <Show when={count.pipe(map((n) => n > 0))}>
      {() => (
        <div class={classes.block}>
          <div class={classes.head}>
            <span class={classes.caption}>
              Superposed <span class={classes.count}>{count}</span>
            </span>
            <Show when={count.pipe(map((n) => n > 1 && !!props.clear))}>
              {() => (
                <button type="button" class={classes.clear} onClick={() => props.clear?.()}>
                  Clear all
                </button>
              )}
            </Show>
          </div>
          <ScrollArea class={classes.list} contentClass={classes.listContent}>
            <For each={props.made} key={(e) => e.id}>
              {(entry) => (
                <div class={classes.row} title={entry.line}>
                  <span class={classes.lead}>
                    <Icon name="superpose" size={13} />
                  </span>
                  <span class={classes.name}>{entry.label}</span>
                  <span class={classes.detail}>{entry.detail}</span>
                  <span class={classes.actions}>
                    <button
                      type="button"
                      class={classes.action}
                      title="Remove the superposition"
                      aria-label={`Remove the superposition of ${entry.label}`}
                      onClick={() => props.remove(entry.id)}
                    >
                      <Icon name="x" size={13} />
                    </button>
                  </span>
                </div>
              )}
            </For>
          </ScrollArea>
        </div>
      )}
    </Show>
  );
}
