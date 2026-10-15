/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { map } from 'rxjs';
import { Icon } from '../../components/icon/icon.js';
import { For } from '../../dom/control-flow.js';
import type { StatusModel } from '../../shell/status.js';
import { classes } from './_generated/log-view.module.css.js';
import { KindIcons } from './canvas-chrome.js';

const time = new Intl.DateTimeFormat(undefined, { hour: '2-digit', minute: '2-digit', second: '2-digit' });

/** The Log tab: every message the status line showed, newest first, with its kind and time. */
export function LogView(props: { status: StatusModel }) {
  const newestFirst = props.status.history.pipe(map((h) => h.slice().reverse()));
  // Not live (a log is, by default): each message is said once, by the status line's region (README §3.7).
  return (
    <div class={classes.log} role="log" aria-live="off" aria-label="Log">
      <For each={newestFirst} key={(m) => m.id} fallback={() => <p class={classes.empty}>No messages yet.</p>}>
        {(m) => (
          <div class={classes.entry} data-kind={m.kind}>
            <Icon name={KindIcons[m.kind]} size={14} class={classes.kind} />
            <span class={classes.text}>{m.text}</span>
            <time class={classes.time} dateTime={m.time.toISOString()}>
              {time.format(m.time)}
            </time>
          </div>
        )}
      </For>
    </div>
  );
}
