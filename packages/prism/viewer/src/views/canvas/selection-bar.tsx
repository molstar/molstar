/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Icon } from '@molstar/prism-ui/components/icon/icon';
import { Show } from '@molstar/prism-ui/dom/control-flow';
import type { SelectionModel } from '../../model/selection.js';
import { classes } from './_generated/selection-bar.module.css.js';

/**
 * The selection, bottom right of the canvas while there is one: what is selected, then what to do with it (color,
 * theme, component, focus, …) as icon buttons, and clear. Mol*'s selection toolbar actions, without the toolbar.
 */
export function SelectionBar(props: { selection: SelectionModel }) {
  const { selection } = props;
  return (
    <Show when={selection.summary}>
      {(summary) => (
        <div class={classes.bar} role="toolbar" aria-label="Selection">
          <span class={classes.summary}>{summary}</span>
          {selection.actions.map((a) => (
            <button
              type="button"
              class={classes.action}
              aria-label={a.label}
              data-tip={a.label}
              data-tip-side="top"
              onClick={() => a.run()}
            >
              <Icon name={a.icon} size={16} />
            </button>
          ))}
          <button
            type="button"
            class={classes.action}
            aria-label="Clear selection"
            data-tip="Clear selection"
            data-tip-side="top"
            onClick={() => selection.clear()}
          >
            <Icon name="x" size={16} />
          </button>
        </div>
      )}
    </Show>
  );
}
