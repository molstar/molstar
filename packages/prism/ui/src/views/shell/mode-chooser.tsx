/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { map, type Observable } from 'rxjs';
import { Icon, type IconName } from '../../components/icon/icon.js';
import { For } from '../../dom/control-flow.js';
import { classes } from './_generated/mode-chooser.module.css.js';

/** A mode as the chooser offers it (design/modes.md). */
export interface ModeChoice {
  id: string;
  title: string;
  icon: IconName;
  /** What it is for, in a line ("Load structures and explore them"). */
  description?: string;
}

/**
 * The app's modes to choose from: `welcome`, a card over the canvas the first time the app opens (each mode a large
 * choice); `rows`, compact, where nothing is loaded yet (Home). The current mode is marked.
 */
export function ModeChooser(props: {
  modes: readonly ModeChoice[];
  active?: Observable<string>;
  onChoose(id: string): void;
  variant: 'welcome' | 'rows';
}) {
  const options = (
    <div class={classes.options} role="radiogroup" aria-label="Mode">
      <For each={props.modes}>
        {(mode) => (
          <button
            type="button"
            role="radio"
            class={classes.option}
            aria-checked={props.active?.pipe(map((a) => a === mode.id)) ?? false}
            onClick={() => props.onChoose(mode.id)}
          >
            <Icon name={mode.icon} size={props.variant === 'welcome' ? 24 : 16} class={classes.icon} />
            <span class={classes.title}>{mode.title}</span>
            {mode.description ? <span class={classes.description}>{mode.description}</span> : null}
          </button>
        )}
      </For>
    </div>
  );
  if (props.variant === 'rows') return <div class={classes.rows}>{options}</div>;
  return (
    <div class={classes.welcome}>
      <div class={classes.card} role="dialog" aria-label="Choose how to start">
        <h2 class={classes.heading}>How do you want to start?</h2>
        {options}
        <p class={classes.note}>You can switch any time, from the rail or from Home.</p>
      </div>
    </div>
  );
}
