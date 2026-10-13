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
import { MolstarWordmark } from '../../components/icon/wordmark.js';
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
  variant: 'splash' | 'welcome' | 'rows';
  /** What choosing means here ("what is loaded now goes"). */
  notice?: string;
  /** Not switching after all (over a running app). */
  onCancel?: () => void;
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
            <Icon name={mode.icon} size={props.variant === 'rows' ? 16 : 24} class={classes.icon} />
            <span class={classes.title}>{mode.title}</span>
            {mode.description ? <span class={classes.description}>{mode.description}</span> : null}
          </button>
        )}
      </For>
    </div>
  );
  if (props.variant === 'rows') return <div class={classes.rows}>{options}</div>;
  const splash = props.variant === 'splash';
  return (
    <div class={splash ? classes.splash : classes.welcome}>
      {splash ? <MolstarWordmark scale={3} flash="reveal" /> : null}
      <div class={classes.card} role="dialog" aria-label={splash ? 'Choose how to start' : 'Switch mode'}>
        <h2 class={classes.heading}>{splash ? 'How do you want to start?' : 'Switch mode'}</h2>
        {options}
        <p class={classes.note}>{props.notice ?? 'You can switch later, from the command menu or from Home.'}</p>
        {props.onCancel ? (
          <button type="button" class={classes.cancel} onClick={() => props.onCancel?.()}>
            Cancel
          </button>
        ) : null}
      </div>
    </div>
  );
}
