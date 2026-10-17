/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { ReadonlyState } from '../../reactive/readonly-state.js';
import { derive } from '../../reactive/state.js';
import { Icon, type IconName } from '../icon/icon.js';
import { classes } from './_generated/choices.module.css.js';

export interface Choice {
  id: string;
  label: string;
  icon?: IconName;
  title?: string;
}

/**
 * A labelled row of one-click actions (quick styles: Default, Cartoon, Spacefill, Surface; wiggle: Dynamics,
 * Uncertainty, Clear). Equal widths; `active` marks the one in effect, if any.
 */
export function Choices(props: {
  label?: string;
  choices: readonly Choice[];
  active?: ReadonlyState<string | undefined>;
  onPick(id: string): void;
}) {
  const { active } = props;
  return (
    <div class={classes.group}>
      {props.label ? <span class={classes.label}>{props.label}</span> : null}
      <div class={classes.row} role="group" aria-label={props.label}>
        {props.choices.map((c) => (
          <button
            type="button"
            class={classes.choice}
            title={c.title}
            aria-pressed={active && derive(active, (a) => a === c.id)}
            onClick={() => props.onPick(c.id)}
          >
            {c.icon ? <Icon name={c.icon} size={15} /> : null}
            <span class={classes.text} data-tip-overflow="">
              {c.label}
            </span>
          </button>
        ))}
      </div>
    </div>
  );
}
