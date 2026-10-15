/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import styles from './choices.module.css';
import { Tooltip } from './tooltip.js';

export interface Choice {
  id: string;
  label: string;
  /** Its tooltip, when the label alone does not say what it does. */
  title?: string;
}

export interface ChoicesProps {
  /** What the row picks, above it (and the group's accessible name). */
  label?: string;
  choices: readonly Choice[];
  /** The one in effect, if any. */
  active?: string;
  onPick(id: string): void;
}

/**
 * A labelled row of one-click actions at equal widths, wrapping into rows of equal widths; `active` marks the one in
 * effect. Prism's `Choices` (components/choices), in React, in its tokens.
 */
export function Choices({ label, choices, active, onPick }: ChoicesProps) {
  return (
    <div className={styles.group} data-msp-box="">
      {label ? <span className={styles.label}>{label}</span> : null}
      <fieldset className={styles.row} aria-label={label}>
        {choices.map((c) => {
          const button = (
            <button
              key={c.id}
              type="button"
              className={styles.choice}
              aria-pressed={c.id === active}
              onClick={() => onPick(c.id)}
            >
              <span className={styles.text}>{c.label}</span>
            </button>
          );
          return c.title ? (
            <Tooltip key={c.id} content={c.title}>
              {button}
            </Tooltip>
          ) : (
            button
          );
        })}
      </fieldset>
    </div>
  );
}
