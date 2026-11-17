/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { chordKeys } from '../../commands/keys.js';
import type { KeyChord } from '../../commands/types.js';
import { classes } from './_generated/footer.module.css.js';

/**
 * A key or a chord as caps: one `kbd` per key inside a `kbd` for the whole (the HTML way to write a combination), so
 * its text reads as the chord (`⌥J`). A string is one cap (`↵`, `esc`).
 */
export function Keys(props: { keys: KeyChord | string; mac?: boolean; class?: string }) {
  const caps = typeof props.keys === 'string' ? [props.keys] : chordKeys(props.keys, props.mac);
  return (
    <kbd class={props.class ? `${classes.keys} ${props.class}` : classes.keys}>
      {caps.map((k) => (
        <kbd class={classes.cap}>{k}</kbd>
      ))}
    </kbd>
  );
}

/** The keys for what is active, under the menu: a cap and what it does. */
export function Footer(props: { hints: readonly (readonly [keys: string, label: string])[] }) {
  return (
    <div class={classes.footer}>
      {props.hints.map(([k, label]) => (
        <span class={classes.hint}>
          <Keys keys={k} /> {label}
        </span>
      ))}
    </div>
  );
}
