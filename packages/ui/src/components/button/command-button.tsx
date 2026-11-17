/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { derive } from '@molstar/core/reactive';
import { formatChord, isMac } from '../../commands/keys.js';
import type { CommandManager } from '../../commands/manager.js';
import type { CommandDefinition, CommandValues, KeyChord } from '../../commands/types.js';
import { type Child, onCleanup, Show } from '../../dom/index.js';
import { Icon } from '../icon/icon.js';
import { Button, type ButtonProps } from './button.js';

/**
 * A button that runs a command by id (commands §7): its icon, its title and key (after the keymap) in the tooltip,
 * greyed with the reason while it cannot run (`when`), gone while no such command is registered (an override removed
 * it). It opens the command through the manager: one without params runs, one with params opens its invocation (the
 * host shows it in the command menu).
 */
export function CommandButton(
  props: Pick<ButtonProps, 'variant' | 'size' | 'block' | 'class'> & {
    manager: CommandManager<any>;
    command: string;
    /** Values it opens with (the row it is for). */
    values?: CommandValues;
    /** A label beside the icon: the command's title (`true`) or other words; icon only without. */
    label?: true | Child;
    /** Told when the command refused to open (a run in flight, `when`). */
    onBlocked?(reason: string): void;
  },
) {
  const { manager } = props;
  // Wanted while shown: a lazy command it names is made in the registry's microtask, never inside this derive
  // (commands §4.2: a render never registers); the derive only looks it up.
  onCleanup(manager.want());
  const command = derive(manager.commands, () => manager.registered(props.command));
  return (
    <Show when={command}>
      {(c: CommandDefinition<any>) => {
        const look = derive(
          [manager.context, manager.keymap.layers],
          () => {
            const blocked = manager.blocked(c);
            const chord = manager.shortcut(c);
            const title = c.title.replace(/…$/, '');
            return {
              blocked,
              tip: blocked ?? (chord ? `${title} (${formatChord(chord)})` : title),
              keys: chord ? ariaChord(chord, isMac()) : undefined,
            };
          },
          (a, b) => a.blocked === b.blocked && a.tip === b.tip && a.keys === b.keys,
        );
        // The tooltip is the reason or the title with its key; the accessible name is the title alone (when there is
        // no label to name it), the reason its description, the key `aria-keyshortcuts`.
        const name = props.label === undefined ? c.title.replace(/…$/, '') : undefined;
        return (
          <Button
            {...(props.variant ? { variant: props.variant } : {})}
            {...(props.size ? { size: props.size } : {})}
            {...(props.block ? { block: true } : {})}
            {...(props.class ? { class: props.class } : {})}
            {...(c.icon ? { icon: <Icon name={c.icon} /> } : {})}
            title={derive(look, (s) => s.tip)}
            {...(name ? { 'aria-label': name } : {})}
            aria-description={derive(look, (s) => s.blocked)}
            aria-keyshortcuts={derive(look, (s) => s.keys)}
            disabled={derive(look, (s) => !!s.blocked)}
            onClick={() => {
              const r = manager.open(c.id, undefined, props.values);
              if (r.kind === 'blocked') props.onBlocked?.(r.reason);
            }}
          >
            {props.label === true ? c.title : props.label}
          </Button>
        );
      }}
    </Show>
  );
}

const AriaKeys: Readonly<Record<string, string>> = { ' ': 'Space', esc: 'Escape' };

/** A chord in `aria-keyshortcuts`' syntax (modifiers first, DOM key names): `Control+Shift+K`, `Meta+K`. */
function ariaChord(chord: KeyChord, mac: boolean): string {
  const keys: string[] = [];
  if (chord.mod) keys.push(mac ? 'Meta' : 'Control');
  if (chord.alt) keys.push('Alt');
  if (chord.shift) keys.push('Shift');
  const key = chord.key;
  keys.push(AriaKeys[key.toLowerCase()] ?? (key.length === 1 ? key.toUpperCase() : key));
  return keys.join('+');
}
