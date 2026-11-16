/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { derive } from '@molstar/core/reactive';
import { formatChord } from '../../commands/keys.js';
import type { CommandManager } from '../../commands/manager.js';
import type { CommandDefinition, CommandValues } from '../../commands/types.js';
import { type Child, Show } from '../../dom/index.js';
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
  const command = derive(manager.commands, () => manager.get(props.command));
  return (
    <Show when={command}>
      {(c: CommandDefinition<any>) => {
        const look = derive(
          [manager.context, manager.keymap.layers],
          () => {
            const blocked = manager.blocked(c);
            const chord = manager.shortcut(c);
            const title = c.title.replace(/…$/, '');
            return { blocked, tip: blocked ?? (chord ? `${title} (${formatChord(chord)})` : title) };
          },
          (a, b) => a.blocked === b.blocked && a.tip === b.tip,
        );
        return (
          <Button
            {...(props.variant ? { variant: props.variant } : {})}
            {...(props.size ? { size: props.size } : {})}
            {...(props.block ? { block: true } : {})}
            {...(props.class ? { class: props.class } : {})}
            {...(c.icon ? { icon: <Icon name={c.icon} /> } : {})}
            title={derive(look, (s) => s.tip)}
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
