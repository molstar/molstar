/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { combineLatest, distinctUntilChanged, map } from 'rxjs';
import { formatChord, isMac } from '../../commands/keys.js';
import type { CommandManager } from '../../commands/manager.js';
import type { CommandDefinition, CommandValues } from '../../commands/types.js';
import { Show } from '../../dom/control-flow.js';
import type { Child } from '../../dom/jsx-types.js';
import { Icon } from '../icon/icon.js';
import { Button, type ButtonProps } from './button.js';

/**
 * A button that runs a command by id (commands.md §7): its icon, its title and key in the tooltip, greyed with the
 * reason while the command cannot run (`when`), and gone while no such command is registered (an override removed it).
 * It opens the command through the manager: one without params runs, one with params opens its form (the host shows
 * it in the command menu).
 */
export function CommandButton(
  props: Pick<ButtonProps, 'variant' | 'size' | 'block' | 'class'> & {
    manager: CommandManager<any>;
    command: string;
    /** Values it opens with (the row it is for). */
    values?: CommandValues;
    /** A label beside the icon: the command's title (`true`) or other words; icon only without. */
    label?: true | Child;
    /** Told when the command refused to open. */
    onBlocked?(reason: string): void;
  },
) {
  const { manager } = props;
  const command = manager.commands.pipe(
    map((all) => all.find((c) => c.id === props.command)),
    distinctUntilChanged(),
  );
  return (
    <Show when={command}>
      {(c: CommandDefinition<any>) => {
        const state = combineLatest([manager.context, manager.keymap.layers]).pipe(
          map(() => {
            const blocked = manager.blocked(c);
            const chord = manager.shortcut(c);
            const title = c.title.replace(/…$/, '');
            return { blocked, tip: blocked ?? (chord ? `${title} (${formatChord(chord, isMac())})` : title) };
          }),
          distinctUntilChanged((a, b) => a.blocked === b.blocked && a.tip === b.tip),
        );
        return (
          <Button
            {...(props.variant ? { variant: props.variant } : {})}
            {...(props.size ? { size: props.size } : {})}
            {...(props.block ? { block: true } : {})}
            {...(props.class ? { class: props.class } : {})}
            {...(c.icon ? { icon: <Icon name={c.icon} /> } : {})}
            title={state.pipe(map((s) => s.tip))}
            disabled={state.pipe(map((s) => !!s.blocked))}
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
