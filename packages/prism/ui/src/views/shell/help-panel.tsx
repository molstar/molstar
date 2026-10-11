/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { map } from 'rxjs';
import { chordKeys } from '../../commands/keys.js';
import type { CommandManager } from '../../commands/manager.js';
import type { CommandDefinition } from '../../commands/types.js';
import { Icon } from '../../components/icon/icon.js';
import { ScrollArea } from '../../components/scroll-area/scroll-area.js';
import { For, Show } from '../../dom/control-flow.js';
import { classes } from './_generated/help-panel.module.css.js';

/** A mouse or touch binding, as the Help card lists it: what it does, and the input (keys and a gesture). */
export interface InputBinding {
  action: string;
  /** The input, as keys and gestures in order: ['Shift', 'Click'], ['Right drag']. */
  input: readonly string[];
  /** Another way to do the same (a trackpad gesture, a touch). */
  or?: readonly string[];
}

/** Bindings that belong together ("Camera", "Selection"): a mode contributes its own. */
export interface InputBindingGroup {
  title: string;
  bindings: readonly InputBinding[];
}

/** The keys that work everywhere, whatever the mode. */
const GlobalKeys: InputBinding[] = [
  { action: 'Command menu', input: ['⌘', 'K'], or: ['⇧', '⌘', 'P'] },
  { action: 'Help (this card)', input: ['?'] },
  { action: 'Close a menu, a card or the command menu', input: ['Esc'] },
];

/**
 * Help: how to work the canvas (the mode's mouse and touch bindings), the keys (global ones, then every command with a
 * shortcut, read from the command registry as it is), and shortcut conflicts (two commands on one chord). A tab at the
 * rail's foot (it stays open while you try things on the canvas); also the command menu's `?`.
 */
export function HelpPanel(props: { commands: CommandManager<any>; bindings: readonly InputBindingGroup[] }) {
  const shortcuts = props.commands.commands.pipe(map((all) => all.filter((c) => c.shortcut)));
  const conflicts = shortcuts.pipe(map(conflictsOf));
  return (
    <ScrollArea edge class={classes.help} contentClass={classes.content} aria-label="Help">
      <Show when={conflicts.pipe(map((c) => (c.length ? c : undefined)))}>
        {(list) => (
          <section class={classes.conflicts} aria-label="Shortcut conflicts">
            <Icon name="triangle-alert" size={14} />
            <div>
              {list.map(([keys, names]) => (
                <p class={classes.conflict}>
                  <Keys input={keys} /> is used by {names.join(' and ')}: only the first runs.
                </p>
              ))}
            </div>
          </section>
        )}
      </Show>
      <For each={props.bindings}>{(group) => <Group title={group.title} bindings={group.bindings} />}</For>
      <Group title="Everywhere" bindings={GlobalKeys} />
      <Show when={shortcuts.pipe(map((s) => (s.length ? s : undefined)))}>
        {(list) => (
          <Group title="Commands" bindings={list.map((c) => ({ action: c.title, input: chordKeys(c.shortcut!) }))} />
        )}
      </Show>
    </ScrollArea>
  );
}

function Group(props: { title: string; bindings: readonly InputBinding[] }) {
  return (
    <section class={classes.group} aria-label={props.title}>
      <h3 class={classes.heading}>{props.title}</h3>
      <dl class={classes.list}>
        {props.bindings.map((b) => (
          <div class={classes.row}>
            <dt class={classes.action}>{b.action}</dt>
            <dd class={classes.input}>
              <Keys input={b.input} />
              {b.or ? (
                <>
                  <span class={classes.or}>or</span>
                  <Keys input={b.or} />
                </>
              ) : null}
            </dd>
          </div>
        ))}
      </dl>
    </section>
  );
}

function Keys(props: { input: readonly string[] }) {
  return (
    <span class={classes.keys}>
      {props.input.map((k) => (
        <kbd class={classes.key}>{k}</kbd>
      ))}
    </span>
  );
}

/** Chords used by more than one command: [the keys, the commands' titles]. */
function conflictsOf(commands: readonly CommandDefinition<any>[]): [string[], string[]][] {
  const byChord = new Map<string, CommandDefinition<any>[]>();
  for (const c of commands) {
    const keys = chordKeys(c.shortcut!);
    const id = keys.join('+');
    byChord.set(id, [...(byChord.get(id) ?? []), c]);
  }
  return [...byChord.values()]
    .filter((list) => list.length > 1)
    .map((list) => [chordKeys(list[0].shortcut!), list.map((c) => c.title)]);
}
