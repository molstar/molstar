/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { chordFromEvent, chordKeys, formatChord, isMac } from '../../commands/keys.js';
import type { CommandManager } from '../../commands/manager.js';
import type { CommandDefinition, KeyChord } from '../../commands/types.js';
import { Button } from '../../components/button/button.js';
import { Icon } from '../../components/icon/icon.js';
import { ScrollArea } from '../../components/scroll-area/scroll-area.js';
import { Select } from '../../components/select/select.js';
import { For, Show } from '../../dom/control-flow.js';
import type { Child } from '../../dom/jsx-types.js';
import { type KeyedValue, keyedSelector } from '../../reactive/keyed-selector.js';
import { observableState } from '../../reactive/observable-state.js';
import type { ReadonlyState } from '../../reactive/readonly-state.js';
import { derive } from '../../reactive/state.js';
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

/** The keys that work everywhere, whatever the mode: in this platform's keys (⌘ on macOS, Ctrl elsewhere). */
const globalKeys = (): InputBinding[] => [
  {
    action: 'Command menu',
    input: chordKeys({ key: 'k', mod: true }),
    or: chordKeys({ key: 'p', mod: true, shift: true }),
  },
  { action: 'Help (this card)', input: ['?'] },
  { action: 'Close a menu, a card or the command menu', input: ['Esc'] },
];

/**
 * Help: how to work the canvas (the mode's mouse and touch bindings), the keys (global ones, then every command with a
 * key, after the keymap: the app's and the person's), and shortcut conflicts (two commands on one key in one scope). A
 * tab at the rail's foot (it stays open while you try things on the canvas); also the command menu's `?`. "Change keys"
 * lets the person record, clear or reset a command's key, or give one to a command without (commands.md §6).
 */
export function HelpPanel(props: {
  commands: CommandManager<any>;
  bindings: readonly InputBindingGroup[];
  /** Above the bindings: a mode's prose (the viewer's basics). */
  intro?: () => Child;
}) {
  const manager = props.commands;
  const keymap = manager.keymap;
  const editing = observableState(false);
  /** The command whose new key is being recorded. */
  const recording = observableState<string | undefined>(undefined);
  const recordingKey = keyedSelector(recording.asReadonly);
  /** What the last change did, said once (a polite live region). */
  const said = observableState('');
  // What is listed: commands with a key, those the person remapped or cleared (to reset them), the one being recorded.
  const listed = derive([manager.commands, keymap.layers, recording.asReadonly], (all, _, rec) =>
    all.filter((c) => manager.shortcut(c) || keymap.isRemapped(c.id) || c.id === rec),
  );
  const conflicts = derive([manager.commands, keymap.layers], (all) => {
    const c = keymap.conflicts(all);
    return c.length ? c : undefined;
  });
  const record = (c: CommandDefinition<any>, chord: KeyChord) => {
    keymap.remap(c.id, chord);
    recording.set(undefined);
    const others = keymap.sharing(c.id, chord, manager.commands.value);
    said.set(
      `${c.title}: ${formatChord(chord)}${others.length ? `, also used by ${others.map((o) => o.title).join(' and ')}` : ''}.`,
    );
  };
  return (
    <ScrollArea edge class={classes.help} contentClass={classes.content} aria-label="Help">
      {props.intro?.()}
      <Show when={conflicts}>
        {(list) => (
          <section class={classes.conflicts} aria-label="Shortcut conflicts">
            <Icon name="triangle-alert" size={14} />
            <div>
              {list.map((x) => (
                <p class={classes.conflict}>
                  <Keys input={chordKeys(x.chord)} /> is used by {x.commands.map((c) => c.title).join(' and ')}
                  {x.scope === 'global' ? '' : ` (in ${x.scope})`}: only the first that can run takes it.
                </p>
              ))}
            </div>
          </section>
        )}
      </Show>
      <For each={props.bindings}>{(group) => <Group title={group.title} bindings={group.bindings} />}</For>
      <Group title="Everywhere" bindings={globalKeys()} />
      <section class={classes.group} aria-label="Commands">
        <div class={classes.headingRow}>
          <h3 class={classes.heading}>Commands</h3>
          <Button
            size="xs"
            variant="ghost"
            pressed={editing.asReadonly}
            onClick={() => {
              editing.set(!editing.value);
              recording.set(undefined);
            }}
          >
            Change keys
          </Button>
        </div>
        <dl class={classes.list}>
          <For each={listed} key={(c) => c.id}>
            {(c) => (
              <CommandKeys
                manager={manager}
                command={c}
                editing={editing.asReadonly}
                recording={recordingKey.key(c.id)}
                onRecord={() => recording.set(c.id)}
                onRecorded={(chord) => (chord ? record(c, chord) : recording.set(undefined))}
                onChanged={(text) => said.set(text)}
              />
            )}
          </For>
        </dl>
        <Show when={derive(editing.asReadonly, (e) => (e ? true : undefined))}>
          {() => (
            <div class={classes.editBar}>
              <Select
                size="sm"
                aria-label="Give a key to a command"
                options={[
                  ['', 'Give a key to…'],
                  ...manager.commands.value
                    .filter((c) => !manager.shortcut(c))
                    .map((c) => [c.id, c.title] as const)
                    .sort((a, b) => a[1].localeCompare(b[1])),
                ]}
                value=""
                onChange={(id) => id && recording.set(id)}
              />
              <Button
                size="xs"
                variant="ghost"
                onClick={() => (keymap.resetAll(), said.set('Every key is back to its default.'))}
              >
                Reset all
              </Button>
            </div>
          )}
        </Show>
        <span class={classes.srOnly} role="status" aria-live="polite">
          {said.asReadonly}
        </span>
      </section>
    </ScrollArea>
  );
}

/** A command's row: its keys; while editing, Change (records the next keys), Clear and Reset. */
function CommandKeys(props: {
  manager: CommandManager<any>;
  command: CommandDefinition<any>;
  editing: ReadonlyState<boolean>;
  recording: KeyedValue<boolean>;
  onRecord(): void;
  /** The keys pressed, or none when recording was cancelled. */
  onRecorded(chord: KeyChord | undefined): void;
  onChanged(text: string): void;
}) {
  const { manager, command: c } = props;
  const keymap = manager.keymap;
  const state = derive(keymap.layers, () => ({ chord: manager.shortcut(c), remapped: keymap.isRemapped(c.id) }));
  const mac = isMac();
  return (
    <div class={classes.row}>
      <dt class={classes.action}>{c.title}</dt>
      <dd class={classes.input}>
        <Show
          when={derive(props.recording, (r) => (r ? true : undefined))}
          fallback={() => (
            <Show when={derive(state, (s) => s.chord)} fallback={() => <span class={classes.or}>no key</span>}>
              {(chord) => <Keys input={chordKeys(chord, mac)} />}
            </Show>
          )}
        >
          {() => (
            <span
              class={classes.recorder}
              tabIndex={0}
              role="button"
              aria-label={`Press the new keys for ${c.title}; Escape cancels`}
              // After whatever opened it settles its own focus (a closing choice gives it back to its button).
              ref={(el) => setTimeout(() => el.isConnected && el.focus())}
              onBlur={() => props.onRecorded(undefined)}
              onKeyDown={(e) => {
                // Nothing behind the recorder takes these keys (the viewer's shortcuts, the menu's Escape).
                e.preventDefault();
                e.stopPropagation();
                if (e.key === 'Escape') return props.onRecorded(undefined);
                const chord = chordFromEvent(e, mac);
                if (chord) props.onRecorded(chord);
              }}
            >
              Press keys…
            </span>
          )}
        </Show>
        <Show when={derive(props.editing, (e) => (e ? true : undefined))}>
          {() => (
            <span class={classes.edit}>
              <Button size="xs" variant="ghost" title={`Change the key of ${c.title}`} onClick={() => props.onRecord()}>
                Change
              </Button>
              <Show when={derive(state, (s) => (s.chord ? true : undefined))}>
                {() => (
                  <Button
                    size="xs"
                    variant="ghost"
                    title={`Remove the key of ${c.title}`}
                    onClick={() => {
                      keymap.remap(c.id, null);
                      props.onChanged(`${c.title}: no key.`);
                    }}
                  >
                    Clear
                  </Button>
                )}
              </Show>
              <Show when={derive(state, (s) => (s.remapped ? true : undefined))}>
                {() => (
                  <Button
                    size="xs"
                    variant="ghost"
                    title={`Back to the default key of ${c.title}`}
                    onClick={() => {
                      keymap.reset(c.id);
                      const d = keymap.defaultShortcut(c);
                      props.onChanged(`${c.title}: ${d ? formatChord(d, mac) : 'no key'} again.`);
                    }}
                  >
                    Reset
                  </Button>
                )}
              </Show>
            </span>
          )}
        </Show>
      </dd>
    </div>
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
