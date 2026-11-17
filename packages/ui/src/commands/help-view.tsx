/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { derive, KeyedSelector, type ReadonlyState } from '@molstar/core/reactive';
import { Button } from '../components/button/button.js';
import { Icon } from '../components/icon/icon.js';
import { Select } from '../components/select/select.js';
import { type Child, For, Show, timeout } from '../dom/index.js';
import { ViewHeading, ViewStack } from '../shell/view-frame.js';
import { classes } from './_generated/help-view.module.css.js';
import { chordFromEvent, chordKeys } from './keys.js';
import { type InputBinding, type InputBindingGroup, KeysModel } from './keys-model.js';
import type { CommandManager } from './manager.js';
import type { CommandDefinition } from './types.js';

export type { InputBinding, InputBindingGroup } from './keys-model.js';

/** Help's words (commands §6.6): strings, or functions for what names a command or a key. */
export interface HelpLabels {
  /** The conflicts' section (its accessible name). */
  readonly conflicts: string;
  /** After a conflict's keys: the commands on them and their scope (undefined: global). */
  readonly conflict: (commands: readonly string[], scope: string | undefined) => string;
  readonly everywhere: string;
  readonly commands: string;
  readonly changeKeys: string;
  /** The keyless commands' choice: its accessible name and its placeholder. */
  readonly giveKey: string;
  readonly giveKeyPlaceholder: string;
  readonly resetAll: string;
  readonly noKey: string;
  readonly or: string;
  readonly change: string;
  readonly changeTip: (command: string) => string;
  readonly clear: string;
  readonly clearTip: (command: string) => string;
  readonly reset: string;
  readonly resetTip: (command: string) => string;
  /** The recorder's accessible name. */
  readonly recorder: (command: string) => string;
  readonly recorderText: string;
  /** The default keys that work everywhere (an app's `global` replaces them). */
  readonly closeAction: string;
  readonly nextControl: string;
}

export const DefaultHelpLabels: HelpLabels = {
  conflicts: 'Shortcut conflicts',
  conflict: (commands, scope) =>
    ` is used by ${commands.join(' and ')}${scope ? ` (in ${scope})` : ''}: only the first that can run takes it.`,
  everywhere: 'Everywhere',
  commands: 'Commands',
  changeKeys: 'Change keys',
  giveKey: 'Give a key to a command',
  giveKeyPlaceholder: 'Give a key to…',
  resetAll: 'Reset all',
  noKey: 'no key',
  or: 'or',
  change: 'Change',
  changeTip: (c) => `Change the key of ${c}`,
  clear: 'Clear',
  clearTip: (c) => `Remove the key of ${c}`,
  reset: 'Reset',
  resetTip: (c) => `Back to the default key of ${c}`,
  recorder: (c) => `Press the new keys for ${c}; Escape cancels`,
  recorderText: 'Press keys…',
  closeAction: 'Close a menu, a card or a dialog',
  nextControl: 'Next control',
};

/** The keys that work everywhere unless an app says otherwise. */
const defaultGlobalKeys = (l: HelpLabels): readonly InputBinding[] => [
  { action: l.closeAction, input: ['Esc'] },
  { action: l.nextControl, input: ['Tab'], or: ['Shift', 'Tab'] },
];

/**
 * Help (commands §6.6): an app's prose (`intro`), shortcut conflicts first, the app's input bindings (mouse, touch), the
 * keys that work everywhere, then every command with a key after the keymap. "Change keys" records, clears or resets a
 * command's key, or gives one to a command without; a polite live region says what each did. The state and every
 * change are a `KeysModel`'s (DOM-free); this renders it, in a view's body (no scroll of its own, styles §8.4).
 */
export function HelpView(props: {
  commands: CommandManager<any>;
  bindings?: readonly InputBindingGroup[];
  /** The keys that work everywhere (the app's command menu key, Help's own); default Esc and Tab. */
  global?: readonly InputBinding[];
  intro?: () => Child;
  /** The model, when the host keeps it (an agent or a test drives it); else the view's own. */
  keys?: KeysModel<any>;
  labels?: Partial<HelpLabels>;
}) {
  const l = { ...DefaultHelpLabels, ...props.labels };
  const model = props.keys ?? new KeysModel(props.commands);
  const conflicts = derive(model.conflicts, (c) => (c.length ? c : undefined));
  const recording = new KeyedSelector(model.recording);
  // Controlled: the command being recorded while it is keyless, else the placeholder again (after a key is taken or
  // the recording is cancelled).
  const given = derive([model.recording, model.keyless], (id, keyless) =>
    id && keyless.some((c) => c.id === id) ? id : '',
  );
  return (
    <ViewStack class={classes.help}>
      {props.intro?.()}
      <Show when={conflicts}>
        {(list) => (
          <section class={classes.conflicts} aria-label={l.conflicts}>
            <Icon name="triangle-alert" />
            <div>
              {list.map((x) => (
                <p class={classes.conflict}>
                  <Keys input={chordKeys(x.chord)} />
                  {l.conflict(
                    x.commands.map((c) => c.title),
                    x.scope === 'global' ? undefined : x.scope,
                  )}
                </p>
              ))}
            </div>
          </section>
        )}
      </Show>
      {(props.bindings ?? []).map((group) => (
        <Group title={group.title} bindings={group.bindings} or={l.or} />
      ))}
      <Group title={l.everywhere} bindings={props.global ?? defaultGlobalKeys(l)} or={l.or} />
      <section class={classes.group} aria-label={l.commands}>
        <div class={classes.headingRow}>
          <ViewHeading>{l.commands}</ViewHeading>
          <Button
            size="xs"
            variant="ghost"
            pressed={model.editing}
            onClick={() => model.setEditing(!model.editing.value)}
          >
            {l.changeKeys}
          </Button>
        </div>
        <dl class={classes.list}>
          <For each={model.listed} by={(c) => c.id}>
            {(c) => <CommandKeys model={model} command={c} recording={recording.key(c.id)} labels={l} />}
          </For>
        </dl>
        <Show when={derive(model.editing, (e) => (e ? true : undefined))}>
          {() => (
            <div class={classes.editBar}>
              <Show when={derive(model.keyless, (k) => (k.length ? k : undefined))}>
                {(keyless) => (
                  <Select
                    size="sm"
                    aria-label={l.giveKey}
                    options={[['', l.giveKeyPlaceholder], ...keyless.map((c) => [c.id, c.title] as const)]}
                    value={given}
                    onChange={(id) => id && model.record(id)}
                  />
                )}
              </Show>
              <Button size="xs" variant="ghost" onClick={() => model.resetAll()}>
                {l.resetAll}
              </Button>
            </div>
          )}
        </Show>
        <span class={classes.srOnly} role="status" aria-live="polite">
          {model.said}
        </span>
      </section>
    </ViewStack>
  );
}

/** A command's row: its keys; while changing keys, Change (records the next keys), Clear and Reset. */
function CommandKeys(props: {
  model: KeysModel<any>;
  command: CommandDefinition<any>;
  recording: ReadonlyState<boolean>;
  labels: HelpLabels;
}) {
  const { model, command: c, labels: l } = props;
  const keymap = model.manager.keymap;
  const state = derive(keymap.layers, () => ({ chord: model.manager.shortcut(c), remapped: keymap.isRemapped(c.id) }));
  const chord = derive(state, (s) => s.chord);
  return (
    <div class={`msui-row msui-row-plain ${classes.row}`}>
      <dt class={classes.action}>{c.title}</dt>
      <dd class={classes.input}>
        <Show
          when={derive(props.recording, (r) => (r ? true : undefined))}
          fallback={() => (
            <Show when={chord} fallback={() => <span class={classes.or}>{l.noKey}</span>}>
              {(k) => <Keys input={chordKeys(k)} />}
            </Show>
          )}
        >
          {() => <Recorder model={model} command={c} labels={l} />}
        </Show>
        <Show when={derive(model.editing, (e) => (e ? true : undefined))}>
          {() => (
            <span class={classes.edit}>
              <Button size="xs" variant="ghost" title={l.changeTip(c.title)} onClick={() => model.record(c.id)}>
                {l.change}
              </Button>
              <Show when={derive(state, (s) => (s.chord ? true : undefined))}>
                {() => (
                  <Button size="xs" variant="ghost" title={l.clearTip(c.title)} onClick={() => model.clear(c.id)}>
                    {l.clear}
                  </Button>
                )}
              </Show>
              <Show when={derive(state, (s) => (s.remapped ? true : undefined))}>
                {() => (
                  <Button size="xs" variant="ghost" title={l.resetTip(c.title)} onClick={() => model.reset(c.id)}>
                    {l.reset}
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

/**
 * Waiting for the new keys: a key-shaped well that takes the focus and every key but Tab (no chord: it moves the focus
 * on, and leaving cancels); Esc cancels too.
 */
function Recorder(props: { model: KeysModel<any>; command: CommandDefinition<any>; labels: HelpLabels }) {
  const { model, command: c, labels: l } = props;
  const el = (
    <span
      class={classes.recorder}
      tabIndex={0}
      role="button"
      aria-label={l.recorder(c.title)}
      onBlur={() => model.cancel()}
      onKeyDown={(e: KeyboardEvent) => {
        // Tab is no chord (`chordFromEvent`): it leaves, and the blur cancels; a keyboard user is never trapped.
        if (e.key === 'Tab') return;
        // Nothing behind the recorder takes these keys (an app's shortcuts, a menu's Escape).
        e.preventDefault();
        e.stopPropagation();
        if (e.key === 'Escape') return model.cancel();
        const chord = chordFromEvent(e);
        if (chord) model.recorded(chord, c.id);
      }}
    >
      {l.recorderText}
    </span>
  ) as HTMLElement;
  // After whatever opened it settles its own focus (a closing choice gives it back to its button).
  timeout(() => el.isConnected && el.focus(), 0);
  return el;
}

function Group(props: { title: string; bindings: readonly InputBinding[]; or: string }) {
  return (
    <section class={classes.group} aria-label={props.title}>
      <ViewHeading>{props.title}</ViewHeading>
      <dl class={classes.list}>
        {props.bindings.map((b) => (
          <div class={`msui-row msui-row-plain ${classes.row}`}>
            <dt class={classes.action}>{b.action}</dt>
            <dd class={classes.input}>
              <Keys input={b.input} />
              {b.or ? (
                <>
                  <span class={classes.or}>{props.or}</span>
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
    <kbd class={classes.keys}>
      {props.input.map((k) => (
        <kbd class={classes.key}>{k}</kbd>
      ))}
    </kbd>
  );
}
