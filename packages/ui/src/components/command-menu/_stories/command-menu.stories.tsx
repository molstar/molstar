/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { story } from '../../../_stories/support.js';
import type { CommandValues } from '../../../commands/types.js';
import { For, onCleanup, Show } from '../../../dom/index.js';
import { state } from '../../../state.js';
import { CommandButton } from '../../button/command-button.js';
import { installTooltips } from '../../tooltip/tooltip.js';
import { CommandMenu, type CommandMenuProps } from '../command-menu.js';
import { colorSwatches, fakeAssistant, prefixes, sampleManager, showAsForm } from './sample.js';

/**
 * The command menu over a sample command set (a fake context, no viewer). Its columns are the foundations' (styles §5):
 * the head's icon in the slot at the inset over the rows' icons; what is typed, the command's title and the rows' text
 * on one column; captions xs on it too; rows 24px.
 */
const meta: Meta = { title: 'Components/Command menu' };
export default meta;

const stage =
  'position: relative; height: 600px; margin: var(--msui-space-4); border-radius: var(--msui-radius-lg); background: var(--msui-surface); overflow: hidden;';
const logStyle =
  'position: absolute; left: var(--msui-space-4); bottom: var(--msui-space-3); margin: 0; color: var(--msui-text-muted); font-size: var(--msui-text-xs);';

interface Setup {
  loaded?: boolean;
  /** Opens this command's sentence first, with these values. */
  open?: [id: string, values?: CommandValues, tail?: string];
  menu?: Partial<Omit<CommandMenuProps, 'manager' | 'onClose'>>;
}

/** The menu over a stage; it reopens from a button once closed, and the stage logs what ran. */
function Stage(props: Setup) {
  const { manager, log } = sampleManager({ loaded: props.loaded ?? true });
  const open = state(true);
  const start = () => {
    if (props.open) manager.open(props.open[0], props.open[2], props.open[1]);
  };
  start();
  return (
    <div style={stage} ref={(el: HTMLElement) => onCleanup(installTooltips(el))}>
      <button
        type="button"
        class="msui-control"
        style="margin: var(--msui-space-3);"
        onClick={() => {
          start();
          open.set(true);
        }}
      >
        Open the command menu
      </button>
      <Show when={open}>
        {() => <CommandMenu manager={manager} onClose={() => open.set(false)} {...(props.menu ?? {})} />}
      </Show>
      <p style={logStyle}>
        <For each={log}>{(l) => <span style="display: block;">{l}</span>}</For>
      </p>
    </div>
  );
}

const make = (setup: Setup): StoryObj => ({ render: () => story(() => <Stage {...setup} />) });

/** Empty: Recent (once something ran), the groups (Browse), then the commands by category; → opens a group. */
export const Browsing = make({ menu: { prefixes } });

/** Typed: every command at any depth, the matched letters marked; a direct alias ("show surface") first. */
export const Search = make({ menu: { initial: 'col', prefixes } });

/** Inside a group with a shared argument: the mode as toggles; the highlighted row counted ("3 residues"). */
export const GroupWithSharedMode = make({ menu: { group: 'select' } });

/**
 * A sentence: the caret in the active field's token (a selection, with its completions and hint), a select, a number
 * with its presets and range, an optional field showing its default; the Options fold closed under the list.
 */
export const Sentence = make({ open: ['look.show'] });

/** The sentence with values set; the fold says what differs from the defaults. */
export const OptionsFold = make({
  open: ['look.show', { what: 'ligands', type: 'ball-and-stick', quality: 'high', opacity: 0.6 }],
});

/** A command's own form: the representations as tiles for its `type` field, under the sentence. */
export const OwnForm = make({ open: ['look.show', { what: 'A' }], menu: { forms: { 'look.show': showAsForm } } });

/** A field editor by format: every `color` field gets the swatch grid in place of the choices (←/→ move). */
export const FieldEditorByFormat = make({
  open: ['look.color', { what: 'ligands' }],
  menu: { fieldViews: { color: colorSwatches } },
});

/** `?` help: the prefixes (picking one types it), then the commands with keys. */
export const PrefixesAndHelp = make({ menu: { initial: '?', prefixes } });

/** A prefix's own list (`@` the scene). */
export const Prefix = make({ menu: { initial: '@', prefixes } });

/** Nothing loaded: rows that cannot run say why, greyed; picking one says it above the list. */
export const BlockedRows = make({ loaded: false, menu: { initial: 'o' } });

/** A run that fails keeps the sentence open with its error; ⌥↵ offers the command's place in the UI. */
export const Failure = make({ open: ['look.failing'] });

/** The Ask row: ↵ asks when no command matches; "paint the ligands red" proposes Color…, filled in. */
export const Assistant = make({ menu: { initial: 'paint the ligands red', assistant: fakeAssistant() } });

/** No model set up: a quiet way to the setup. */
export const AssistantNotSetUp = make({ menu: { assistant: fakeAssistant(false) } });

/** Opened for one command from elsewhere: Esc closes the menu. */
export const Direct = make({ open: ['file.download'], menu: { direct: true } });

/** Buttons by command id: the tooltip with the key; disabled with the reason; a label or the icon alone. */
export const CommandButtons: StoryObj = {
  render: () =>
    story(() => {
      const { manager, context } = sampleManager({ loaded: false });
      return (
        <div
          style="display: flex; gap: var(--msui-space-2); align-items: center; padding: var(--msui-space-5);"
          ref={(el: HTMLElement) => onCleanup(installTooltips(el))}
        >
          <CommandButton manager={manager} command="view.reset" label />
          <CommandButton manager={manager} command="edit.remove" label />
          <CommandButton manager={manager} command="file.open" />
          <CommandButton manager={manager} command="view.reset" size="sm" variant="ghost" />
          <button
            type="button"
            class="msui-control"
            onClick={() => context.set({ ...context.value, loaded: !context.value.loaded })}
          >
            Toggle loaded
          </button>
        </div>
      );
    }),
};
