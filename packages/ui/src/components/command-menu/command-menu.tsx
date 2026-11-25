/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { effect, onCleanup, rememberFocus, Show, trapTab } from '../../dom/index.js';
import { state } from '../../state.js';
import { classes } from './_generated/command-menu.module.css.js';
import { ListMode } from './list.js';
import { SentenceMode } from './sentence-mode.js';
import type { CommandMenuProps } from './types.js';

export type { CommandMenuItem, CommandMenuPrefix } from '../../commands/types.js';
export type {
  CommandFieldView,
  CommandFieldViewProps,
  CommandFormView,
  CommandMenuAssistant,
  CommandMenuProps,
  CommandViewProps,
} from './types.js';

/**
 * The command menu (commands §5): a modal dialog over its positioned container (a scrim fills it; the content behind
 * stays visible for previews). Empty, it browses (Recent, groups, categories); typed, it searches; a command with
 * params opens as its sentence, which is the editor: the caret in the active field's token, the field's choices under
 * it, the params past the sentence in an Options fold. A command that ran closes it; a failure keeps the sentence open
 * with its error. Closing gives the focus back to what opened it.
 */
export function CommandMenu(props: CommandMenuProps) {
  const { manager } = props;
  // Open, the menu wants every command: lazy ones are made now and follow changes while it is open (commands §4.2).
  onCleanup(manager.want());
  // The words the assistant proposed the open sentence from, until it closes.
  const proposed = state<string | undefined>(undefined);
  effect(manager.active, (a) => {
    if (!a) proposed.set(undefined);
  });
  effect(manager.outcomes, (o) => {
    if (o.kind === 'done') props.onClose();
  });
  onCleanup(() => manager.close());
  const restore = rememberFocus();
  let menu: HTMLElement | undefined;
  onCleanup(() => restore(menu));
  return (
    <div
      class={classes.scrim}
      onPointerDown={(e) => {
        if (e.target === e.currentTarget) props.onClose();
      }}
    >
      <div
        class={classes.menu}
        role="dialog"
        aria-modal="true"
        aria-label="Commands"
        ref={(el) => {
          menu = el;
        }}
        onKeyDown={(e) => {
          trapTab(e, e.currentTarget);
          // Escape anywhere in the menu not taken by what has the focus (a select's list) closes it.
          if (e.key !== 'Escape' || e.defaultPrevented) return;
          e.preventDefault();
          e.stopPropagation();
          props.onClose();
        }}
      >
        <Show when={manager.active} fallback={() => <ListMode {...props} onPropose={(t) => proposed.set(t)} />}>
          {(invocation) => (
            <SentenceMode
              manager={manager}
              invocation={invocation}
              proposed={proposed.asReadonly}
              onClose={props.onClose}
              {...(props.direct ? { direct: true } : {})}
              {...(props.forms?.[invocation.command.id] ? { form: props.forms[invocation.command.id] } : {})}
              {...(props.fieldViews ? { fieldViews: props.fieldViews } : {})}
            />
          )}
        </Show>
      </div>
    </div>
  );
}
