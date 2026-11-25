/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The sentence is the editor (commands §5.3): the command, then its fields as tokens, the caret inside the active
// field's token. The one input (the combobox) moves into the token being edited; the tokens around it are buttons.

import { derive } from '@molstar/core/reactive';
import { sameShallow } from '@molstar/core/util/object';
import type { CommandInvocation } from '../../commands/invocation.js';
import type { SentenceField } from '../../commands/sentence.js';
import { deepActiveElement, For, setAttr, setShown, setText, watch } from '../../dom/index.js';
import type { FormFieldView } from '../../forms/flat-fields.js';
import { Button } from '../button/button.js';
import { Icon } from '../icon/icon.js';
import { classes as menu } from './_generated/command-menu.module.css.js';
import { classes } from './_generated/sentence.module.css.js';
import { Keys } from './footer.js';

/** How an empty or set field reads: its value, its name (required, empty) or its default (optional, empty). */
export function tokenText(f: SentenceField): { text: string; state: 'set' | 'missing' | 'default' } {
  if (!f.empty) return { text: f.display, state: 'set' };
  if (f.required) return { text: f.label, state: 'missing' };
  return { text: f.defaultDisplay || f.label, state: 'default' };
}

/** Moves `el` into `parent` (if not there yet), keeping the focus it had: a move out of the document blurs. */
function moveInto(parent: HTMLElement, el: HTMLElement, focused: HTMLElement) {
  if (el.parentNode === parent) return;
  const had = deepActiveElement(focused.ownerDocument) === focused;
  parent.append(el);
  if (had) focused.focus({ preventScroll: true });
}

/**
 * The head of an open command: its icon in the slot (the list's icon column), its title, the tokens, then where it
 * lives in the UI (⌥↵) and the run button (↵). `edit` is the input's box (the input and its sizer), placed in the active
 * token, or after the tokens when no field is active.
 */
export function SentenceHead(props: {
  invocation: CommandInvocation<any>;
  edit: HTMLElement;
  input: HTMLInputElement;
  onFocusField(id: string): void;
  onRun(): void;
  onView?: () => void;
}) {
  const { invocation, edit, input } = props;
  const command = invocation.command;
  const fields = derive(invocation.state, (s) => s.fields);
  const tail = (<span class={classes.tail} />) as HTMLElement;
  // No active field (every param is past the sentence, or it has none): the input waits after the tokens.
  watch([derive(invocation.state, (s) => s.fields.some((f) => f.id === s.field))], () => {
    const s = invocation.state.value;
    if (!s.fields.some((f) => f.id === s.field)) moveInto(tail, edit, input);
  });
  return (
    <div class={menu.head} data-sentence="">
      <span class={menu.headSlot} data-command="">
        {command.icon ? <Icon name={command.icon} /> : null}
      </span>
      <div class={classes.sentence} role="group" aria-label={`${command.title.replace(/…$/, '')}: arguments`}>
        <span class={classes.command}>{command.title.replace(/…$/, '')}</span>
        <For each={fields} by={(f) => f.id}>
          {(f) => <Token field={f} invocation={invocation} edit={edit} input={input} onFocus={props.onFocusField} />}
        </For>
        {tail}
        <span class={classes.end}>
          {props.onView && command.view ? (
            <Button variant="ghost" title={`Open ${command.view.title} (⌥↵)`} onClick={() => props.onView?.()}>
              {command.view.title}
            </Button>
          ) : null}
          <Button variant="primary" class={classes.run} onClick={() => props.onRun()}>
            {invocation.state.value.action}
            <Keys keys="↵" class={classes.runKeys} />
          </Button>
        </span>
      </div>
    </div>
  );
}

/**
 * One field: its quiet before-word, then its token. Not active, a button showing its value (or its name, dashed, when
 * required and empty; its default, muted, when optional and empty); active, the input inside it.
 */
function Token(props: {
  field: FormFieldView;
  invocation: CommandInvocation<any>;
  edit: HTMLElement;
  input: HTMLInputElement;
  onFocus(id: string): void;
}) {
  const { field, invocation, edit, input } = props;
  const view = derive(invocation.sentence, (all) => all.find((x) => x.id === field.id), sameShallow);
  const text = document.createTextNode('');
  const button = (
    <button
      type="button"
      class={classes.value}
      tabIndex={-1}
      onMouseDown={(e) => e.preventDefault()}
      onClick={() => props.onFocus(field.id)}
    >
      {text}
    </button>
  ) as HTMLElement;
  const slot = (<span class={classes.slot} />) as HTMLElement;
  const box = (
    <span class={classes.token}>
      {button}
      {slot}
    </span>
  ) as HTMLElement;
  watch([view], () => {
    const f = view.value;
    if (!f) return;
    const t = tokenText(f);
    setAttr(box, 'data-state', t.state);
    setAttr(box, 'data-active', f.active);
    setAttr(box, 'data-error', !!f.error);
    setAttr(box, 'data-proposed', f.proposed);
    setText(text, t.text);
    setAttr(button, 'aria-label', `${f.label}: ${t.state === 'missing' ? 'empty' : t.text}`);
    setShown(button, !f.active);
    if (f.active) moveInto(slot, edit, input);
  });
  return (
    <span class={classes.field} data-field={field.id}>
      {field.before ? <span class={classes.before}>{field.before}</span> : null}
      {box}
    </span>
  );
}
