/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// An open command (commands §5.3–5.5): the sentence it is edited in, the active field's choices (or its editor), the
// command's own form, the Options fold, the keys. The input keeps the focus and the keys.

import { derive, type ReadonlyState } from '@molstar/core/reactive';
import type { CommandInvocation } from '../../commands/invocation.js';
import type { CommandManager } from '../../commands/manager.js';
import type { FieldChoice } from '../../commands/sentence.js';
import { deepActiveElement, focusables, onCleanup, Show, setAttr, uniqueId, watch } from '../../dom/index.js';
import { Icon } from '../icon/icon.js';
import { classes as menu } from './_generated/command-menu.module.css.js';
import { classes } from './_generated/sentence.module.css.js';
import { ChoiceList, FieldLine, RangeRow } from './choices.js';
import { Footer } from './footer.js';
import { OptionsFold } from './options-fold.js';
import { SentenceHead, tokenText } from './sentence.js';
import type { CommandFieldView, CommandFormView, CommandViewProps } from './types.js';

const NoChoices: readonly FieldChoice[] = [];

export function SentenceMode(props: {
  manager: CommandManager<any>;
  invocation: CommandInvocation<any>;
  /** The words the assistant proposed it from. */
  proposed?: ReadonlyState<string | undefined>;
  onClose(): void;
  /** Opened for this command alone: leaving it closes the menu. */
  direct?: boolean;
  form?: CommandFormView;
  fieldViews?: Readonly<Record<string, CommandFieldView>>;
}) {
  const { manager, invocation } = props;
  const state = invocation.state;
  const command = invocation.command;
  const listId = uniqueId('msui-command-choices');
  const optionId = (i: number) => `${listId}-${i}`;
  // Keys seen before the menu's: the own form's, and the active field editor's.
  let formKeys: ((e: KeyboardEvent) => boolean) | undefined;
  let fieldKeys: ((e: KeyboardEvent) => boolean) | undefined;
  // Who draws the active field's choices: the own form (a field it names), a field editor (by format), or the list.
  const drawer = derive(state, (s): string | undefined => {
    if (s.field !== undefined && props.form?.fields?.includes(s.field)) return 'form';
    const format = s.fields.find((f) => f.id === s.field)?.format;
    return format && props.fieldViews?.[format] ? `field:${s.field}:${format}` : undefined;
  });
  const items = derive([drawer, invocation.choices], (d, c) => (d || !c ? NoChoices : c.items));

  // What is below the sentence and takes the focus by Tab: the own form's controls, the Options fold.
  let below: HTMLElement | undefined;
  const belowStops = () => (below ? focusables(below).filter((el) => el.tabIndex >= 0) : []);
  /** Tab past the last field (Shift+Tab before the first), nothing to take: into what is below. False: the sentence's. */
  const intoBelow = (back: boolean): boolean => {
    const s = state.value;
    if (s.canAssist || (s.text.trim() && items.value.length) || s.pick >= 0) return false;
    const i = s.fields.findIndex((f) => f.id === s.field);
    if (i >= 0 && (back ? i > 0 : i < s.fields.length - 1)) return false;
    const stops = belowStops();
    const to = back ? stops[stops.length - 1] : stops[0];
    if (!to) return false;
    to.focus();
    return true;
  };
  /** Tab past the last control below (Shift+Tab before the first): the input, on the first (last) field. */
  const outOfBelow = (e: KeyboardEvent) => {
    if (e.key !== 'Tab' || e.defaultPrevented || e.altKey || e.metaKey || e.ctrlKey || !below) return;
    const stops = belowStops();
    const at = deepActiveElement(below.ownerDocument);
    if (at !== (e.shiftKey ? stops[0] : stops[stops.length - 1])) return;
    e.preventDefault();
    const fields = state.value.fields;
    const f = e.shiftKey ? fields[fields.length - 1] : fields[0];
    if (f) invocation.focus(f.id);
    input.focus();
  };
  const leave = () => {
    manager.close();
    if (props.direct) props.onClose();
  };
  const goToView = () => {
    const view = command.view;
    if (!view) return;
    manager.close();
    props.onClose();
    view.open(manager.context.value);
  };
  const run = () => void invocation.enter();
  const focusField = (id: string) => {
    invocation.focus(id);
    input.focus();
  };

  const input = (
    <input
      class={classes.input}
      type="text"
      role="combobox"
      aria-controls={listId}
      aria-autocomplete="list"
      autocomplete="off"
      spellcheck={false}
      size={1}
      onInput={(e) => invocation.type(e.currentTarget.value)}
      onKeyDown={(e) => {
        if (formKeys?.(e) || fieldKeys?.(e)) {
          e.preventDefault();
          return;
        }
        const el = e.currentTarget;
        const s = state.value;
        const i = s.fields.findIndex((f) => f.id === s.field);
        if (e.key === 'Enter' && e.altKey && command.view) {
          e.preventDefault();
          goToView();
        } else if (e.key === 'Enter') {
          e.preventDefault();
          run();
        } else if (e.key === 'Tab') {
          e.preventDefault();
          if (!intoBelow(e.shiftKey)) void invocation.tab(e.shiftKey ? -1 : 1);
        } else if (e.key === 'ArrowDown' || e.key === 'ArrowUp') {
          e.preventDefault();
          invocation.move(e.key === 'ArrowDown' ? 1 : -1);
        } else if (e.key === 'ArrowLeft' && atEdge(el, 'start') && i > 0) {
          // At the text's start: the field before (commands §5.3).
          e.preventDefault();
          invocation.step(-1);
        } else if (e.key === 'ArrowRight' && atEdge(el, 'end') && i >= 0 && i < s.fields.length - 1) {
          e.preventDefault();
          invocation.step(1);
        } else if (e.key === 'Backspace' && !el.value) {
          e.preventDefault();
          if (!invocation.back()) leave();
        } else if (e.key === 'Escape') {
          e.preventDefault();
          e.stopPropagation();
          leave();
        }
      }}
    />
  ) as HTMLInputElement;
  // The input's box: the input over a sizer that takes its text (or what the token reads as), so it is as wide.
  const edit = (<span class={classes.edit}>{input}</span>) as HTMLElement;
  queueMicrotask(() => input.focus());

  // The input as a function of the state and the list: one watch.
  watch([state, items, drawer], () => {
    const s = state.value;
    const f = invocation.sentence.value.find((x) => x.id === s.field);
    const listed = items.value.length > 0;
    const shown = f ? tokenText(f) : undefined;
    setAttr(input, 'aria-expanded', String(listed));
    setAttr(
      input,
      'aria-activedescendant',
      listed && s.pick >= 0 && s.pick < items.value.length ? optionId(s.pick) : undefined,
    );
    setAttr(input, 'aria-label', f ? f.label : command.title);
    setAttr(input, 'placeholder', shown ? shown.text : `${s.action}: Enter`);
    setAttr(edit, 'data-state', shown?.state);
    setAttr(edit, 'data-text', s.text || (shown ? shown.text : `${s.action}: Enter`));
    setAttr(input, 'readonly', s.status !== 'editing');
    if (input.value !== s.text) input.value = s.text;
  });

  const viewProps: CommandViewProps = {
    manager,
    invocation,
    run,
    keys: (handler) => {
      formKeys = handler;
    },
  };
  const fieldView = derive(drawer, (d) => (d?.startsWith('field:') ? d : undefined));

  return (
    <>
      <Show when={props.proposed}>
        {(words) => (
          <p class={classes.proposed}>
            <span class={menu.headSlot}>
              <Icon name="sparkles" />
            </span>
            Proposed from “{words}”: check it, then Enter to run
          </p>
        )}
      </Show>
      <SentenceHead
        invocation={invocation}
        edit={edit}
        input={input}
        onFocusField={focusField}
        onRun={run}
        {...(command.view ? { onView: goToView } : {})}
      />
      <FieldLine invocation={invocation} choices={invocation.choices} />
      <RangeRow invocation={invocation} />
      <ChoiceList invocation={invocation} items={items} id={listId} optionId={optionId} onTaken={() => input.focus()} />
      <Show when={fieldView}>
        {(key) => {
          const [, field, format] = key.split(':') as [string, string, string];
          const view = props.fieldViews![format]!;
          onCleanup(() => {
            fieldKeys = undefined;
          });
          return (
            <div class={classes.fieldView} onMouseDown={(e) => e.preventDefault()}>
              {view.render({
                manager,
                invocation,
                run,
                field,
                choices: invocation.choices,
                set: (value) => {
                  invocation.set(field, value);
                  input.focus();
                },
                keys: (handler) => {
                  fieldKeys = handler;
                },
              })}
            </div>
          );
        }}
      </Show>
      <div
        class={classes.below}
        ref={(el: HTMLElement) => {
          below = el;
        }}
        onKeyDown={outOfBelow}
      >
        {props.form ? <div class={classes.form}>{props.form.render(viewProps)}</div> : null}
        {props.form?.options === false ? null : <OptionsFold invocation={invocation} />}
      </div>
      <Footer
        hints={[
          ['↑↓', 'choose'],
          ['tab', 'next field'],
          ['↵', 'run'],
          ...(command.view ? [['⌥↵', `open ${command.view.title}`] as const] : []),
          ['⌫', 'back'],
          props.direct ? ['esc', 'close'] : ['esc', 'commands'],
        ]}
      />
    </>
  );
}

/** Is the caret at the text's start (end), nothing selected? */
function atEdge(el: HTMLInputElement, edge: 'start' | 'end'): boolean {
  const { selectionStart: a, selectionEnd: b } = el;
  if (a === null || a !== b) return false;
  return edge === 'start' ? a === 0 : a === el.value.length;
}
