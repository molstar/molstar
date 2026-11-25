/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// What the active field takes (commands §5.3): the line above the list (an error, a hint, the assistant's state), a
// slider's range, and the choices as options of a listbox. Always the active field's.

import { derive, type ReadonlyState } from '@molstar/core/reactive';
import type { CommandInvocation } from '../../commands/invocation.js';
import type { FieldChoice, FieldChoices } from '../../commands/sentence.js';
import { For, Show, setAttr, setShown, setText, watch } from '../../dom/index.js';
import { Icon } from '../icon/icon.js';
import { ScrollArea } from '../scroll-area/scroll-area.js';
import { classes as menu } from './_generated/command-menu.module.css.js';
import { classes } from './_generated/sentence.module.css.js';

/**
 * One line above the list: the field's error (an alert), else what the assistant is doing or could do, its note, or a
 * format's hint. Hidden when there is nothing to say.
 */
export function FieldLine(props: {
  invocation: CommandInvocation<any>;
  choices: ReadonlyState<FieldChoices | undefined>;
}) {
  const { invocation } = props;
  const s = invocation.state;
  const error = derive(s, (x) => x.error);
  const info = derive([s, props.choices], (x, c) => {
    if (x.error) return '';
    const words = x.text.trim();
    if (x.status === 'resolving') return `Reading “${words}”…`;
    if (x.canAssist) return `Tab: the assistant reads “${words}”`;
    if (x.note) return x.note;
    return c?.hint ?? '';
  });
  const text = document.createTextNode('');
  const line = (<p class={classes.line}>{text}</p>) as HTMLElement;
  watch([info], () => {
    setText(text, info.value);
    setShown(line, !!info.value);
  });
  return (
    <>
      <Show when={error}>
        {(e) => (
          <p class={classes.line} data-tone="error" role="alert">
            {e.message}
            {e.hint ? <span class={classes.lineHint}> {e.hint}</span> : null}
          </p>
        )}
      </Show>
      {line}
    </>
  );
}

/** A slider's range (a number with min and max): drag sets the value and keeps the field active. */
export function RangeRow(props: { invocation: CommandInvocation<any> }) {
  const { invocation } = props;
  // Keyed by the layer's range (the same object while the field and its bounds are): a value set while dragging keeps
  // the input under the pointer.
  return (
    <Show when={invocation.range}>
      {(r) => {
        const field = invocation.sentence.value.find((f) => f.id === r.field);
        const input = (
          <input
            class={classes.range}
            type="range"
            min={String(r.min)}
            max={String(r.max)}
            step={r.step === undefined ? 'any' : String(r.step)}
            aria-label={field?.label ?? r.field}
            onInput={(e) => invocation.set(r.field, Number(e.currentTarget.value), { stay: true })}
          />
        ) as HTMLInputElement;
        const shown = document.createTextNode('');
        const value = derive(invocation.sentence, (all) => {
          const f = all.find((x) => x.id === r.field);
          return f ? (f.value ?? f.defaultValue) : undefined;
        });
        watch([value], () => {
          const v = typeof value.value === 'number' ? value.value : r.min;
          if (input.valueAsNumber !== v) input.value = String(v);
          setText(shown, String(v));
        });
        return (
          <div class={classes.rangeRow}>
            <span class={`${classes.rangeEnd} msui-num`}>{String(r.min)}</span>
            {input}
            <span class={`${classes.rangeEnd} msui-num`}>{String(r.max)}</span>
            <span class={`${classes.rangeValue} msui-num`}>{shown}</span>
          </div>
        );
      }}
    </Show>
  );
}

/**
 * The active field's choices as a listbox (the input owns it: `aria-controls`, `aria-activedescendant`): the current
 * value checked in the icon slot, the label, its detail. A click takes it; the pointer highlights it. A press keeps the
 * focus in the input.
 */
export function ChoiceList(props: {
  invocation: CommandInvocation<any>;
  items: ReadonlyState<readonly FieldChoice[]>;
  id: string;
  optionId(i: number): string;
  onTaken(): void;
}) {
  const { invocation } = props;
  const pick = derive(invocation.state, (s) => s.pick);
  return (
    <ScrollArea
      class={menu.list}
      contentClass={menu.listContent}
      id={props.id}
      role="listbox"
      aria-label="Choices"
      onPointerDown={(e: PointerEvent) => e.preventDefault()}
      onMouseDown={(e: MouseEvent) => e.preventDefault()}
    >
      <For each={props.items}>
        {(c, i) => {
          const option = (
            <div
              class={menu.row}
              role="option"
              onClick={() => {
                invocation.accept(i.value);
                props.onTaken();
              }}
              onPointerMove={() => invocation.point(i.value)}
            >
              <span class={menu.rowSlot}>{c.current ? <Icon name="check" /> : null}</span>
              <span class={menu.rowTitle} data-current={c.current ? '' : undefined}>
                {c.label}
              </span>
              <span class={menu.rowMeta}>{c.detail ?? ''}</span>
              <span />
            </div>
          ) as HTMLElement;
          // Its place and whether it is the pick: only a move of the pick or of the row wakes it, not a keystroke.
          watch([pick, i], () => {
            setAttr(option, 'id', props.optionId(i.value));
            setAttr(option, 'aria-selected', String(pick.value === i.value));
            if (pick.value === i.value) queueMicrotask(() => option.scrollIntoView?.({ block: 'nearest' }));
          });
          return option;
        }}
      </For>
    </ScrollArea>
  );
}
