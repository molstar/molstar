/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { Show } from '../../dom/control-flow.js';
import { watch } from '../../dom/scope.js';
import { setAttr } from '../../dom/writes.js';
import type { ObservableState } from '../../reactive/observable-state.js';
import type { ReadonlyState } from '../../reactive/readonly-state.js';
import { derive } from '../../reactive/state.js';
import { Icon } from '../icon/icon.js';
import { Floating, OptionList } from '../select/select.js';
import { classes } from './_generated/params-form.module.css.js';

/**
 * PD's presets (a group's `presets`, an object list's `presets`, a line graph's `getPresets`): named values to start
 * from, as Mol*'s Presets menu. A quiet button opens them as a popover; a pick applies one.
 */
export function PresetsButton(props: { open: ObservableState<boolean>; label: string }) {
  const { open } = props;
  const el = (
    <button
      type="button"
      class={classes.presets}
      ref={(el: HTMLButtonElement) => anchors.set(open, el)}
      aria-label={`${props.label} presets`}
      // The list's opener: a press on it toggles the list (it is not a press outside that closes it first).
      aria-haspopup="listbox"
      onClick={(e) => {
        // In a group's header: the click is the presets', not the header's (it would fold the group).
        e.stopPropagation();
        open.set(!open.value);
      }}
    >
      <Icon name="swatch-book" size={12} />
      Presets
    </button>
  ) as HTMLElement;
  watch([open], () => setAttr(el, 'aria-expanded', String(open.value)));
  return el;
}

/** Where each presets list opens from: its button (`PresetsButton`), by the state they share. */
const anchors = new WeakMap<ObservableState<boolean>, HTMLElement>();

/**
 * The presets as a popover under what opened them (its `PresetsButton`, or `anchor`), so a long list never pushes the
 * form down; the one the value matches is marked.
 */
export function PresetsList<T>(props: {
  open: ObservableState<boolean>;
  presets: () => readonly PD.SelectOption<T>[];
  /** The preset the value matches, by index (undefined when none: Custom). */
  current: ReadonlyState<number | undefined>;
  onPick(value: T): void;
  /** What it opens from, when not its `PresetsButton` (a select-like row). */
  anchor?: () => HTMLElement | undefined;
}) {
  const { open } = props;
  return (
    <Show when={open}>
      {() => {
        // Read when opened: a line graph's presets come from the data at hand (a volume's statistics).
        const presets = props.presets();
        const anchor = props.anchor?.() ?? anchors.get(open);
        const list = (
          <OptionList
            aria-label={anchor?.getAttribute('aria-label') ?? 'Presets'}
            options={presets.map(([, label], i) => [String(i), label] as const)}
            value={derive(props.current, (i) => (i === undefined ? undefined : String(i)))}
            onPick={(i) => {
              open.set(false);
              const preset = presets[Number(i)];
              if (preset) props.onPick(preset[0]);
            }}
            onClose={() => open.set(false)}
            floating
          />
        );
        return anchor ? <Floating anchor={anchor}>{list}</Floating> : <div class={classes.presetsList}>{list}</div>;
      }}
    </Show>
  );
}

const same = (a: unknown, b: unknown) => JSON.stringify(a) === JSON.stringify(b);

/** The preset a value matches: all of it (a list, a curve), or, for a group, each field the preset sets. */
export function presetIndex<T>(
  presets: readonly PD.SelectOption<T>[],
  value: unknown,
  partial = false,
): number | undefined {
  const i = presets.findIndex(([preset]) =>
    partial && preset && typeof preset === 'object' && value && typeof value === 'object'
      ? Object.entries(preset).every(([k, v]) => same((value as Record<string, unknown>)[k], v))
      : same(preset, value),
  );
  return i < 0 ? undefined : i;
}
