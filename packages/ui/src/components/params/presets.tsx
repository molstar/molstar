/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { derive, type ReadonlyState, type State } from '@molstar/core/reactive';
import type { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { type Child, Show, setAttr, watch } from '../../dom/index.js';
import { state } from '../../state.js';
import { Button } from '../button/button.js';
import { Icon } from '../icon/icon.js';
import { Floating, OptionList } from '../select/select.js';
import { classes } from './_generated/params-form.module.css.js';

/**
 * PD's presets (a group's `presets`, an object list's `presets`, a line graph's `getPresets`): named values to start
 * from, as Mol*'s Presets menu. A quiet button (its label is the UI's: never cut); opened, the presets float over it
 * (a select's options: nothing in a form moves), the one the value matches marked.
 */
export function PresetsButton<T>(props: {
  /** What the presets are for: the button is named "<label> presets". */
  label: string;
  /** Read when opened: a line graph's presets come from the data at hand (a volume's statistics). */
  presets: () => readonly PD.SelectOption<T>[];
  /** The preset the value matches, by index (undefined when none: the values are custom). */
  current: ReadonlyState<number | undefined>;
  onPick(value: T): void;
}): Child {
  const open = state(false);
  const button = (
    <Button
      variant="ghost"
      size="xs"
      class={classes.presets}
      icon={<Icon name="swatch-book" />}
      onClick={(e) => {
        // In a group's header: the click is the presets', not the header's.
        e.stopPropagation();
        open.set(!open.value);
      }}
    >
      Presets
    </Button>
  ) as HTMLButtonElement;
  button.setAttribute('aria-label', `${props.label} presets`);
  // The list's opener: a press on it toggles the list (it is not a press outside that closes it first).
  button.setAttribute('aria-haspopup', 'listbox');
  watch([open], () => setAttr(button, 'aria-expanded', String(open.value)));
  return (
    <>
      {button}
      <PresetsList
        open={open}
        anchor={button}
        label={`${props.label} presets`}
        presets={props.presets}
        current={props.current}
        onPick={props.onPick}
      />
    </>
  );
}

/** The presets as options over what opened them (`anchor`), the one the value matches marked; a pick applies it. */
export function PresetsList<T>(props: {
  open: State<boolean>;
  anchor: HTMLElement;
  label: string;
  presets: () => readonly PD.SelectOption<T>[];
  current: ReadonlyState<number | undefined>;
  onPick(value: T): void;
}): Child {
  const { open } = props;
  return (
    <Show when={open}>
      {() => {
        const presets = props.presets();
        return (
          <Floating anchor={props.anchor} placement="over">
            <OptionList
              aria-label={props.label}
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
          </Floating>
        );
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
