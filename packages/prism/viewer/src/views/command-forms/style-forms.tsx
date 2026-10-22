/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { colorToHex, type CommandFormView, type CommandFormViewProps, Select } from '@molstar/prism-ui';
import type { InvocationState } from '@molstar/prism-ui/commands/invocation';
import { Icon } from '@molstar/prism-ui/components/icon/icon';
import { For, Show } from '@molstar/prism-ui/dom/control-flow';
import type { Child } from '@molstar/prism-ui/dom/jsx-types';
import { effect } from '@molstar/prism-ui/dom/scope';
import type { ReadonlyState } from '@molstar/prism-ui/reactive/readonly-state';
import { derive } from '@molstar/prism-ui/reactive/state';
import { decodeColor } from '@molstar/core/util/color/utils';
import type { StyleChoice, StyleFormsModel, TypeOption } from '../../model/command-forms.js';
import { QuickColorRow, QuickRange } from '../common/quick-controls.js';
import { classes } from './_generated/style-forms.module.css.js';

// The command menu's own forms for a component's look (commands.md §5.1, api.md §9 "Own forms"): Show as… as the types
// the component supports, drawn as the inspector's representation rows (icon and name), with the picked type's main
// options under them; Color by theme… as the themes that apply, each painted with its colors on the structure, and the
// uniform color's swatches, by the registry's categories; a click on a theme (or a swatch) applies it, as the quick
// settings' tiles do, a click on a type picks it (its options follow). The menu's input keeps the focus: typing filters
// the choices, the arrows move the pick, Enter runs, Tab goes to the sentence's other fields and, past the last, to the
// options and swatches under the choices (and back), ⌫ and Esc as in any form.

/** The own forms, by command id: given to the command menu (`CommandMenuProps.forms`). */
export function styleCommandForms(model: StyleFormsModel): Record<string, CommandFormView> {
  return {
    'representation.set': { fields: ['type'], render: (p) => <ShowAsForm {...p} model={model} /> },
    'representation.color-theme': { fields: ['theme'], render: (p) => <ColorThemeForm {...p} model={model} /> },
  };
}

interface Group {
  readonly category: string;
  readonly items: readonly StyleChoice[];
}

const NoChoices: readonly StyleChoice[] = [];

/** Choices under their categories, the categories in the order they first come; one group without them. */
function groupsOf(list: readonly StyleChoice[]): Group[] {
  const by = new Map<string, StyleChoice[]>();
  for (const c of list) {
    const k = c.category ?? '';
    const g = by.get(k);
    if (g) g.push(c);
    else by.set(k, [c]);
  }
  return [...by].map(([category, items]) => ({ category, items }));
}

const sameChoices = (a: readonly StyleChoice[], b: readonly StyleChoice[]) =>
  a === b || (a.length === b.length && a.every((x, i) => x === b[i]));

type FormProps = CommandFormViewProps & { model: StyleFormsModel };

/** What the choices are read for: the component and the structure the sentence names. */
const subjectOf = (s: InvocationState) => `${s.values.component ?? ''}\u0000${s.values.structure ?? ''}`;

export function ShowAsForm(props: FormProps) {
  const { invocation, model } = props;
  const state = invocation.state;
  const types = derive(derive(state, subjectOf), () => model.types(invocation.values()));
  // Add: a type it does not show yet, unless one was asked for (typed, by a preset, or in the values it opened with).
  // The invocation does not say where a value came from, so the type's default (the registry's first) counts as not
  // asked: an asked type other than it is kept, shown or not; asked for, the default is still replaced while shown.
  const first = types.value.find((t) => !t.shown);
  const asked = state.value.values.type;
  const fallback = (invocation.params.type as { defaultValue?: unknown } | undefined)?.defaultValue;
  const notAsked = asked === undefined || asked === fallback;
  if (first && notAsked && types.value.find((t) => t.id === asked)?.shown !== false) {
    invocation.set('type', first.id, { stay: true });
  }
  const type = derive(state, (s) => (typeof s.values.type === 'string' ? s.values.type : ''));
  // The options are the type's: another type starts from its defaults.
  let was = type.value;
  effect(type, (t) => {
    if (t !== was) invocation.setParam('params', undefined);
    was = t;
  });
  const options = derive(
    derive(state, (s) => `${s.values.type ?? ''}\u0000${subjectOf(s)}`),
    () => model.typeOptions(type.value, invocation.values()),
  );
  const patch = derive(state, (s) => (s.params.params ?? {}) as Readonly<Record<string, unknown>>);
  const setOption = (name: string, value: unknown) => invocation.setParam('params', { ...patch.value, [name]: value });
  return (
    <Picker
      {...props}
      field="type"
      label="Representation"
      choices={types}
      columns={3}
      tile={(c) => (
        <>
          <Icon name={c.icon ?? 'shapes'} size={15} class={classes.icon} />
          <span class={classes.label} data-tip-overflow="">
            {c.label}
          </span>
          {c.shown ? <span class={classes.meta}>shown</span> : null}
        </>
      )}
      after={() => (
        <Show when={derive(options, (o) => o.length > 0)}>
          {() => (
            <section class={classes.options} aria-label="Options">
              <For each={options} key={(o) => `${type.value}:${o.name}`}>
                {(o) => <OptionControl option={o} patch={patch} onChange={setOption} />}
              </For>
            </section>
          )}
        </Show>
      )}
    />
  );
}

function OptionControl(props: {
  option: TypeOption;
  patch: ReadonlyState<Readonly<Record<string, unknown>>>;
  onChange(name: string, value: unknown): void;
}) {
  const o = props.option;
  if (o.kind === 'range') {
    const digits = o.step >= 1 ? 0 : o.step >= 0.1 ? 1 : 2;
    return (
      <QuickRange
        label={o.label}
        min={o.min}
        max={o.max}
        step={o.step}
        value={derive(props.patch, (p) => (typeof p[o.name] === 'number' ? (p[o.name] as number) : o.default))}
        format={(v) => (o.name === 'alpha' ? `${Math.round(v * 100)}%` : v.toFixed(digits))}
        onInput={(v) => props.onChange(o.name, v)}
        onCommit={() => {}}
      />
    );
  }
  return (
    <div class={classes.choice}>
      <span class={classes.choiceLabel}>{o.label}</span>
      <Select
        size="sm"
        over
        aria-label={o.label}
        options={o.options}
        value={derive(props.patch, (p) => (typeof p[o.name] === 'string' ? (p[o.name] as string) : o.default))}
        onChange={(v) => props.onChange(o.name, v)}
      />
    </div>
  );
}

export function ColorThemeForm(props: FormProps) {
  const { invocation, model } = props;
  const state = invocation.state;
  const themes = derive(
    derive(state, (s) => `${s.values.structure ?? ''}`),
    () => model.themes(invocation.values()),
  );
  const uniform = derive(state, (s) => s.values.theme === 'uniform');
  // The color is the uniform theme's: another theme runs without it (the operation takes a color or a theme, not both,
  // api.md §2 `representation`), as Show as… starts another type from its defaults.
  effect(uniform, (u) => {
    if (!u && state.value.params.color !== undefined) invocation.setParam('color', undefined);
  });
  const color = derive(state, (s) => (typeof s.params.color === 'string' ? decodeColor(s.params.color) : undefined));
  const pickColor = (c: number) => {
    invocation.set('theme', 'uniform', { stay: true });
    invocation.setParam('color', colorToHex(c));
  };
  return (
    <Picker
      {...props}
      field="theme"
      label="Color theme"
      choices={themes}
      columns={3}
      // A theme has nothing more to set: a click applies it. One color opens its swatches (a swatch applies).
      applyOnClick={(c) => c.id !== 'uniform'}
      tile={(c) => (
        <>
          <span
            class={classes.swatch}
            classList={{ [classes.unpainted]: !c.swatch && c.id !== 'uniform' }}
            style={{
              background:
                c.id === 'uniform' ? derive(color, (v) => (v === undefined ? '' : colorToHex(v))) : (c.swatch ?? ''),
            }}
            aria-hidden="true"
          >
            {c.id === 'uniform' ? <span class={classes.one} /> : null}
          </span>
          <span class={classes.label} data-tip-overflow="">
            {c.label}
          </span>
        </>
      )}
      after={() => (
        <Show when={uniform} reveal>
          {() => (
            <section class={classes.options} aria-label="Color">
              <QuickColorRow
                label="Color"
                value={color}
                onPick={(c) => {
                  pickColor(c);
                  props.run();
                }}
                onInput={pickColor}
              />
            </section>
          )}
        </Show>
      )}
    />
  );
}

/**
 * The choices of one field as a grid of ghost buttons, the picked one marked (`aria-checked`) and, while words are
 * typed into the field, the one Enter takes highlighted. The keys go to the menu's input: ←→ (with nothing typed) and
 * ↑↓ move the pick; Enter takes the highlighted choice and runs.
 */
function Picker(
  props: FormProps & {
    field: string;
    label: string;
    choices: ReadonlyState<readonly StyleChoice[]>;
    columns: number;
    tile(c: StyleChoice): Child;
    after?(): Child;
    /** A click on this choice runs the command (else it only picks it; a double click runs). */
    applyOnClick?(c: StyleChoice): boolean;
  },
) {
  const { invocation, field, choices, columns } = props;
  const state = invocation.state;
  // Typing goes into this field.
  invocation.focus(field);
  const active = derive(state, (s) => s.field === field);
  const typed = derive(state, (s) => (s.field === field ? s.text.trim() : ''));
  // What the typed words match, in the field's order (the best first), among the choices that apply.
  const visible = derive(
    [choices, typed, derive(state, (s) => (s.field === field ? s.completions : undefined))],
    (all, t, c) => {
      if (!t || !c) return all;
      const byId = new Map(all.map((x) => [x.id, x]));
      return c.flatMap((x) => byId.get(x.value) ?? []);
    },
    sameChoices,
  );
  const current = derive(state, (s) => (typeof s.values[field] === 'string' ? (s.values[field] as string) : undefined));
  // The choice Enter takes: the highlighted completion while words are typed (the first that applies), else none.
  const highlighted = derive([visible, typed, state], (v, t, s) => {
    if (!t) return undefined;
    const picked = s.completions[s.pick]?.value;
    return (v.find((x) => x.id === picked) ?? v[0])?.id;
  });

  // The choices under their categories (in the order they first come), each group its own grid.
  const groups = derive(visible, groupsOf);
  // The rows as drawn (a group's choices in rows of `columns`): ↑↓ go to the row above or below, the same column.
  const rowsOf = (gs: readonly Group[]) =>
    gs.flatMap((g) =>
      Array.from({ length: Math.ceil(g.items.length / columns) }, (_, r) =>
        g.items.slice(r * columns, (r + 1) * columns),
      ),
    );
  const move = (by: 'row' | 'item', delta: 1 | -1) => {
    const list = visible.value;
    if (!list.length) return;
    const id = highlighted.value ?? current.value;
    let next: StyleChoice | undefined;
    if (by === 'item') {
      const at = list.findIndex((c) => c.id === id);
      next = list[at < 0 ? (delta > 0 ? 0 : list.length - 1) : Math.min(list.length - 1, Math.max(0, at + delta))];
    } else {
      const rows = rowsOf(groups.value);
      const r = rows.findIndex((row) => row.some((c) => c.id === id));
      if (r < 0) next = delta > 0 ? list[0] : list[list.length - 1];
      else {
        const row = rows[Math.min(rows.length - 1, Math.max(0, r + delta))];
        const col = rows[r].findIndex((c) => c.id === id);
        next = row[Math.min(col, row.length - 1)];
      }
    }
    if (!next) return;
    const to = next.id;
    if (typed.value) invocation.point(state.value.completions.findIndex((c) => c.value === to));
    else invocation.set(field, to, { stay: true });
  };
  props.keys((e) => {
    if (!active.value || e.altKey || e.metaKey || e.ctrlKey) return false;
    switch (e.key) {
      case 'ArrowDown':
      case 'ArrowUp':
        move('row', e.key === 'ArrowDown' ? 1 : -1);
        return true;
      case 'ArrowLeft':
      case 'ArrowRight':
        if (typed.value) return false;
        move('item', e.key === 'ArrowRight' ? 1 : -1);
        return true;
      case 'Enter': {
        if (!typed.value) return false;
        const id = highlighted.value;
        if (!id) return false;
        if (invocation.set(field, id, { stay: true })) props.run();
        return true;
      }
      default:
        return false;
    }
  });

  return (
    <div
      class={classes.form}
      onKeyDown={(e) => {
        // Enter on an option's control (a slider, a select's button) runs, as in the sentence.
        if (e.key !== 'Enter' || e.defaultPrevented || (e.target as HTMLElement).closest('button, [role="listbox"]'))
          return;
        e.preventDefault();
        props.run();
      }}
    >
      <div role="radiogroup" aria-label={props.label}>
        <For each={groups} key={(g) => g.category} fallback={() => <p class={classes.empty}>Nothing matches.</p>}>
          {(g) => (
            <div role="group" aria-label={g.category || props.label}>
              {/* Ungrouped choices (the types) are captioned by what they are. */}
              <div class={classes.caption}>{g.category || props.label}</div>
              <div class={classes.grid} style={{ '--columns': String(columns) }}>
                <For
                  each={derive(
                    groups,
                    (gs) => gs.find((x) => x.category === g.category)?.items ?? NoChoices,
                    sameChoices,
                  )}
                  key={(c) => c.id}
                >
                  {(c) => (
                    <button
                      type="button"
                      role="radio"
                      tabIndex={-1}
                      class={classes.tile}
                      aria-checked={derive(current, (v) => v === c.id)}
                      data-highlight={derive(highlighted, (h) => h === c.id)}
                      onClick={() => {
                        if (invocation.set(field, c.id, { stay: true }) && props.applyOnClick?.(c)) props.run();
                      }}
                      onDblClick={() => {
                        if (invocation.set(field, c.id, { stay: true })) props.run();
                      }}
                    >
                      {props.tile(c)}
                    </button>
                  )}
                </For>
              </div>
            </div>
          )}
        </For>
      </div>
      {props.after?.()}
    </div>
  );
}
