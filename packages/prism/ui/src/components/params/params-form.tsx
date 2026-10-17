/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { Asset } from '@molstar/core/util/assets';
import { camelCaseToWords } from '@molstar/core/util/string';
import { For, Show } from '../../dom/control-flow.js';
import type { Child } from '../../dom/jsx-types.js';
import { uniqueId } from '../../dom/focus.js';
import { effect, frameState, watch } from '../../dom/scope.js';
import { setAttr, setStyle, setText } from '../../dom/writes.js';
import { type ParamSection, ParamsFormModel, paramSections } from '../../forms/params-form.js';
import { type ObservableState, observableState } from '../../reactive/observable-state.js';
import type { ReadonlyState } from '../../reactive/readonly-state.js';
import { Floating, OptionList, SelectButton } from '../select/select.js';
import { ColorPicker, colorToHex, colorToRgb } from '../color/color-picker.js';
import { Icon } from '../icon/icon.js';
import { classes } from './_generated/params-form.module.css.js';
import {
  ColorListRow,
  IntervalControl,
  LineGraphRow,
  Mat4Row,
  MultiSelectRow,
  UrlControl,
  Vec3Control,
} from './extra-controls.js';
import type { FieldBinding } from './field-binding.js';
import { helpOf, withHelp } from './help.js';
import { PresetsButton, PresetsList, presetIndex } from './presets.js';

import { derive } from '../../reactive/state.js';
/** A field's value as a control sees it: a stream, the current value, and how to change it. */
export type { FieldBinding };

/**
 * The form for a `ParamsFormModel`, laid out as the classic UI's (v5 `ParameterControls`), at every level: essentials
 * first, the rest under "Advanced"; a choice (PD.Mapped) with its chosen case's params behind a sliders button; groups
 * folded unless flat. Each row binds only to its own field.
 */
export function ParamsForm(props: {
  form: ParamsFormModel<any>;
  /**
   * Shows only the rows whose label or description contains this text (or that hold such rows); groups open while it
   * is set. A search over settings.
   */
  filter?: ReadonlyState<string>;
}) {
  const { form } = props;
  return (
    <div class={classes.form}>
      <ParamControls
        params={form.params}
        field={(name) => ({
          // States at runtime (`derive`), typed as streams by the model.
          value: form.field(name) as ReadonlyState<any>,
          get: () => form.values.value[name],
          set: (v) => form.set(name, v),
        })}
        visible={(name) => (form.params[name].hideIf ? (form.visible(name) as ReadonlyState<boolean>) : undefined)}
        filter={normalized(props.filter)}
        place={{ form, path: '' }}
      />
    </div>
  );
}

/** Where rows are: their form (which keeps the folds' state) and their path in it (the folds' keys). */
interface Place {
  form: ParamsFormModel<any>;
  path: string;
}

/** Whether any of the params has a row (not hidden, not a plain value). */
function hasRows(params: PD.Params): boolean {
  return Object.keys(params).some((k) => !params[k].isHidden && params[k].type !== 'value');
}

/**
 * The rows of a set of params (a form's, a group's, a choice's case's), as v5's `ParameterControls`: the essential ones
 * (`isEssential`), then the rest folded under "Advanced", only when there are both; within each, a category's params
 * (PD `category`) under a fold of their own. `isHidden` params have no row, a `hideIf` one hides by the values it is
 * given; a fold with nothing to show hides.
 */
function ParamControls(props: {
  params: PD.Params;
  field(name: string): FieldBinding;
  /** Whether a row shows now (its `hideIf`); undefined when it always does. */
  visible(name: string): ReadonlyState<boolean> | undefined;
  filter?: ReadonlyState<string>;
  place: Place;
}): Child {
  const { params, filter, place } = props;
  const { essentials, advanced } = paramSections(params);
  const force = searching(filter);
  const visibles = new Map<string, ReadonlyState<boolean> | undefined>();
  const visible = (name: string) => {
    if (!visibles.has(name)) visibles.set(name, props.visible(name));
    return visibles.get(name);
  };
  const shown = (name: string): ReadonlyState<boolean> | undefined => {
    const v = visible(name);
    if (!filter) return v;
    return v
      ? derive([v, filter], (a, q) => a && matches(name, params[name], q))
      : derive(filter, (q) => matches(name, params[name], q));
  };
  // Hidden while none of its rows shows (all hidden by their hideIf, none matching a search); undefined: always shown
  // (a row that always shows among them).
  const hiddenWhenEmpty = (names: string[]) => {
    const states = names.map(shown);
    if (states.some((x) => !x)) return undefined;
    return derive(states as ReadonlyState<boolean>[], (...vs) => !vs.some(Boolean));
  };
  const row = (name: string) => (
    <ParamRow
      name={name}
      param={params[name]}
      binding={props.field(name)}
      visible={visible(name)}
      filter={filter}
      place={{ form: place.form, path: `${place.path}/${name}` }}
    />
  );
  const sections = (list: ParamSection[], key: string): Child[] =>
    list.flatMap((s) =>
      s.category
        ? [
            <div hidden={hiddenWhenEmpty(s.names)}>
              <Collapsible
                title={s.category}
                state={place.form.fold(`${place.path}#${key}${s.category}`)}
                force={force}
              >
                {() => s.names.map(row)}
              </Collapsible>
            </div>,
          ]
        : s.names.map(row),
    );
  if (!essentials.length || !advanced.length) return <>{sections(essentials.length ? essentials : advanced, '')}</>;
  return (
    <>
      {sections(essentials, '')}
      <div hidden={hiddenWhenEmpty(advanced.flatMap((s) => s.names))}>
        <Collapsible title="Advanced" state={place.form.fold(`${place.path}#advanced`)} force={force}>
          {() => sections(advanced, 'advanced/')}
        </Collapsible>
      </div>
    </>
  );
}

/** The filter text, trimmed and lower-cased; undefined when there is no filter. */
function normalized(filter?: ReadonlyState<string>): ReadonlyState<string> | undefined {
  return filter && derive(filter, (q) => q.trim().toLowerCase());
}

function searching(filter?: ReadonlyState<string>): ReadonlyState<boolean> | undefined {
  return filter && derive(filter, (q) => q.length > 0);
}

function labelOf(name: string, param: PD.Any): string {
  return param.label || camelCaseToWords(name);
}

/** Whether any of the params (or anything inside them) matches a search text: a settings group with no match hides. */
export function paramsMatch(params: PD.Params, text: string): boolean {
  const q = text.trim().toLowerCase();
  return Object.keys(params).some((k) => !params[k].isHidden && matches(k, params[k], q));
}

/** Whether a param or anything inside it matches the (normalized) filter text. */
function matches(name: string, param: PD.Any, q: string): boolean {
  if (!q) return true;
  if (labelOf(name, param).toLowerCase().includes(q)) return true;
  if (param.description?.toLowerCase().includes(q)) return true;
  if (param.type === 'group') return Object.keys(param.params).some((k) => matches(k, param.params[k], q));
  if (param.type === 'mapped') {
    return param.select.options.some(([option]) => {
      const inner = param.map(option);
      return inner.type === 'group'
        ? Object.keys(inner.params).some((k) => matches(k, inner.params[k], q))
        : matches(option, inner, q);
    });
  }
  return false;
}

function ParamRow(props: {
  name: string;
  param: PD.Any;
  binding: FieldBinding;
  visible?: ReadonlyState<boolean>;
  filter?: ReadonlyState<string>;
  place: Place;
}) {
  const { name, param, binding, filter, place } = props;
  const label = labelOf(name, param);
  const visible = props.visible;
  const shown = filter
    ? visible
      ? derive([visible, filter], (v, q) => v && matches(name, param, q))
      : derive(filter, (q) => matches(name, param, q))
    : visible;
  const hidden = shown && derive(shown, (v) => !v);
  // A group or a choice whose own name matches shows all of its rows.
  const inner = filter && derive(filter, (q) => (label.toLowerCase().includes(q) ? '' : q));
  if (param.type === 'group') {
    return (
      <div hidden={hidden}>
        <GroupControl
          label={label}
          param={param}
          binding={binding}
          filter={inner}
          force={searching(filter)}
          place={place}
        />
      </div>
    );
  }
  // Its help (a description, a legend; for a theme choice, the chosen theme's): a ? by the label, opening under the row.
  // The label's text has an id: the row's control is labelled by it (accessibility.md: every field has a name).
  const labelId = uniqueId('msp-field');
  const { head, block } = withHelp(label, helpOf(param, binding.value), labelId);
  const wrap = (row: Child) => (
    <div hidden={hidden}>
      {row}
      {block}
    </div>
  );
  if (param.type === 'mapped') {
    return wrap(
      isOnOff(param) ? (
        <OnOffControl
          label={label}
          labelId={labelId}
          head={head}
          description={param.description}
          param={param}
          binding={binding}
          place={place}
        />
      ) : (
        <MappedControl
          label={label}
          head={head}
          param={param}
          binding={binding}
          filter={inner}
          force={searching(filter)}
          place={place}
        />
      ),
    );
  }
  if (param.type === 'select') {
    return wrap(
      <SelectRow label={label} head={head} description={param.description} options={param.options} binding={binding} />,
    );
  }
  if (param.type === 'color') {
    return wrap(
      <ColorRow label={label} labelId={labelId} head={head} description={param.description} binding={binding} />,
    );
  }
  // Fields that open something under their row.
  if (param.type === 'multi-select')
    return wrap(<MultiSelectRow head={head} label={label} param={param} binding={binding} />);
  if (param.type === 'mat4') return wrap(<Mat4Row head={head} label={label} param={param} binding={binding} />);
  if (param.type === 'line-graph')
    return wrap(<LineGraphRow head={head} label={label} param={param} binding={binding} />);
  if (param.type === 'color-list')
    return wrap(<ColorListRow head={head} label={label} param={param} binding={binding} />);
  if (param.type === 'object-list')
    return wrap(<ObjectListRow head={head} label={label} param={param} binding={binding} />);
  if (param.type === 'conditioned')
    return wrap(<ConditionedControl head={head} label={label} param={param} binding={binding} place={place} />);
  if (param.type === 'script') return wrap(<ScriptControl head={head} label={label} binding={binding} />);
  // Converted: the control of what it converts to, the value converted both ways.
  if (param.type === 'converted') {
    const converted: FieldBinding = {
      value: derive(binding.value, (v) => param.fromValue(v)),
      get: () => param.fromValue(binding.get()),
      set: (v) => binding.set(param.toValue(v)),
    };
    const as = { ...param.converted, label: param.converted.label ?? label } as PD.Any;
    return (
      <ParamRow name={name} param={as} binding={converted} visible={props.visible} filter={filter} place={place} />
    );
  }
  // A plain value (data carried along, not edited): no row, as in the classic UI.
  if (param.type === 'value') return null;
  return wrap(
    <div class={classes.row}>
      <span class={classes.label}>{head}</span>
      <div class={classes.control}>{control(param, binding, labelId, label)}</div>
    </div>,
  );
}

/** A color: the swatch and its channels (decimal) as a button; the picker opens under the row, full width. */
function ColorRow(props: {
  label: string;
  labelId: string;
  head?: Child;
  description?: string;
  binding: FieldBinding;
}) {
  const { binding, labelId } = props;
  const open = observableState(false);
  let swatch!: HTMLElement;
  const channels = document.createTextNode('');
  const content = (
    <>
      <div class={classes.row}>
        <span class={classes.label}>{props.head ?? <span id={labelId}>{props.label}</span>}</span>
        <div class={classes.control}>
          <button
            type="button"
            class={classes.color}
            aria-labelledby={`${labelId} ${labelId}-value`}
            aria-expanded={open}
            onClick={() => open.set(!open.value)}
          >
            <span class={classes.swatch} ref={(el: HTMLElement) => (swatch = el)} />
            <span class={classes.channels} id={`${labelId}-value`} data-tip-overflow="">
              {channels}
            </span>
          </button>
        </div>
      </div>
      <Show when={open} reveal>
        {() => <ColorPicker value={binding.value} onChange={(c) => binding.set(c)} />}
      </Show>
    </>
  );
  watch([binding.value], () => {
    const v = binding.get() as number;
    setStyle(swatch, 'background', colorToHex(v));
    setText(channels, colorToRgb(v));
  });
  return content;
}

/**
 * Files (PD.File, PD.FileList): a field like a dropdown's, "Select files…" until some are chosen, then their names; a
 * click picks, files dropped on it are taken too, ✕ clears.
 */
function FileControl(props: { param: PD.FileParam | PD.FileListParam; binding: FieldBinding; labelId?: string }) {
  const { param, binding, labelId } = props;
  const multiple = param.type === 'file-list';
  const dragging = observableState(false);
  const take = (files: readonly File[]) => {
    if (!files.length) return;
    binding.set(multiple ? files.map((f) => Asset.File(f)) : Asset.File(files[0]));
  };
  const input = (
    <input
      type="file"
      hidden
      multiple={multiple}
      accept={param.accept}
      onChange={(e) => {
        take([...(e.currentTarget.files ?? [])]);
        e.currentTarget.value = '';
      }}
    />
  ) as HTMLInputElement;
  const names = derive(binding.value, (v: Asset.File | Asset.File[] | null) =>
    (v ? (Array.isArray(v) ? v : [v]) : []).map((f) => f.name),
  );
  return (
    <div class={classes.file} data-dragging={derive(dragging.asReadonly, (d) => (d ? '' : undefined))}>
      <button
        type="button"
        class={classes.filePick}
        aria-labelledby={labelId ? `${labelId} ${labelId}-value` : undefined}
        onClick={() => input.click()}
        onDragOver={(e) => {
          e.preventDefault();
          dragging.set(true);
        }}
        onDragLeave={() => dragging.set(false)}
        onDrop={(e) => {
          e.preventDefault();
          e.stopPropagation();
          dragging.set(false);
          take([...(e.dataTransfer?.files ?? [])]);
        }}
      >
        <span
          class={classes.fileNames}
          id={labelId ? `${labelId}-value` : undefined}
          data-empty={derive(names, (n) => (n.length ? undefined : ''))}
        >
          {derive(names, (n) =>
            n.length === 0
              ? multiple
                ? 'Select files…'
                : 'Select a file…'
              : n.length === 1
                ? n[0]
                : `${n.length} files · ${n.join(', ')}`,
          )}
        </span>
      </button>
      <Show when={derive(names, (n) => n.length > 0)}>
        {() => (
          <button
            type="button"
            class={classes.fileClear}
            aria-label="Clear"
            title="Clear"
            onClick={() => binding.set(null)}
          >
            <Icon name="x" size={13} />
          </button>
        )}
      </Show>
      {input}
    </div>
  );
}

/** A plain row's control, labelled by the row's label (`labelId`; `label` its text, for the parts of a compound one). */
function control(param: PD.Any, binding: FieldBinding, labelId: string, label: string): Child {
  switch (param.type) {
    case 'boolean':
      return <Toggle binding={binding} labelId={labelId} />;
    case 'number':
      return <NumberControl param={param} binding={binding} labelId={labelId} />;
    case 'text':
      // Several lines (a description, a script): a field that grows by hand, Enter makes a new line.
      return (
        <TextControl binding={binding} labelId={labelId} placeholder={param.placeholder} multiline={param.multiline} />
      );
    case 'file':
    case 'file-list':
      return <FileControl param={param} binding={binding} labelId={labelId} />;
    case 'vec3':
      return <Vec3Control param={param} binding={binding} label={label} />;
    case 'interval':
      return <IntervalControl param={param} binding={binding} label={label} />;
    case 'url':
      return <UrlControl binding={binding} labelId={labelId} />;
    default:
      return <span class={classes.unsupported}>{param.type}</span>;
  }
}

/**
 * A field's control as view = f(state) (rendering.md §2): its DOM made once, one `watch` on the field's stream whose
 * render reads the current value (`binding.get()`) and writes what differs.
 */
function TextControl(props: { binding: FieldBinding; labelId?: string; placeholder?: string; multiline?: boolean }) {
  const { binding } = props;
  // Several lines (a description, a script): a field that grows by hand, Enter makes a new line.
  const el = (
    props.multiline ? (
      <textarea
        class={classes.textarea}
        rows={3}
        aria-labelledby={props.labelId}
        placeholder={props.placeholder}
        onChange={(e) => binding.set(e.currentTarget.value)}
      />
    ) : (
      <input
        class={classes.input}
        aria-labelledby={props.labelId}
        placeholder={props.placeholder}
        onChange={(e) => binding.set(e.currentTarget.value)}
      />
    )
  ) as HTMLInputElement | HTMLTextAreaElement;
  watch([binding.value], () => setValue(el, binding.get()));
  return el;
}

/** An input's value ('' for none), written only when it differs (a write moves the caret). */
function setValue(el: HTMLInputElement | HTMLTextAreaElement, value: unknown) {
  const v = value === null || value === undefined ? '' : String(value);
  if (el.value !== v) el.value = v;
}

function Toggle(props: { binding: FieldBinding; labelId?: string }) {
  const { binding } = props;
  const el = (
    <button
      type="button"
      role="switch"
      class={classes.switch}
      aria-labelledby={props.labelId}
      onClick={() => binding.set(!binding.get())}
    >
      <span class={classes.knob} />
    </button>
  ) as HTMLElement;
  watch([binding.value], () => setAttr(el, 'aria-checked', String(!!binding.get())));
  return el;
}

function NumberControl(props: { param: PD.Numeric; binding: FieldBinding; labelId?: string }) {
  const { param, binding, labelId } = props;
  if (param.min === undefined || param.max === undefined) {
    // Text, not a native number field: that one shows decimals in the browser's locale ("0,5") while every other
    // number reads "0.5"; a comma typed is read as the decimal point.
    const input = (
      <input
        class={classes.input}
        type="text"
        inputMode="decimal"
        aria-labelledby={labelId}
        onChange={(e) => {
          const v = parseNumber(e.currentTarget.value);
          if (Number.isFinite(v)) binding.set(v);
          e.currentTarget.value = format(Number.isFinite(v) ? v : (binding.get() as number), param.step);
        }}
      />
    ) as HTMLInputElement;
    watch([binding.value], () => setValue(input, format(binding.get() as number, param.step)));
    return input;
  }
  const min = param.min;
  const max = param.max;
  // A drag in the frame lane (rendering.md §2.3): input events come several a frame, the field is set once a frame (it
  // may update a structure's representation).
  const dragged = frameState(Number.NaN, (a, b) => Object.is(a, b));
  effect(dragged, (v) => {
    if (!Number.isNaN(v)) binding.set(v);
  });
  let range!: HTMLInputElement;
  let number!: HTMLInputElement;
  const el = (
    <div class={classes.slider}>
      <input
        class={classes.range}
        type="range"
        aria-labelledby={labelId}
        min={String(min)}
        max={String(max)}
        step={String(param.step ?? (max - min) / 100)}
        ref={(e: HTMLInputElement) => (range = e)}
        onInput={(e) => dragged.offer(Number(e.currentTarget.value))}
      />
      <input
        class={classes.number}
        aria-labelledby={labelId}
        ref={(e: HTMLInputElement) => (number = e)}
        onChange={(e) => {
          const v = parseNumber(e.currentTarget.value);
          const shown = Number.isFinite(v) ? Math.min(max, Math.max(min, v)) : (binding.get() as number);
          if (Number.isFinite(v)) binding.set(shown);
          // The box says what holds, whether or not the value changed (a clamp to the same value, text that is no number).
          e.currentTarget.value = format(shown, param.step);
        }}
      />
    </div>
  );
  watch([binding.value], () => {
    const v = binding.get() as number;
    setValue(range, v);
    // The filled part of the track, for the quiet styling (no native accent fill).
    setStyle(range, '--msp-progress', `${(100 * (v - min)) / (max - min || 1)}%`);
    setValue(number, format(v, param.step));
  });
  return el;
}

function format(v: number, step?: number) {
  if (!Number.isFinite(v)) return '';
  // No step: the value as it is (at most 4 decimals), not rounded to a whole number.
  if (step === undefined) return String(Number(v.toFixed(4)));
  const digits = step < 1 ? Math.min(4, Math.ceil(-Math.log10(step))) : 0;
  return v.toFixed(digits);
}

/** A typed number: "0.5", "0,5" (a comma as the decimal point), " -1 "; NaN when it is not one. */
export function parseNumber(text: string): number {
  const t = text.trim().replace(',', '.');
  return t === '' ? Number.NaN : Number(t);
}

/**
 * A choice as a row: the label and the current option; opened, the options in a popover under the field (grouped by
 * Mol*'s categories). Values may be any type: the list works with indices.
 */
function SelectRow(props: {
  label: string;
  /** The label's content (with its help); the label text without. */
  head?: Child;
  description?: string;
  options: readonly PD.SelectOption<any>[];
  binding: FieldBinding;
  /** After the select, in the row (a choice's params toggle). */
  after?: Child;
}) {
  const { options, binding } = props;
  const open = observableState(false);
  const index = derive(binding.value, (v) => String(options.findIndex((o) => o[0] === v)));
  // Mol*'s option categories are the list's groups.
  const indexed = options.map(([, label, group], i) => [String(i), label, group] as const);
  let control: HTMLElement | undefined;
  return (
    <>
      <div class={classes.row}>
        <span class={classes.label}>{props.head ?? props.label}</span>
        <div class={classes.control} ref={(el: HTMLElement) => (control = el)}>
          <SelectButton
            size="sm"
            open={open}
            aria-label={props.label}
            label={derive(index, (i) => indexed[Number(i)]?.[1] ?? '')}
          />
          {props.after}
        </div>
      </div>
      {/* A popover under the field (the form does not move for it), at least as wide as the field. */}
      <Show when={open}>
        {() => {
          const list = (
            <OptionList
              options={indexed}
              value={index}
              onPick={(i) => {
                open.set(false);
                binding.set(options[Number(i)][0]);
              }}
              onClose={() => open.set(false)}
              floating
              aria-label={props.label}
            />
          );
          return control ? <Floating anchor={control}>{list}</Floating> : list;
        }}
      </Show>
    </>
  );
}

function GroupControl(props: {
  label: string;
  param: PD.Group<any>;
  binding: FieldBinding;
  /** For its rows (cleared when the group's own name matches: all of them show). */
  filter?: ReadonlyState<string>;
  /** Open while searching. */
  force?: ReadonlyState<boolean>;
  place: Place;
}) {
  const { param, binding, filter, place } = props;
  // As v5's: a group with nothing to show has no row.
  if (!hasRows(param.params)) return null;
  const presets = param.presets;
  const fields = (params: PD.Params) => <GroupFields params={params} binding={binding} filter={filter} place={place} />;
  // Its fold, open as the param says at first, then as the user leaves it (for the form's life).
  const folded = place.form.fold(`${place.path}#open`, !!param.isExpanded);
  // A preset sets the fields it names, over the group's current values (Mol*'s presets name them all).
  const apply = (value: unknown) => binding.set({ ...(binding.get() ?? {}), ...(value as object) });
  const current = presets ? derive(binding.value, (v) => presetIndex(presets, v, true)) : undefined;
  // A flat group has no header of its own: its params are rows like the others, its presets a row before them.
  if (param.isFlat)
    return (
      <>
        {presets && current ? (
          <PresetsRow label={props.label} presets={presets} current={current} onPick={apply} />
        ) : null}
        {fields(param.params)}
      </>
    );
  // A pivot (Mol*'s: a provider's id) is the group's main field: a row of its own, the rest folded under the group.
  const pivot = typeof param.pivot === 'string' && param.params[param.pivot] ? param.pivot : undefined;
  if (pivot && (!presets || !current)) {
    const { [pivot]: main, ...rest } = param.params;
    return (
      <>
        {fields({ [pivot]: main })}
        {hasRows(rest) ? (
          <Collapsible title={props.label} state={folded} force={props.force}>
            {() => fields(rest)}
          </Collapsible>
        ) : null}
      </>
    );
  }
  if (!presets || !current)
    return (
      <Collapsible title={props.label} state={folded} force={props.force}>
        {() => fields(param.params)}
      </Collapsible>
    );
  // With presets: Presets in the header (as Mol*'s), the matching one named after the title; a pick opens the group
  // so what it set shows.
  const presetsOpen = observableState(false);
  const opened = observableState(false);
  return (
    <Collapsible
      title={props.label}
      state={folded}
      force={props.force}
      openWhen={opened.asReadonly}
      meta={derive(current, (i) => (i === undefined ? '' : presets[i][1]))}
      action={<PresetsButton open={presetsOpen} label={props.label} />}
      below={() => (
        <PresetsList
          open={presetsOpen}
          presets={() => presets}
          current={current}
          onPick={(v) => {
            apply(v);
            opened.set(true);
          }}
        />
      )}
    >
      {() => fields(param.params)}
    </Collapsible>
  );
}

/** A flat group's presets: a row like a select's, reading the matching preset (or Custom). */
function PresetsRow(props: {
  label: string;
  presets: readonly PD.SelectOption<unknown>[];
  current: ReadonlyState<number | undefined>;
  onPick(value: unknown): void;
}) {
  const open = observableState(false);
  let control: HTMLElement | undefined;
  return (
    <>
      <div class={classes.row}>
        <span class={classes.label}>Presets</span>
        <div class={classes.control} ref={(el: HTMLElement) => (control = el)}>
          <SelectButton
            size="sm"
            open={open}
            aria-label={`${props.label} presets`}
            label={derive(props.current, (i) => (i === undefined ? 'Custom' : props.presets[i][1]))}
          />
        </div>
      </div>
      <PresetsList
        open={open}
        presets={() => props.presets}
        current={props.current}
        onPick={props.onPick}
        anchor={() => control}
      />
    </>
  );
}

/** Mol*'s on/off pattern (occlusion, shadows, outline): a mapped param whose choices are exactly 'on' and 'off'. */
function isOnOff(param: PD.Mapped<any>): boolean {
  const names = param.select.options.map((o) => o[0]).sort();
  return names.length === 2 && names[0] === 'off' && names[1] === 'on';
}

/**
 * An on/off mapped param: a switch, and while on, its own settings behind a sliders button, expanding inline under the
 * row (occlusion's samples, radius, bias, blur…).
 */
function OnOffControl(props: {
  label: string;
  labelId: string;
  head?: Child;
  description?: string;
  param: PD.Mapped<any>;
  binding: FieldBinding;
  place: Place;
}) {
  const { param, binding, place } = props;
  const on = derive(binding.value, (v) => v?.name === 'on');
  const open = place.form.fold(`${place.path}#params`);
  const set = (name: 'on' | 'off') => binding.set({ name, params: defaults(param.map(name)) });
  const toggle: FieldBinding = {
    value: on,
    get: () => binding.get()?.name === 'on',
    set: (v) => {
      set(v ? 'on' : 'off');
      if (!v) open.set(false);
    },
  };
  const inner = param.map('on');
  const hasSettings = inner.type === 'group' && hasRows(inner.params);
  const paramsBinding: FieldBinding = {
    value: derive(binding.value, (v) => v?.params),
    get: () => binding.get()?.params,
    set: (p) => binding.set({ name: 'on', params: p }),
  };
  return (
    <>
      <div class={classes.row}>
        <span class={classes.label}>{props.head ?? <span id={props.labelId}>{props.label}</span>}</span>
        <div class={classes.control}>
          <Toggle binding={toggle} labelId={props.labelId} />
          {hasSettings ? (
            <Show when={on}>
              {() => (
                <button
                  type="button"
                  class={classes.settings}
                  title={`${props.label} settings`}
                  aria-label={`${props.label} settings`}
                  aria-expanded={open}
                  onClick={() => open.set(!open.value)}
                >
                  <Icon name="sliders-horizontal" size={14} />
                </button>
              )}
            </Show>
          ) : null}
        </div>
      </div>
      {hasSettings ? (
        <Show when={derive([on, open.asReadonly], (a, b) => a && b)}>
          {() => (
            <div class={classes.subSettings}>
              <GroupFields
                params={(inner as PD.Group<any>).params}
                binding={paramsBinding}
                place={{ form: place.form, path: `${place.path}:on` }}
              />
            </div>
          )}
        </Show>
      ) : null}
    </>
  );
}

/**
 * A mapped param (PD.Mapped, MappedStatic), as v5's `MappedControl`: the choice; when the chosen case's params are a
 * group, a sliders button after it opens them (with the form's rules inside: essentials, then Advanced); a flat group's
 * params are rows at the same level, any other case's value one row under the choice. Its open state is the form's
 * (kept from one case to the next, as v5's).
 */
function MappedControl(props: {
  label: string;
  head?: Child;
  param: PD.Mapped<any>;
  binding: FieldBinding;
  /** For the chosen case's rows. */
  filter?: ReadonlyState<string>;
  /** Open while searching. */
  force?: ReadonlyState<boolean>;
  place: Place;
}) {
  const { param, binding, place } = props;
  const name = derive(binding.value, (v) => v?.name);
  const nameBinding: FieldBinding = {
    value: name,
    get: () => binding.get()?.name,
    set: (n) => {
      const params = param.map(n);
      binding.set({ name: n, params: defaults(params) });
    },
  };
  /** Does the case have params behind the button (a group that is not flat, with a row to show)? */
  const toggled = (n: string | undefined) => {
    if (n === undefined) return false;
    const inner = param.map(n);
    return inner.type === 'group' && !inner.isFlat && hasRows(inner.params);
  };
  const hasToggle = derive(name, toggled);
  const open = place.form.fold(`${place.path}#params`);
  const shown = props.force ? derive([open.asReadonly, props.force], (o, f) => o || f) : open.asReadonly;
  const toggle = (
    <Show when={hasToggle}>
      {() => (
        <button
          type="button"
          class={`${classes.settings} ${classes.afterSelect}`}
          title={`${props.label} settings`}
          aria-label={`${props.label} settings`}
          aria-expanded={shown}
          onClick={() => open.set(!open.value)}
        >
          <Icon name="sliders-horizontal" size={14} />
        </button>
      )}
    </Show>
  );
  return (
    <div class={classes.mapped}>
      <SelectRow
        label={props.label}
        head={props.head}
        options={param.select.options}
        binding={nameBinding}
        after={toggle}
      />
      <Show when={name}>
        {(n) => {
          const inner = param.map(n);
          const paramsBinding: FieldBinding = {
            value: derive(binding.value, (v) => v?.params),
            get: () => binding.get()?.params,
            set: (p) => binding.set({ name: n, params: p }),
          };
          const at: Place = { form: place.form, path: `${place.path}:${n}` };
          if (inner.type === 'group') {
            // A flat group (PD's isFlat) reads as more rows at the same level.
            if (inner.isFlat)
              return <GroupFields params={inner.params} binding={paramsBinding} filter={props.filter} place={at} />;
            if (!hasRows(inner.params)) return null;
            return (
              <Show when={shown} reveal>
                {() => (
                  <GroupFields params={inner.params} binding={paramsBinding} filter={props.filter} place={at} nested />
                )}
              </Show>
            );
          }
          // Another kind of value: its own row, named by the case (as v5's).
          return <ParamRow name={n} param={inner} binding={paramsBinding} filter={props.filter} place={at} />;
        }}
      </Show>
    </div>
  );
}

/**
 * The params of a group as rows (the form's rules: essentials, Advanced, categories, hideIf on the group's values),
 * bound to the group's value object.
 */
function GroupFields(props: {
  params: PD.Params;
  binding: FieldBinding;
  nested?: boolean;
  filter?: ReadonlyState<string>;
  place: Place;
}) {
  const { params, binding } = props;
  return (
    <div class={props.nested ? classes.nested : classes.fields}>
      <ParamControls
        params={params}
        field={(k) => ({
          value: derive(binding.value, (v) => v?.[k]),
          get: () => binding.get()?.[k],
          set: (cv) => binding.set({ ...binding.get(), [k]: cv }),
        })}
        visible={(k) => {
          const hideIf = params[k].hideIf;
          if (!hideIf) return undefined;
          return derive(binding.value, (v) => !hideIf(v ?? {}));
        }}
        filter={props.filter}
        place={props.place}
      />
    </div>
  );
}

/** Conditioned: the condition as a choice, then the value's own control for it (nested under the choice). */
function ConditionedControl(props: {
  head: Child;
  label: string;
  param: PD.Conditioned<any, any, any>;
  binding: FieldBinding;
  place: Place;
}) {
  const { param, binding, place } = props;
  const condition = derive(binding.value, (v) => String(param.conditionForValue(v)));
  const conditionBinding: FieldBinding = {
    value: condition,
    get: () => String(param.conditionForValue(binding.get())),
    set: (c) => binding.set(param.conditionedValue(binding.get(), c)),
  };
  return (
    <div class={classes.mapped}>
      <SelectRow label={props.label} head={props.head} options={param.select.options} binding={conditionBinding} />
      <Show when={condition}>
        {(c) => (
          <div class={classes.nested}>
            <ParamRow
              name={c}
              param={param.conditionParams[c]}
              binding={binding}
              place={{ form: place.form, path: `${place.path}:${c}` }}
            />
          </div>
        )}
      </Show>
    </div>
  );
}

const ScriptLanguages: PD.SelectOption<string>[] = [
  ['mol-script', 'MolScript'],
  ['pymol', 'PyMOL'],
  ['vmd', 'VMD'],
  ['jmol', 'Jmol'],
];

/** A script: its language as a choice, its expression under it (several lines). */
function ScriptControl(props: { head: Child; label: string; binding: FieldBinding }) {
  const { binding } = props;
  const part = (k: 'language' | 'expression'): FieldBinding => ({
    value: derive(binding.value, (v) => v?.[k] ?? ''),
    get: () => binding.get()?.[k],
    set: (x) => binding.set({ ...binding.get(), [k]: x }),
  });
  const expression = part('expression');
  const expressionId = uniqueId('msp-field');
  return (
    <div class={classes.mapped}>
      <SelectRow label={props.label} head={props.head} options={ScriptLanguages} binding={part('language')} />
      <div class={classes.nested}>
        <div class={classes.row}>
          <span class={classes.label} id={expressionId}>
            Expression
          </span>
          <div class={classes.control}>
            <textarea
              class={classes.textarea}
              aria-labelledby={expressionId}
              rows={3}
              spellcheck={false}
              value={expression.value}
              onChange={(e) => expression.set(e.currentTarget.value)}
            />
          </div>
        </div>
      </div>
    </div>
  );
}

/**
 * An object list (a list of like things: a volume's channels, a set of points): the count in the row; under it one quiet
 * list of the items by their labels (edit and remove on hover; an item's form unfolds inside the list), Add at its end.
 */
function ObjectListRow(props: { head: Child; label: string; param: PD.ObjectList<any>; binding: FieldBinding }) {
  const { param, binding } = props;
  const presets = param.presets;
  const presetsOpen = observableState(false);
  const items = derive(binding.value, (v: unknown[] | undefined) => v ?? []);
  const replace = (i: number, item: unknown) => {
    const list = [...(binding.get() ?? [])];
    list[i] = item;
    binding.set(list);
  };
  return (
    <>
      <div class={classes.row}>
        <span class={classes.label}>{props.head}</span>
        <div class={classes.control}>
          <span class={classes.count}>{derive(items, (l) => (l.length === 1 ? '1 item' : `${l.length} items`))}</span>
        </div>
      </div>
      {/* The items in one quiet list under the row; Add ends it. */}
      <div class={classes.items}>
        <For each={derive(items, (l) => l.map((_, i) => i))}>
          {(i) => <ObjectListItem param={param} binding={binding} index={i} replace={replace} />}
        </For>
        <div class={classes.itemActions}>
          <button
            type="button"
            class={classes.itemAdd}
            aria-label={`Add to ${props.label}`}
            onClick={() => binding.set([...(binding.get() ?? []), param.ctor()])}
          >
            <Icon name="plus" size={13} />
            Add
          </button>
          {presets ? <PresetsButton open={presetsOpen} label={props.label} /> : null}
        </div>
        {presets ? (
          <PresetsList
            open={presetsOpen}
            presets={() => presets}
            current={derive(binding.value, (v) => presetIndex(presets, v))}
            onPick={(v) => binding.set(v)}
          />
        ) : null}
      </div>
    </>
  );
}

function ObjectListItem(props: {
  param: PD.ObjectList<any>;
  binding: FieldBinding;
  index: number;
  replace: (i: number, item: unknown) => void;
}) {
  const { param, binding, index } = props;
  const open = observableState(false);
  return (
    <div class={classes.item} data-open={derive(open.asReadonly, (o) => (o ? '' : undefined))}>
      {/* A click anywhere on the row opens or closes the item's params (the sliders say so; remove does not). */}
      <div class={classes.itemRow} onClick={() => open.set(!open.value)}>
        <span class={classes.itemLabel}>
          {derive(binding.value, (v: unknown[]) => (v?.[index] !== undefined ? param.getLabel(v[index]) : ''))}
        </span>
        <button
          type="button"
          class={classes.settings}
          title="Edit"
          aria-label="Edit"
          aria-expanded={open}
          onClick={(e) => {
            e.stopPropagation();
            open.set(!open.value);
          }}
        >
          <Icon name="sliders-horizontal" size={14} />
        </button>
        <button
          type="button"
          class={classes.settings}
          title="Remove"
          aria-label="Remove"
          onClick={(e) => {
            e.stopPropagation();
            binding.set((binding.get() ?? []).filter((_: unknown, j: number) => j !== index));
          }}
        >
          <Icon name="x" size={13} />
        </button>
      </div>
      <Show when={open} reveal>
        {() => {
          // The item's own form: each change replaces the item in the list (a new list: the state stores it as is).
          const form = new ParamsFormModel(param.element, binding.get()?.[index], {
            onChange: (v) => props.replace(index, v),
          });
          return (
            <div class={classes.itemForm}>
              <ParamsForm form={form} />
            </div>
          );
        }}
      </Show>
    </div>
  );
}

function defaults(param: PD.Any): any {
  if (param.type === 'group') {
    const out: Record<string, unknown> = {};
    for (const k of Object.keys(param.params)) out[k] = param.params[k].defaultValue;
    return out;
  }
  return param.defaultValue;
}

/** A titled block that opens and closes; its content is rendered while open. */
function Collapsible(props: {
  title: string;
  open?: boolean;
  /** Its open state, when kept outside (the form's folds); else its own, from `open`. */
  state?: ObservableState<boolean>;
  force?: ReadonlyState<boolean>;
  /** Quiet words after the title (the preset the values match). */
  meta?: ReadonlyState<string>;
  /** A control at the header's end (Presets). */
  action?: Child;
  /** Under the header, open or folded (the presets' list). */
  below?: () => Child;
  /** Opens it when this emits true (a preset picked: what it set shows). */
  openWhen?: ReadonlyState<boolean>;
  children: () => Child;
}) {
  const own = props.state ?? observableState(!!props.open);
  const openWhen = props.openWhen;
  if (openWhen) {
    // When it turns true (not again at a resume: a group closed since stays closed).
    let was = openWhen.value;
    effect(openWhen, (o) => {
      if (o && !was) own.set(true);
      was = o;
    });
  }
  // Forced open while a filter is set (its matches are inside); the user's choice again once it clears.
  const open = props.force ? derive([own.asReadonly, props.force], (o, f) => o || f) : own.asReadonly;
  return (
    <div class={classes.collapsible}>
      <div class={classes.collapsibleBar}>
        <button
          type="button"
          class={classes.collapsibleHeader}
          aria-expanded={open}
          onClick={() => own.set(!own.value)}
        >
          {/* The name first, aligned with the labels; the chevron after it. */}
          {props.title}
          <Icon name="chevron-right" size={12} class={classes.chevron} />
          {props.meta ? <span class={classes.collapsibleMeta}>{props.meta}</span> : null}
        </button>
        {props.action}
      </div>
      {props.below?.()}
      <Show when={open} reveal>
        {() => <div class={classes.nested}>{props.children()}</div>}
      </Show>
    </div>
  );
}
