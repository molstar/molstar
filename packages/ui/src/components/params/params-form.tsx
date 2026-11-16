/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { derive, type ReadonlyState, type State } from '@molstar/core/reactive';
import { Asset } from '@molstar/core/util/assets';
import type { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { camelCaseToWords } from '@molstar/core/util/string';
import {
  type Child,
  effect,
  For,
  frameState,
  Show,
  setAttr,
  setStyle,
  setText,
  uniqueId,
  type WatchSource,
  watch,
} from '../../dom/index.js';
import { type ParamSection, ParamsFormModel, paramSections } from '../../forms/params-form.js';
import { state } from '../../state.js';
import { Button } from '../button/button.js';
import { ColorPicker, colorToHex, colorToRgb } from '../color/color-picker.js';
import { Icon } from '../icon/icon.js';
import { Floating, OptionList, SelectButton } from '../select/select.js';
import { classes } from './_generated/params-form.module.css.js';
import {
  ColorListRow,
  IntervalControl,
  LineGraphRow,
  Mat4Row,
  MultiSelectRow,
  PanelRow,
  Row,
  UrlControl,
  Vec3Control,
} from './extra-controls.js';
import { type FieldBinding, keyBinding } from './field-binding.js';
import { helpOf, withHelp } from './help.js';
import { PresetsButton, PresetsList, presetIndex } from './presets.js';

export type { FieldBinding };

/**
 * The form for a `ParamsFormModel`, laid out as the classic UI's (v5 `ParameterControls`), at every level: essentials
 * first, the rest under "Advanced"; a choice (PD.Mapped) with its chosen case's params behind a sliders button; groups
 * folded unless flat. Each row is a view of its own field (a derived state of it: a change of another field wakes it
 * not), its control's properties written by one `watch`.
 *
 * Alignment (styles §5): the form pads `--msui-inner`, so in a surface its labels start at the inset; every row is on
 * the form's two columns (the label column, a gap, the values to the inner padding), its control 20px in a 24px row,
 * a control's text at `--msui-control-inset`. The label column is one for the whole form: it fits the widest label at
 * any depth, at least `--msui-label` and at most half the form (every block a subgrid of the form's columns); a label
 * past the cap stays on one line and fades, its full text in its tooltip. A nested level (a group's rows, a choice's
 * params, Advanced, a category) indents its labels by `--msui-form-indent` inside that column, so values share one
 * column at every depth. Selects and lists open over their field and panels float under it: nothing opens inline but
 * the form's own rows (folds, a choice's params) and a field's help.
 */
export function ParamsForm(props: {
  form: ParamsFormModel<any>;
  /**
   * Shows only the rows whose label or description contains this text (or that hold such rows); groups open while it
   * is set. A search over settings.
   */
  filter?: ReadonlyState<string>;
}) {
  return (
    <div class={classes.form}>
      <FormRows form={props.form} filter={normalized(props.filter)} level={0} />
    </div>
  );
}

/** A form's rows at a level (the top of a form, an object list's item). */
function FormRows(props: { form: ParamsFormModel<any>; filter?: ReadonlyState<string>; level: number }) {
  const { form } = props;
  return (
    <ParamControls
      params={form.params}
      field={(name) => ({ value: form.field(name), set: (v) => form.set(name, v) })}
      visible={(name) => (form.params[name].hideIf ? form.visible(name) : undefined)}
      filter={props.filter}
      place={{ form, path: '', level: props.level }}
    />
  );
}

/** Where rows are: their form (which keeps the folds' state), their path in it (the folds' keys), their level. */
interface Place {
  form: ParamsFormModel<any>;
  path: string;
  /** The nesting depth: labels indent by this many `--msui-form-indent`. */
  level: number;
}

const at = (place: Place, path: string, deeper = false): Place => ({
  form: place.form,
  path,
  level: deeper ? place.level + 1 : place.level,
});

/** A block of rows one level deeper: its labels step in, the label column shrinks as much (styles §5). */
function Nested(props: { level: number; children: Child }) {
  return (
    <div class={classes.fields} style={{ '--msui-params-level': String(props.level) }}>
      {props.children}
    </div>
  );
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
  // A section is hidden while none of its rows shows (all hidden by their hideIf, none matching a search): one watch of
  // their hideIfs and the search. None when a row among them always shows and there is no search.
  const hiddenWhenEmpty = (el: HTMLElement, names: string[]) => {
    if (!filter && names.some((n) => !visible(n))) return el;
    const sources: WatchSource[] = [];
    for (const n of names) {
      const v = visible(n);
      if (v) sources.push(v);
    }
    if (filter) sources.push(filter);
    watch(sources, () => {
      const q = filter?.value ?? '';
      setAttr(el, 'hidden', !names.some((n) => visible(n)?.value !== false && matches(n, params[n], q)));
    });
    return el;
  };
  const row = (name: string, base: Place) => (
    <ParamRow
      name={name}
      param={params[name]}
      binding={props.field(name)}
      visible={visible(name)}
      filter={filter}
      place={at(base, `${place.path}/${name}`)}
    />
  );
  // A section's rows at `base`'s level; a category's under its fold, a level deeper.
  const sections = (list: ParamSection[], key: string, base: Place): Child[] =>
    list.flatMap((s) => {
      if (!s.category) return s.names.map((n) => row(n, base));
      const inside = at(base, base.path, true);
      return [
        hiddenWhenEmpty(
          (
            <div class={classes.fields}>
              <Collapsible
                title={s.category}
                state={place.form.fold(`${place.path}#${key}${s.category}`)}
                force={force}
                level={base.level}
              >
                {() => s.names.map((n) => row(n, inside))}
              </Collapsible>
            </div>
          ) as HTMLElement,
          s.names,
        ),
      ];
    });
  if (!essentials.length || !advanced.length)
    return <>{sections(essentials.length ? essentials : advanced, '', place)}</>;
  return (
    <>
      {sections(essentials, '', place)}
      {hiddenWhenEmpty(
        (
          <div class={classes.fields}>
            <Collapsible
              title="Advanced"
              state={place.form.fold(`${place.path}#advanced`)}
              force={force}
              level={place.level}
            >
              {() => sections(advanced, 'advanced/', at(place, place.path, true))}
            </Collapsible>
          </div>
        ) as HTMLElement,
        advanced.flatMap((s) => s.names),
      )}
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
}): Child {
  const { name, param, binding, filter, place, visible } = props;
  // A plain value (data carried along, not edited): no row, as in the classic UI.
  if (param.type === 'value') return null;
  // Converted: the control of what it converts to, the value converted both ways.
  if (param.type === 'converted') {
    const converted: FieldBinding = {
      value: derive(binding.value, (v) => param.fromValue(v)),
      set: (v) => binding.set(param.toValue(v)),
    };
    const as = { ...param.converted, label: param.converted.label ?? labelOf(name, param) } as PD.Any;
    return <ParamRow name={name} param={as} binding={converted} visible={visible} filter={filter} place={place} />;
  }
  const label = labelOf(name, param);
  // Hidden by its hideIf and the search: one watch of both (none when neither applies).
  const sources: WatchSource[] = [];
  if (visible) sources.push(visible);
  if (filter) sources.push(filter);
  const wrapper = (...children: Child[]) => {
    const el = (<div class={classes.field}>{children}</div>) as HTMLElement;
    if (sources.length)
      watch(sources, () =>
        setAttr(el, 'hidden', !((visible?.value ?? true) && matches(name, param, filter?.value ?? ''))),
      );
    return el;
  };
  // A group or a choice whose own name matches shows all of its rows.
  const inner = filter && derive(filter, (q) => (label.toLowerCase().includes(q) ? '' : q));
  if (param.type === 'group') {
    return wrapper(
      <GroupControl
        label={label}
        param={param}
        binding={binding}
        filter={inner}
        force={searching(filter)}
        place={place}
      />,
    );
  }
  // Its help (a description, a legend; for a theme choice, the chosen theme's): a ? after the label, opening under the
  // row. The label's text has an id: the row's control is labelled by it (every field has a name).
  const labelId = uniqueId('msui-field');
  const { head, block } = withHelp(label, helpOf(param, binding.value), labelId);
  const wrap = (row: Child) => wrapper(row, block);
  switch (param.type) {
    case 'mapped':
      return wrap(
        isOnOff(param) ? (
          <OnOffControl label={label} labelId={labelId} head={head} param={param} binding={binding} place={place} />
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
    case 'select':
      return wrap(<SelectRow label={label} head={head} options={param.options} binding={binding} />);
    case 'color':
      return wrap(<ColorRow label={label} labelId={labelId} head={head} binding={binding} />);
    case 'multi-select':
      return wrap(<MultiSelectRow head={head} label={label} param={param} binding={binding} />);
    case 'mat4':
      return wrap(<Mat4Row head={head} label={label} param={param} binding={binding} />);
    case 'line-graph':
      return wrap(<LineGraphRow head={head} label={label} param={param} binding={binding} />);
    case 'color-list':
      return wrap(<ColorListRow head={head} label={label} param={param} binding={binding} />);
    case 'object-list':
      return wrap(<ObjectListRow head={head} label={label} param={param} binding={binding} place={place} />);
    case 'conditioned':
      return wrap(<ConditionedControl head={head} label={label} param={param} binding={binding} place={place} />);
    case 'script':
      return wrap(<ScriptControl head={head} label={label} binding={binding} place={place} />);
    default:
      return wrap(<Row head={head}>{control(param, binding, labelId, label)}</Row>);
  }
}

/** A color: the swatch and its channels (decimal) in a field; the picker in a panel under it. */
function ColorRow(props: { label: string; labelId: string; head: Child; binding: FieldBinding }) {
  const { binding, labelId } = props;
  const swatch = (<span class={classes.swatch} />) as HTMLElement;
  const channels = document.createTextNode('');
  watch([binding.value], () => {
    const v = binding.value.value as number;
    setStyle(swatch, 'background', colorToHex(v));
    setText(channels, colorToRgb(v));
  });
  return (
    <PanelRow
      head={props.head}
      label={props.label}
      labelledBy={`${labelId} ${labelId}-value`}
      field={
        <>
          {swatch}
          <span class={classes.channels} id={`${labelId}-value`}>
            {channels}
          </span>
        </>
      }
      panel={() => <ColorPicker value={binding.value} onChange={(c) => binding.set(c)} />}
    />
  );
}

/**
 * Files (PD.File, PD.FileList): a field like a select's, "Select files…" until some are chosen, then their names; a
 * click picks, files dropped on it are taken too, ✕ clears.
 */
function FileControl(props: { param: PD.FileParam | PD.FileListParam; binding: FieldBinding; labelId?: string }) {
  const { param, binding, labelId } = props;
  const multiple = param.type === 'file-list';
  const dragging = state(false);
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
  const namesText = document.createTextNode('');
  const nameField = (
    <span class={classes.fileNames} id={labelId ? `${labelId}-value` : undefined}>
      {namesText}
    </span>
  ) as HTMLElement;
  const root = (
    <div class={classes.file}>
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
        {nameField}
      </button>
      <Show when={derive(names, (n) => n.length > 0)}>
        {() => (
          <button
            type="button"
            class={classes.fileClear}
            aria-label="Clear"
            data-tip="Clear"
            onClick={() => binding.set(null)}
          >
            <Icon name="x" />
          </button>
        )}
      </Show>
      {input}
    </div>
  ) as HTMLElement;
  watch([names, dragging], () => {
    const n = names.value;
    setAttr(root, 'data-dragging', dragging.value);
    setAttr(nameField, 'data-empty', n.length === 0);
    setText(
      namesText,
      n.length === 0
        ? multiple
          ? 'Select files…'
          : 'Select a file…'
        : n.length === 1
          ? n[0]
          : `${n.length} files · ${n.join(', ')}`,
    );
  });
  return root;
}

/** A plain row's control, labelled by the row's label (`labelId`; `label` its text, for the parts of a compound one). */
function control(param: PD.Any, binding: FieldBinding, labelId: string, label: string): Child {
  switch (param.type) {
    case 'boolean':
      return <Toggle binding={binding} labelId={labelId} />;
    case 'number':
      return <NumberControl param={param} binding={binding} labelId={labelId} />;
    case 'text':
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

/** Text: one line, or several (a description, a script: a field that grows by hand, Enter makes a new line). */
function TextControl(props: { binding: FieldBinding; labelId?: string; placeholder?: string; multiline?: boolean }) {
  const { binding } = props;
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
  watch([binding.value], () => setValue(el, binding.value.value));
  return el;
}

/** An input's value ('' for none), written only when it differs (a write moves the caret). */
function setValue(el: HTMLInputElement | HTMLTextAreaElement, value: unknown) {
  const v = value === null || value === undefined ? '' : String(value);
  if (el.value !== v) el.value = v;
}

function Toggle(props: { binding: FieldBinding<boolean>; labelId?: string }) {
  const { binding } = props;
  const el = (
    <button
      type="button"
      role="switch"
      class={classes.switch}
      aria-labelledby={props.labelId}
      onClick={() => binding.set(!binding.value.value)}
    >
      <span class={classes.knob} />
    </button>
  ) as HTMLElement;
  watch([binding.value], () => setAttr(el, 'aria-checked', String(!!binding.value.value)));
  return el;
}

function NumberControl(props: { param: PD.Numeric; binding: FieldBinding<number>; labelId?: string }) {
  const { param, binding, labelId } = props;
  if (param.min === undefined || param.max === undefined) {
    // Text, not a native number field: that one shows decimals in the browser's locale ("0,5") while every other
    // number reads "0.5"; a comma typed is read as the decimal point.
    const input = (
      <input
        class={`${classes.input} ${classes.numeric}`}
        type="text"
        inputMode="decimal"
        aria-labelledby={labelId}
        onChange={(e) => {
          const v = parseNumber(e.currentTarget.value);
          if (Number.isFinite(v)) binding.set(v);
          e.currentTarget.value = format(Number.isFinite(v) ? v : binding.value.value, param.step);
        }}
      />
    ) as HTMLInputElement;
    watch([binding.value], () => setValue(input, format(binding.value.value, param.step)));
    return input;
  }
  const min = param.min;
  const max = param.max;
  // A drag in the frame lane (dom §3): input events come several a frame, the field is set once a frame (it may update
  // a structure's representation).
  const dragged = frameState(Number.NaN, (a, b) => Object.is(a, b));
  effect(dragged, (v) => {
    if (!Number.isNaN(v)) binding.set(v);
  });
  const range = (
    <input
      class={classes.range}
      type="range"
      aria-labelledby={labelId}
      min={String(min)}
      max={String(max)}
      step={String(param.step ?? (max - min) / 100)}
      onInput={(e) => dragged.offer(Number(e.currentTarget.value))}
    />
  ) as HTMLInputElement;
  const number = (
    <input
      class={classes.number}
      type="text"
      inputMode="decimal"
      aria-labelledby={labelId}
      onChange={(e) => {
        const v = parseNumber(e.currentTarget.value);
        const shown = Number.isFinite(v) ? Math.min(max, Math.max(min, v)) : binding.value.value;
        if (Number.isFinite(v)) binding.set(shown);
        // The box says what holds, whether or not the value changed (a clamp to the same value, text that is no number).
        e.currentTarget.value = format(shown, param.step);
      }}
    />
  ) as HTMLInputElement;
  watch([binding.value], () => {
    const v = binding.value.value;
    setValue(range, v);
    // The filled part of the track (no native accent fill: a setting is not the action).
    setStyle(range, '--msui-params-progress', `${(100 * (v - min)) / (max - min || 1)}%`);
    setValue(number, format(v, param.step));
  });
  return (
    <div class={classes.slider}>
      {range}
      {number}
    </div>
  );
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
 * A choice as a row: the label and the current option in a select's field; opened, the options over the field (the
 * current one on its text; grouped by Mol*'s categories). Values may be any type: the list works with indices.
 */
function SelectRow(props: {
  label: string;
  head: Child;
  options: readonly PD.SelectOption<any>[];
  binding: FieldBinding;
  /** After the select, in the row (a choice's params toggle). */
  after?: Child;
}) {
  const { options, binding } = props;
  const open = state(false);
  const index = derive(binding.value, (v) => String(options.findIndex((o) => o[0] === v)));
  // Mol*'s option categories are the list's groups.
  const indexed = options.map(([, label, group], i) => [String(i), label, group] as const);
  // The button's text: a node of the row's own, written by the row's watch.
  const labelText = document.createTextNode('');
  const button = (<SelectButton size="sm" open={open} aria-label={props.label} label={labelText} />) as HTMLElement;
  watch([index], () => setText(labelText, indexed[Number(index.value)]?.[1] ?? ''));
  return (
    <>
      <Row head={props.head}>
        {button}
        {props.after}
      </Row>
      <Show when={open}>
        {() => (
          <Floating anchor={button} placement="over">
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
          </Floating>
        )}
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
}): Child {
  const { param, binding, filter, place } = props;
  // As v5's: a group with nothing to show has no row.
  if (!hasRows(param.params)) return null;
  const presets = param.presets;
  // A fold's rows are a level deeper; a flat group's at its own level.
  const fields = (params: PD.Params, deeper: boolean) => (
    <GroupFields params={params} binding={binding} filter={filter} place={at(place, place.path, deeper)} />
  );
  // Its fold, open as the param says at first, then as the user leaves it (for the form's life).
  const folded = place.form.fold(`${place.path}#open`, !!param.isExpanded);
  // A preset sets the fields it names, over the group's current values (Mol*'s presets name them all).
  const apply = (value: unknown) => binding.set({ ...(binding.value.value ?? {}), ...(value as object) });
  const current = presets ? derive(binding.value, (v) => presetIndex(presets, v, true)) : undefined;
  // A flat group has no header of its own: its params are rows like the others, its presets a row before them.
  if (param.isFlat)
    return (
      <>
        {presets && current ? (
          <PresetsRow label={props.label} presets={presets} current={current} onPick={apply} />
        ) : null}
        {fields(param.params, false)}
      </>
    );
  // A pivot (Mol*'s: a provider's id) is the group's main field: a row of its own, the rest folded under the group.
  const pivot = typeof param.pivot === 'string' && param.params[param.pivot] ? param.pivot : undefined;
  if (pivot && (!presets || !current)) {
    const { [pivot]: main, ...rest } = param.params;
    return (
      <>
        {fields({ [pivot]: main }, false)}
        {hasRows(rest) ? (
          <Collapsible title={props.label} state={folded} force={props.force} level={place.level}>
            {() => fields(rest, true)}
          </Collapsible>
        ) : null}
      </>
    );
  }
  if (!presets || !current)
    return (
      <Collapsible title={props.label} state={folded} force={props.force} level={place.level}>
        {() => fields(param.params, true)}
      </Collapsible>
    );
  // With presets: Presets at the header's end (as Mol*'s), the matching one named on the value column; a pick opens the
  // group so what it set shows.
  return (
    <Collapsible
      title={props.label}
      state={folded}
      force={props.force}
      level={place.level}
      meta={derive(current, (i) => (i === undefined ? 'Custom' : presets[i][1]))}
      action={
        <PresetsButton
          label={props.label}
          presets={() => presets}
          current={current}
          onPick={(v) => {
            apply(v);
            folded.set(true);
          }}
        />
      }
    >
      {() => fields(param.params, true)}
    </Collapsible>
  );
}

/** A flat group's presets: a row like a select's, reading the matching preset (or Custom); the presets over it. */
function PresetsRow(props: {
  label: string;
  presets: readonly PD.SelectOption<unknown>[];
  current: ReadonlyState<number | undefined>;
  onPick(value: unknown): void;
}) {
  const open = state(false);
  const labelText = document.createTextNode('');
  const button = (
    <SelectButton size="sm" open={open} aria-label={`${props.label} presets`} label={labelText} />
  ) as HTMLElement;
  watch([props.current], () => {
    const i = props.current.value;
    setText(labelText, i === undefined ? 'Custom' : props.presets[i][1]);
  });
  return (
    <>
      <Row
        head={
          <span class={classes.labelText} data-tip-overflow="">
            Presets
          </span>
        }
      >
        {button}
      </Row>
      <PresetsList
        open={open}
        anchor={button}
        label={`${props.label} presets`}
        presets={() => props.presets}
        current={props.current}
        onPick={props.onPick}
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
 * An on/off mapped param: a switch, and while on, its own settings behind a sliders button after it, opening as rows a
 * level deeper under the row (occlusion's samples, radius, bias…).
 */
function OnOffControl(props: {
  label: string;
  labelId: string;
  head: Child;
  param: PD.Mapped<any>;
  binding: FieldBinding;
  place: Place;
}) {
  const { param, binding, place } = props;
  const on = derive(binding.value, (v) => v?.name === 'on');
  const open = place.form.fold(`${place.path}#params`);
  const toggle: FieldBinding<boolean> = {
    value: on,
    set: (v) => {
      binding.set({ name: v ? 'on' : 'off', params: defaults(param.map(v ? 'on' : 'off')) });
      if (!v) open.set(false);
    },
  };
  const inner = param.map('on');
  const hasSettings = inner.type === 'group' && hasRows(inner.params);
  const paramsBinding: FieldBinding = {
    value: derive(binding.value, (v) => v?.params),
    set: (p) => binding.set({ name: 'on', params: p }),
  };
  return (
    <>
      <Row head={props.head}>
        <Toggle binding={toggle} labelId={props.labelId} />
        {hasSettings ? (
          <Show when={on}>
            {() => <ExpandButton open={open} label={`${props.label} settings`} onClick={() => open.set(!open.value)} />}
          </Show>
        ) : null}
      </Row>
      {hasSettings ? (
        <Show when={derive([on, open], (a, b) => a && b)} reveal>
          {() => (
            <GroupFields
              params={(inner as PD.Group<any>).params}
              binding={paramsBinding}
              place={at(place, `${place.path}:on`, true)}
            />
          )}
        </Show>
      ) : null}
    </>
  );
}

/**
 * A mapped param (PD.Mapped, MappedStatic), as v5's `MappedControl`: the choice; when the chosen case's params are a
 * group, a sliders button after it opens them a level deeper (with the form's rules inside: essentials, then Advanced);
 * a flat group's params are rows at the same level, any other case's value one row under the choice. Its open state is
 * the form's (kept from one case to the next, as v5's).
 */
function MappedControl(props: {
  label: string;
  head: Child;
  param: PD.Mapped<any>;
  binding: FieldBinding;
  /** For the chosen case's rows. */
  filter?: ReadonlyState<string>;
  /** Open while searching. */
  force?: ReadonlyState<boolean>;
  place: Place;
}) {
  const { param, binding, place } = props;
  const name = derive(binding.value, (v) => v?.name as string | undefined);
  const nameBinding: FieldBinding = {
    value: name,
    set: (n) => binding.set({ name: n, params: defaults(param.map(n)) }),
  };
  /** Does the case have params behind the button (a group that is not flat, with a row to show)? */
  const toggled = (n: string | undefined) => {
    if (n === undefined) return false;
    const inner = param.map(n);
    return inner.type === 'group' && !inner.isFlat && hasRows(inner.params);
  };
  const hasToggle = derive(name, toggled);
  const open = place.form.fold(`${place.path}#params`);
  const shown = props.force ? derive([open, props.force], (o, f) => o || f) : open;
  const toggle = (
    <Show when={hasToggle}>
      {() => <ExpandButton open={shown} label={`${props.label} settings`} onClick={() => open.set(!open.value)} />}
    </Show>
  );
  return (
    <>
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
            set: (p) => binding.set({ name: n, params: p }),
          };
          const path = `${place.path}:${n}`;
          if (inner.type === 'group') {
            // A flat group (PD's isFlat) reads as more rows at the same level.
            if (inner.isFlat)
              return (
                <GroupFields
                  params={inner.params}
                  binding={paramsBinding}
                  filter={props.filter}
                  place={at(place, path)}
                />
              );
            if (!hasRows(inner.params)) return null;
            return (
              <Show when={shown} reveal>
                {() => (
                  <GroupFields
                    params={inner.params}
                    binding={paramsBinding}
                    filter={props.filter}
                    place={at(place, path, true)}
                  />
                )}
              </Show>
            );
          }
          // Another kind of value: its own row, named by the case (as v5's), a level deeper.
          return (
            <Nested level={place.level + 1}>
              <ParamRow
                name={n}
                param={inner}
                binding={paramsBinding}
                filter={props.filter}
                place={at(place, path, true)}
              />
            </Nested>
          );
        }}
      </Show>
    </>
  );
}

/**
 * The params of a group as rows (the form's rules: essentials, Advanced, categories, hideIf on the group's values),
 * bound to the group's value object, at their place's level.
 */
function GroupFields(props: {
  params: PD.Params;
  binding: FieldBinding;
  filter?: ReadonlyState<string>;
  place: Place;
}) {
  const { params, binding } = props;
  return (
    <Nested level={props.place.level}>
      <ParamControls
        params={params}
        field={(k) => keyBinding(binding, k)}
        visible={(k) => {
          const hideIf = params[k].hideIf;
          if (!hideIf) return undefined;
          return derive(binding.value, (v) => !hideIf(v ?? {}));
        }}
        filter={props.filter}
        place={props.place}
      />
    </Nested>
  );
}

/** Conditioned: the condition as a choice, then the value's own control for it, a level deeper. */
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
    set: (c) => binding.set(param.conditionedValue(binding.value.value, c)),
  };
  return (
    <>
      <SelectRow label={props.label} head={props.head} options={param.select.options} binding={conditionBinding} />
      <Show when={condition}>
        {(c) => (
          <Nested level={place.level + 1}>
            <ParamRow
              name={c}
              param={param.conditionParams[c]}
              binding={binding}
              place={at(place, `${place.path}:${c}`, true)}
            />
          </Nested>
        )}
      </Show>
    </>
  );
}

const ScriptLanguages: PD.SelectOption<string>[] = [
  ['mol-script', 'MolScript'],
  ['pymol', 'PyMOL'],
  ['vmd', 'VMD'],
  ['jmol', 'Jmol'],
];

/** A script: its language as a choice, its expression under it (several lines), a level deeper. */
function ScriptControl(props: { head: Child; label: string; binding: FieldBinding; place: Place }) {
  const { binding } = props;
  const language = keyBinding(binding, 'language');
  const expression = keyBinding(binding, 'expression');
  const expressionId = uniqueId('msui-field');
  const area = (
    <textarea
      class={classes.textarea}
      aria-labelledby={expressionId}
      rows={3}
      spellcheck={false}
      onChange={(e) => expression.set(e.currentTarget.value)}
    />
  ) as HTMLTextAreaElement;
  watch([expression.value], () => setValue(area, expression.value.value));
  return (
    <>
      <SelectRow label={props.label} head={props.head} options={ScriptLanguages} binding={language} />
      <Nested level={props.place.level + 1}>
        <Row
          head={
            <span class={classes.labelText} id={expressionId} data-tip-overflow="">
              Expression
            </span>
          }
        >
          {area}
        </Row>
      </Nested>
    </>
  );
}

/**
 * An object list (a list of like things: a volume's channels, a set of points): the count on the value column, Add and
 * Presets at the row's end; the items under it a level deeper, each a row of its label (edit and remove on hover; its
 * form unfolds under it, a level deeper again).
 */
function ObjectListRow(props: {
  head: Child;
  label: string;
  param: PD.ObjectList<any>;
  binding: FieldBinding<unknown[] | undefined>;
  place: Place;
}) {
  const { param, binding, place } = props;
  const presets = param.presets;
  const items = derive(binding.value, (v) => v ?? []);
  const replace = (i: number, item: unknown) => {
    const list = [...(binding.value.value ?? [])];
    list[i] = item;
    binding.set(list);
  };
  const count = document.createTextNode('');
  watch([items], () => {
    const n = items.value.length;
    setText(count, n === 1 ? '1 item' : `${n} items`);
  });
  return (
    <>
      <Row head={props.head}>
        <span class={classes.count}>{count}</span>
        <Button
          variant="ghost"
          size="sm"
          title={`Add to ${props.label}`}
          icon={<Icon name="plus" />}
          onClick={() => binding.set([...(binding.value.value ?? []), param.ctor()])}
        />
        {presets ? (
          <PresetsButton
            label={props.label}
            presets={() => presets}
            current={derive(binding.value, (v) => presetIndex(presets, v))}
            onPick={(v) => binding.set(v)}
          />
        ) : null}
      </Row>
      <Nested level={place.level + 1}>
        <For each={derive(items, (l) => l.map((_, i) => i))}>
          {(i) => (
            <ObjectListItem param={param} binding={binding} index={i} replace={replace} level={place.level + 1} />
          )}
        </For>
      </Nested>
    </>
  );
}

function ObjectListItem(props: {
  param: PD.ObjectList<any>;
  binding: FieldBinding<unknown[] | undefined>;
  index: number;
  replace: (i: number, item: unknown) => void;
  level: number;
}) {
  const { param, binding, index } = props;
  const open = state(false);
  const labelText = document.createTextNode('');
  const edit = (
    <Button
      variant="ghost"
      size="sm"
      class={classes.itemAction}
      title="Edit"
      icon={<Icon name="sliders-horizontal" />}
      onClick={(e) => {
        e.stopPropagation();
        open.set(!open.value);
      }}
    />
  ) as HTMLElement;
  const root = (
    <div class={classes.item}>
      {/* A click anywhere on the row opens or closes the item's params (the sliders say so; remove does not). */}
      <div class={classes.itemRow} onClick={() => open.set(!open.value)}>
        <span class={classes.itemLabel}>{labelText}</span>
        {edit}
        <Button
          variant="ghost"
          size="sm"
          class={classes.itemAction}
          title="Remove"
          icon={<Icon name="x" />}
          onClick={(e) => {
            e.stopPropagation();
            binding.set((binding.value.value ?? []).filter((_, j) => j !== index));
          }}
        />
      </div>
      <Show when={open} reveal>
        {() => {
          // The item's own form: each change replaces the item in the list (a new list: the state stores it as is).
          const form = new ParamsFormModel(param.element, binding.value.value?.[index] as any, {
            onChange: (v) => props.replace(index, v),
          });
          return (
            <Nested level={props.level + 1}>
              <FormRows form={form} level={props.level + 1} />
            </Nested>
          );
        }}
      </Show>
    </div>
  ) as HTMLElement;
  // The item's label only: a change of another item in the list wakes this watch no further.
  const label = derive(binding.value, (v) => (v?.[index] !== undefined ? param.getLabel(v[index]) : ''));
  watch([label, open], () => {
    setText(labelText, label.value);
    setAttr(root, 'data-open', open.value);
    setAttr(edit, 'aria-expanded', String(open.value));
  });
  return root;
}

/** The sliders button that opens a choice's or a switch's settings: `aria-expanded` follows `open`. */
function ExpandButton(props: { open: ReadonlyState<boolean>; label: string; onClick(): void }) {
  const el = (
    <Button
      variant="ghost"
      size="sm"
      class={classes.expand}
      title={props.label}
      icon={<Icon name="sliders-horizontal" />}
      onClick={props.onClick}
    />
  ) as HTMLElement;
  watch([props.open], () => setAttr(el, 'aria-expanded', String(props.open.value)));
  return el;
}

function defaults(param: PD.Any): any {
  if (param.type === 'group') {
    const out: Record<string, unknown> = {};
    for (const k of Object.keys(param.params)) out[k] = param.params[k].defaultValue;
    return out;
  }
  return param.defaultValue;
}

/**
 * A fold of the form (Advanced, a category, a group): its header a row on the form's columns (the title on the label
 * column at its level, then a chevron; on the value column the preset its values match and an action at the end); its
 * rows, rendered while open, a level deeper.
 */
function Collapsible(props: {
  title: string;
  /** Its open state (the form's folds). */
  state: State<boolean>;
  /** Forced open (a search: its matches are inside); the user's choice again once it clears. */
  force?: ReadonlyState<boolean>;
  level: number;
  /** A value after the title (the preset the values match). */
  meta?: ReadonlyState<string>;
  /** A control at the header's end (Presets). */
  action?: Child;
  children: () => Child;
}) {
  const own = props.state;
  const open = props.force ? derive([own, props.force], (o, f) => o || f) : own;
  const meta = props.meta;
  const metaText = document.createTextNode('');
  const header = (
    <button type="button" class={classes.foldHeader} onClick={() => own.set(!own.value)}>
      <span class={classes.labelText} data-tip-overflow="">
        {props.title}
      </span>
      <Icon name="chevron-right" class={classes.chevron} />
    </button>
  ) as HTMLElement;
  watch(meta ? [open, meta] : [open], () => {
    setAttr(header, 'aria-expanded', String(open.value));
    if (meta) setText(metaText, meta.value);
  });
  return (
    <div class={classes.fold}>
      <div class={classes.row}>
        <span class={classes.label}>{header}</span>
        {meta || props.action ? (
          <div class={classes.control}>
            {meta ? <span class={classes.meta}>{metaText}</span> : null}
            {props.action}
          </div>
        ) : null}
      </div>
      <Show when={open} reveal>
        {() => <Nested level={props.level + 1}>{props.children()}</Nested>}
      </Show>
    </div>
  );
}
