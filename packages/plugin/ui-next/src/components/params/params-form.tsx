/** @jsxImportSource @molstar/plugin-ui-next/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { Asset } from '@molstar/core/util/assets';
import { camelCaseToWords } from '@molstar/core/util/string';
import { combineLatest, distinctUntilChanged, map, type Observable, of } from 'rxjs';
import { For, Show } from '../../dom/control-flow.js';
import type { Child } from '../../dom/jsx-types.js';
import { ParamsFormModel } from '../../forms/params-form.js';
import { observableState } from '../../reactive/observable-state.js';
import { OptionList, SelectButton } from '../select/select.js';
import { ColorPicker, colorToHex } from '../color/color-picker.js';
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

/** A field's value as a control sees it: a stream, the current value, and how to change it. */
export type { FieldBinding };

/**
 * The form for a `ParamsFormModel`: essentials first, the rest under "Advanced" (as the classic UI). Each row binds
 * only to its own field; groups and mapped params nest.
 */
export function ParamsForm(props: {
  form: ParamsFormModel<any>;
  /**
   * Shows only the rows whose label or description contains this text (or that hold such rows); groups open while it
   * is set. A search over settings.
   */
  filter?: Observable<string>;
}) {
  const { form } = props;
  const filter = normalized(props.filter);
  const { essentials, advanced } = form.layout();
  const row = (name: string) => (
    <ParamRow
      name={name}
      param={form.params[name]}
      visible={form.visible(name)}
      filter={filter}
      binding={{ value: form.field(name), get: () => form.values.value[name], set: (v) => form.set(name, v) }}
    />
  );
  const advancedShown = filter?.pipe(
    map((q) => advanced.some((name) => matches(name, form.params[name], q))),
    distinctUntilChanged(),
  );
  return (
    <div class={classes.form}>
      {essentials.map(row)}
      {advanced.length ? (
        <div hidden={advancedShown?.pipe(map((v) => !v))}>
          <Collapsible title="Advanced" force={searching(filter)}>
            {() => advanced.map(row)}
          </Collapsible>
        </div>
      ) : null}
    </div>
  );
}

/** The filter text, trimmed and lower-cased; undefined when there is no filter. */
function normalized(filter?: Observable<string>): Observable<string> | undefined {
  return filter?.pipe(
    map((q) => q.trim().toLowerCase()),
    distinctUntilChanged(),
  );
}

function searching(filter?: Observable<string>): Observable<boolean> | undefined {
  return filter?.pipe(
    map((q) => q.length > 0),
    distinctUntilChanged(),
  );
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
  visible?: Observable<boolean>;
  filter?: Observable<string>;
}) {
  const { name, param, binding, filter } = props;
  const label = labelOf(name, param);
  const shown = filter
    ? combineLatest([props.visible ?? of(true), filter]).pipe(map(([v, q]) => v && matches(name, param, q)))
    : props.visible;
  const hidden = shown?.pipe(
    map((v) => !v),
    distinctUntilChanged(),
  );
  if (param.type === 'group') {
    // A group whose own name matches shows all of its rows.
    const inner = filter?.pipe(map((q) => (label.toLowerCase().includes(q) ? '' : q)));
    return (
      <div hidden={hidden}>
        <GroupControl label={label} param={param} binding={binding} filter={inner} force={searching(filter)} />
      </div>
    );
  }
  // Its help (a description, a legend; for a theme choice, the chosen theme's): a ? by the label, opening under the row.
  const { head, block } = withHelp(label, helpOf(param, binding.value));
  const wrap = (row: Child) => (
    <div hidden={hidden}>
      {row}
      {block}
    </div>
  );
  if (param.type === 'mapped') {
    return wrap(
      isOnOff(param) ? (
        <OnOffControl label={label} head={head} description={param.description} param={param} binding={binding} />
      ) : (
        <MappedControl label={label} head={head} param={param} binding={binding} />
      ),
    );
  }
  if (param.type === 'select') {
    return wrap(
      <SelectRow label={label} head={head} description={param.description} options={param.options} binding={binding} />,
    );
  }
  if (param.type === 'color') {
    return wrap(<ColorRow label={label} head={head} description={param.description} binding={binding} />);
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
    return wrap(<ConditionedControl head={head} label={label} param={param} binding={binding} />);
  if (param.type === 'script') return wrap(<ScriptControl head={head} label={label} binding={binding} />);
  // Converted: the control of what it converts to, the value converted both ways.
  if (param.type === 'converted') {
    const converted: FieldBinding = {
      value: binding.value.pipe(map((v) => param.fromValue(v))),
      get: () => param.fromValue(binding.get()),
      set: (v) => binding.set(param.toValue(v)),
    };
    const inner = { ...param.converted, label: param.converted.label ?? label } as PD.Any;
    return <ParamRow name={name} param={inner} binding={converted} visible={props.visible} filter={filter} />;
  }
  // A plain value (data carried along, not edited): no row, as in the classic UI.
  if (param.type === 'value') return null;
  return wrap(
    <div class={classes.row} title={param.description ?? label}>
      <span class={classes.label}>{head}</span>
      <div class={classes.control}>{control(param, binding)}</div>
    </div>,
  );
}

/** A color: the swatch and hex as a button; the picker opens under the row, full width. */
function ColorRow(props: { label: string; head?: Child; description?: string; binding: FieldBinding }) {
  const { binding } = props;
  const open = observableState(false);
  const hex = binding.value.pipe(map((v: number) => colorToHex(v)));
  return (
    <>
      <div class={classes.row} title={props.description}>
        <span class={classes.label}>{props.head ?? props.label}</span>
        <div class={classes.control}>
          <button type="button" class={classes.color} aria-expanded={open} onClick={() => open.set(!open.value)}>
            <span class={classes.swatch} style={{ background: hex }} />
            <span class={classes.hex}>{hex}</span>
          </button>
        </div>
      </div>
      <Show when={open}>{() => <ColorPicker value={binding.value} onChange={(c) => binding.set(c)} />}</Show>
    </>
  );
}

/**
 * Files (PD.File, PD.FileList): a field like a dropdown's, "Select files…" until some are chosen, then their names; a
 * click picks, files dropped on it are taken too, ✕ clears.
 */
function FileControl(props: { param: PD.FileParam | PD.FileListParam; binding: FieldBinding }) {
  const { param, binding } = props;
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
  const names = binding.value.pipe(
    map((v: Asset.File | Asset.File[] | null) => (v ? (Array.isArray(v) ? v : [v]) : []).map((f) => f.name)),
  );
  return (
    <div class={classes.file} data-dragging={dragging.asReadonly.pipe(map((d) => (d ? '' : undefined)))}>
      <button
        type="button"
        class={classes.filePick}
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
        <span class={classes.fileNames} data-empty={names.pipe(map((n) => (n.length ? undefined : '')))}>
          {names.pipe(
            map((n) =>
              n.length === 0
                ? multiple
                  ? 'Select files…'
                  : 'Select a file…'
                : n.length === 1
                  ? n[0]
                  : `${n.length} files · ${n.join(', ')}`,
            ),
          )}
        </span>
      </button>
      <Show when={names.pipe(map((n) => n.length > 0))}>
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

function control(param: PD.Any, binding: FieldBinding): Child {
  switch (param.type) {
    case 'boolean':
      return <Toggle binding={binding} />;
    case 'number':
      return <NumberControl param={param} binding={binding} />;
    case 'text':
      // Several lines (a description, a script): a field that grows by hand, Enter makes a new line.
      if (param.multiline) {
        return (
          <textarea
            class={classes.textarea}
            rows={3}
            value={binding.value}
            placeholder={param.placeholder}
            onChange={(e) => binding.set(e.currentTarget.value)}
          />
        );
      }
      return (
        <input
          class={classes.input}
          value={binding.value}
          placeholder={param.placeholder}
          onChange={(e) => binding.set(e.currentTarget.value)}
        />
      );
    case 'file':
    case 'file-list':
      return <FileControl param={param} binding={binding} />;
    case 'vec3':
      return <Vec3Control param={param} binding={binding} />;
    case 'interval':
      return <IntervalControl param={param} binding={binding} />;
    case 'url':
      return <UrlControl binding={binding} />;
    default:
      return <span class={classes.unsupported}>{param.type}</span>;
  }
}

function Toggle(props: { binding: FieldBinding }) {
  const { binding } = props;
  return (
    <button
      type="button"
      role="switch"
      class={classes.switch}
      aria-checked={binding.value.pipe(map((v) => !!v))}
      onClick={() => binding.set(!binding.get())}
    >
      <span class={classes.knob} />
    </button>
  );
}

function NumberControl(props: { param: PD.Numeric; binding: FieldBinding }) {
  const { param, binding } = props;
  const text = binding.value.pipe(map((v: number) => format(v, param.step)));
  if (param.min === undefined || param.max === undefined) {
    return (
      <input
        class={classes.input}
        type="number"
        step={String(param.step ?? 'any')}
        value={text}
        onChange={(e) => {
          const v = Number(e.currentTarget.value);
          if (Number.isFinite(v)) binding.set(v);
        }}
      />
    );
  }
  return (
    <div class={classes.slider}>
      <input
        class={classes.range}
        type="range"
        min={String(param.min)}
        max={String(param.max)}
        step={String(param.step ?? (param.max - param.min) / 100)}
        value={binding.value.pipe(map(String))}
        style={{
          // The filled part of the track, for the quiet styling (no native accent fill).
          '--msp-progress': binding.value.pipe(
            map((v: number) => `${(100 * (v - param.min!)) / (param.max! - param.min! || 1)}%`),
          ),
        }}
        onInput={(e) => binding.set(Number(e.currentTarget.value))}
      />
      <input
        class={classes.number}
        value={text}
        onChange={(e) => {
          const v = Number(e.currentTarget.value);
          if (Number.isFinite(v)) binding.set(Math.min(param.max!, Math.max(param.min!, v)));
        }}
      />
    </div>
  );
}

function format(v: number, step?: number) {
  if (!Number.isFinite(v)) return '';
  const digits = step && step < 1 ? Math.min(4, Math.ceil(-Math.log10(step))) : 0;
  return v.toFixed(digits);
}

/**
 * A choice as a row: the label and the current option; opened, the options expand inline under the whole row (Mol*'s
 * action menu). Values may be any type: the list works with indices.
 */
function SelectRow(props: {
  label: string;
  /** The label's content (with its help); the label text without. */
  head?: Child;
  description?: string;
  options: readonly PD.SelectOption<any>[];
  binding: FieldBinding;
}) {
  const { options, binding } = props;
  const open = observableState(false);
  const index = binding.value.pipe(map((v) => String(options.findIndex((o) => o[0] === v))));
  const indexed = options.map(([, label], i) => [String(i), label] as const);
  return (
    <>
      <div class={classes.row} title={props.description}>
        <span class={classes.label}>{props.head ?? props.label}</span>
        <div class={classes.control}>
          <SelectButton
            size="sm"
            open={open}
            aria-label={props.label}
            label={index.pipe(map((i) => indexed[Number(i)]?.[1] ?? ''))}
          />
        </div>
      </div>
      <Show when={open}>
        {() => (
          <OptionList
            options={indexed}
            value={index}
            onPick={(i) => {
              open.set(false);
              binding.set(options[Number(i)][0]);
            }}
            onClose={() => open.set(false)}
          />
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
  filter?: Observable<string>;
  /** Open while searching. */
  force?: Observable<boolean>;
}) {
  const { param, binding, filter } = props;
  // A flat group has no header of its own: its params are rows like the others.
  if (param.isFlat) return <Fields params={param.params} binding={binding} filter={filter} />;
  return (
    <Collapsible title={props.label} open={!!param.isExpanded} force={props.force}>
      {() => <Fields params={param.params} binding={binding} filter={filter} />}
    </Collapsible>
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
  head?: Child;
  description?: string;
  param: PD.Mapped<any>;
  binding: FieldBinding;
}) {
  const { param, binding } = props;
  const on = binding.value.pipe(
    map((v) => v?.name === 'on'),
    distinctUntilChanged(),
  );
  const open = observableState(false);
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
  const hasSettings = inner.type === 'group' && Object.keys(inner.params).length > 0;
  const paramsBinding: FieldBinding = {
    value: binding.value.pipe(
      map((v) => v?.params),
      distinctUntilChanged(),
    ),
    get: () => binding.get()?.params,
    set: (p) => binding.set({ name: 'on', params: p }),
  };
  return (
    <>
      <div class={classes.row} title={props.description}>
        <span class={classes.label}>{props.head ?? props.label}</span>
        <div class={classes.control}>
          <Toggle binding={toggle} />
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
        <Show when={combineLatest([on, open.asReadonly]).pipe(map(([a, b]) => a && b))}>
          {() => (
            <div class={classes.subSettings}>
              <Fields params={(inner as PD.Group<any>).params} binding={paramsBinding} />
            </div>
          )}
        </Show>
      ) : null}
    </>
  );
}

/** A mapped param: the choice, then the chosen option's own params (rendered again when the choice changes). */
function MappedControl(props: { label: string; head?: Child; param: PD.Mapped<any>; binding: FieldBinding }) {
  const { param, binding } = props;
  const name = binding.value.pipe(
    map((v) => v?.name),
    distinctUntilChanged(),
  );
  const nameBinding: FieldBinding = {
    value: name,
    get: () => binding.get()?.name,
    set: (n) => {
      const params = param.map(n);
      binding.set({ name: n, params: defaults(params) });
    },
  };
  return (
    <div class={classes.mapped}>
      <SelectRow label={props.label} head={props.head} options={param.select.options} binding={nameBinding} />
      <Show when={name}>
        {(n) => {
          const inner = param.map(n);
          const paramsBinding: FieldBinding = {
            value: binding.value.pipe(
              map((v) => v?.params),
              distinctUntilChanged(),
            ),
            get: () => binding.get()?.params,
            set: (p) => binding.set({ name: n, params: p }),
          };
          // A flat group (PD's isFlat) reads as more rows at the same level; others nest under the choice.
          if (inner.type === 'group')
            return <Fields params={inner.params} binding={paramsBinding} nested={!inner.isFlat} />;
          return <ParamRow name="value" param={inner} binding={paramsBinding} />;
        }}
      </Show>
    </div>
  );
}

/** The params of a group as rows, bound to the group's value object. */
function Fields(props: { params: PD.Params; binding: FieldBinding; nested?: boolean; filter?: Observable<string> }) {
  const { params, binding } = props;
  return (
    <div class={props.nested ? classes.nested : classes.fields}>
      {Object.keys(params)
        .filter((k) => !params[k].isHidden)
        .map((k) => {
          const param = params[k];
          const visible = param.hideIf ? binding.value.pipe(map((v) => !param.hideIf!(v ?? {}))) : undefined;
          return (
            <ParamRow
              name={k}
              param={param}
              visible={visible}
              filter={props.filter}
              binding={{
                value: binding.value.pipe(
                  map((v) => v?.[k]),
                  distinctUntilChanged(),
                ),
                get: () => binding.get()?.[k],
                set: (cv) => binding.set({ ...binding.get(), [k]: cv }),
              }}
            />
          );
        })}
    </div>
  );
}

/** Conditioned: the condition as a choice, then the value's own control for it (nested under the choice). */
function ConditionedControl(props: {
  head: Child;
  label: string;
  param: PD.Conditioned<any, any, any>;
  binding: FieldBinding;
}) {
  const { param, binding } = props;
  const condition = binding.value.pipe(
    map((v) => String(param.conditionForValue(v))),
    distinctUntilChanged(),
  );
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
            <ParamRow name={c} param={param.conditionParams[c]} binding={binding} />
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
    value: binding.value.pipe(
      map((v) => v?.[k] ?? ''),
      distinctUntilChanged(),
    ),
    get: () => binding.get()?.[k],
    set: (x) => binding.set({ ...binding.get(), [k]: x }),
  });
  const expression = part('expression');
  return (
    <div class={classes.mapped}>
      <SelectRow label={props.label} head={props.head} options={ScriptLanguages} binding={part('language')} />
      <div class={classes.nested}>
        <div class={classes.row}>
          <span class={classes.label}>Expression</span>
          <div class={classes.control}>
            <textarea
              class={classes.textarea}
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
  const items = binding.value.pipe(map((v: unknown[] | undefined) => v ?? []));
  const replace = (i: number, item: unknown) => {
    const list = [...(binding.get() ?? [])];
    list[i] = item;
    binding.set(list);
  };
  return (
    <>
      <div class={classes.row} title={param.description}>
        <span class={classes.label}>{props.head}</span>
        <div class={classes.control}>
          <span class={classes.count}>{items.pipe(map((l) => (l.length === 1 ? '1 item' : `${l.length} items`)))}</span>
        </div>
      </div>
      {/* The items in one quiet list under the row; Add ends it. */}
      <div class={classes.items}>
        <For each={items.pipe(map((l) => l.map((_, i) => i)))}>
          {(i) => <ObjectListItem param={param} binding={binding} index={i} replace={replace} />}
        </For>
        <button
          type="button"
          class={classes.itemAdd}
          aria-label={`Add to ${props.label}`}
          onClick={() => binding.set([...(binding.get() ?? []), param.ctor()])}
        >
          Add
        </button>
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
    <div class={classes.item} data-open={open.asReadonly.pipe(map((o) => (o ? '' : undefined)))}>
      <div class={classes.itemRow}>
        <span class={classes.itemLabel}>
          {binding.value.pipe(map((v: unknown[]) => (v?.[index] !== undefined ? param.getLabel(v[index]) : '')))}
        </span>
        <button
          type="button"
          class={classes.settings}
          title="Edit"
          aria-label="Edit"
          aria-expanded={open}
          onClick={() => open.set(!open.value)}
        >
          <Icon name="sliders-horizontal" size={14} />
        </button>
        <button
          type="button"
          class={classes.settings}
          title="Remove"
          aria-label="Remove"
          onClick={() => binding.set((binding.get() ?? []).filter((_: unknown, j: number) => j !== index))}
        >
          <Icon name="x" size={13} />
        </button>
      </div>
      <Show when={open}>
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
function Collapsible(props: { title: string; open?: boolean; force?: Observable<boolean>; children: () => Child }) {
  const own = observableState(!!props.open);
  // Forced open while a filter is set (its matches are inside); the user's choice again once it clears.
  const open = props.force
    ? combineLatest([own.asReadonly, props.force]).pipe(map(([o, f]) => o || f))
    : own.asReadonly;
  return (
    <div class={classes.collapsible}>
      <button type="button" class={classes.collapsibleHeader} aria-expanded={open} onClick={() => own.set(!own.value)}>
        {/* The name first, aligned with the labels; the chevron after it. */}
        {props.title}
        <Icon name="chevron-right" size={12} class={classes.chevron} />
      </button>
      <Show when={open}>{() => <div class={classes.nested}>{props.children()}</div>}</Show>
    </div>
  );
}
