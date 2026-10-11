/** @jsxImportSource @molstar/plugin-ui-next/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { camelCaseToWords } from '@molstar/core/util/string';
import { combineLatest, distinctUntilChanged, map, type Observable } from 'rxjs';
import { Show } from '../../dom/control-flow.js';
import type { Child } from '../../dom/jsx-types.js';
import type { ParamsFormModel } from '../../forms/params-form.js';
import { observableState } from '../../reactive/observable-state.js';
import { OptionList, SelectButton } from '../select/select.js';
import { ColorPicker, colorToHex } from '../color/color-picker.js';
import { Icon } from '../icon/icon.js';
import { classes } from './_generated/params-form.module.css.js';

/** A field's value as a control sees it: a stream, the current value, and how to change it. */
interface FieldBinding {
  value: Observable<any>;
  get(): any;
  set(value: any): void;
}

/**
 * The form for a `ParamsFormModel`: essentials first, the rest under "Advanced" (as the classic UI). Each row binds
 * only to its own field; groups and mapped params nest.
 */
export function ParamsForm(props: { form: ParamsFormModel<any> }) {
  const { form } = props;
  const { essentials, advanced } = form.layout();
  const row = (name: string) => (
    <ParamRow
      name={name}
      param={form.params[name]}
      visible={form.visible(name)}
      binding={{ value: form.field(name), get: () => form.values.value[name], set: (v) => form.set(name, v) }}
    />
  );
  return (
    <div class={classes.form}>
      {essentials.map(row)}
      {advanced.length ? <Collapsible title="Advanced">{() => advanced.map(row)}</Collapsible> : null}
    </div>
  );
}

function ParamRow(props: { name: string; param: PD.Any; binding: FieldBinding; visible?: Observable<boolean> }) {
  const { name, param, binding } = props;
  const label = param.label || camelCaseToWords(name);
  const hidden = props.visible?.pipe(map((v) => !v));
  if (param.type === 'group') {
    return (
      <div hidden={hidden}>
        <GroupControl label={label} param={param} binding={binding} />
      </div>
    );
  }
  if (param.type === 'mapped') {
    return (
      <div hidden={hidden}>
        {isOnOff(param) ? (
          <OnOffControl label={label} description={param.description} param={param} binding={binding} />
        ) : (
          <MappedControl label={label} param={param} binding={binding} />
        )}
      </div>
    );
  }
  if (param.type === 'select') {
    return (
      <div hidden={hidden}>
        <SelectRow label={label} description={param.description} options={param.options} binding={binding} />
      </div>
    );
  }
  if (param.type === 'color') {
    return (
      <div hidden={hidden}>
        <ColorRow label={label} description={param.description} binding={binding} />
      </div>
    );
  }
  return (
    <div class={classes.row} hidden={hidden} title={param.description}>
      <span class={classes.label}>{label}</span>
      <div class={classes.control}>{control(param, binding)}</div>
    </div>
  );
}

/** A color: the swatch and hex as a button; the picker opens under the row, full width. */
function ColorRow(props: { label: string; description?: string; binding: FieldBinding }) {
  const { binding } = props;
  const open = observableState(false);
  const hex = binding.value.pipe(map((v: number) => colorToHex(v)));
  return (
    <>
      <div class={classes.row} title={props.description}>
        <span class={classes.label}>{props.label}</span>
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

function control(param: PD.Any, binding: FieldBinding): Child {
  switch (param.type) {
    case 'boolean':
      return <Toggle binding={binding} />;
    case 'number':
      return <NumberControl param={param} binding={binding} />;
    case 'text':
      return (
        <input
          class={classes.input}
          value={binding.value}
          placeholder={param.placeholder}
          onChange={(e) => binding.set(e.currentTarget.value)}
        />
      );
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
        <span class={classes.label}>{props.label}</span>
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

function GroupControl(props: { label: string; param: PD.Group<any>; binding: FieldBinding }) {
  const { param, binding } = props;
  // A flat group has no header of its own: its params are rows like the others.
  if (param.isFlat) return <Fields params={param.params} binding={binding} />;
  return (
    <Collapsible title={props.label} open={!!param.isExpanded}>
      {() => <Fields params={param.params} binding={binding} />}
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
function OnOffControl(props: { label: string; description?: string; param: PD.Mapped<any>; binding: FieldBinding }) {
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
        <span class={classes.label}>{props.label}</span>
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
function MappedControl(props: { label: string; param: PD.Mapped<any>; binding: FieldBinding }) {
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
      <SelectRow label={props.label} options={param.select.options} binding={nameBinding} />
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
function Fields(props: { params: PD.Params; binding: FieldBinding; nested?: boolean }) {
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

function defaults(param: PD.Any): any {
  if (param.type === 'group') {
    const out: Record<string, unknown> = {};
    for (const k of Object.keys(param.params)) out[k] = param.params[k].defaultValue;
    return out;
  }
  return param.defaultValue;
}

/** A titled block that opens and closes; its content is rendered while open. */
function Collapsible(props: { title: string; open?: boolean; children: () => Child }) {
  const open = observableState(!!props.open);
  return (
    <div class={classes.collapsible}>
      <button
        type="button"
        class={classes.collapsibleHeader}
        aria-expanded={open}
        onClick={() => open.set(!open.value)}
      >
        <Icon name="chevron-right" size={14} class={classes.chevron} />
        {props.title}
      </button>
      <Show when={open}>{() => <div class={classes.nested}>{props.children()}</div>}</Show>
    </div>
  );
}
