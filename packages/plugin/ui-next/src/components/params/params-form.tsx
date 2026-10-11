/** @jsxImportSource @molstar/plugin-ui-next/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { camelCaseToWords } from '@molstar/core/util/string';
import { distinctUntilChanged, map, type Observable } from 'rxjs';
import { Show } from '../../dom/control-flow.js';
import type { Child } from '../../dom/jsx-types.js';
import type { ParamsFormModel } from '../../forms/params-form.js';
import { observableState } from '../../reactive/observable-state.js';
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
        <MappedControl label={label} param={param} binding={binding} />
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

function control(param: PD.Any, binding: FieldBinding): Child {
  switch (param.type) {
    case 'boolean':
      return <Toggle binding={binding} />;
    case 'number':
      return <NumberControl param={param} binding={binding} />;
    case 'select':
      return <SelectControl options={param.options} binding={binding} />;
    case 'text':
      return (
        <input
          class={classes.input}
          value={binding.value}
          placeholder={param.placeholder}
          onChange={(e) => binding.set(e.currentTarget.value)}
        />
      );
    case 'color':
      return <ColorControl binding={binding} />;
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
        type="range"
        min={String(param.min)}
        max={String(param.max)}
        step={String(param.step ?? (param.max - param.min) / 100)}
        value={binding.value.pipe(map(String))}
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

function SelectControl(props: { options: readonly PD.SelectOption<any>[]; binding: FieldBinding }) {
  const { options, binding } = props;
  // Values may be any type: the select works with indices.
  const select = (
    <select
      class={classes.select}
      onChange={(e) => binding.set(options[Number(e.currentTarget.value)][0])}
      value={binding.value.pipe(map((v) => String(options.findIndex((o) => o[0] === v))))}
    >
      {options.map(([, label], i) => (
        <option value={String(i)}>{label}</option>
      ))}
    </select>
  ) as HTMLSelectElement;
  return select;
}

function ColorControl(props: { binding: FieldBinding }) {
  const { binding } = props;
  const hex = binding.value.pipe(map((v: number) => `#${(v >>> 0).toString(16).padStart(6, '0').slice(-6)}`));
  return (
    <label class={classes.color}>
      <span class={classes.swatch} style={{ background: hex }} />
      <span class={classes.hex}>{hex}</span>
      <input
        type="color"
        value={hex}
        onInput={(e) => binding.set(Number.parseInt(e.currentTarget.value.slice(1), 16))}
      />
    </label>
  );
}

function GroupControl(props: { label: string; param: PD.Group<any>; binding: FieldBinding }) {
  const { param, binding } = props;
  return (
    <Collapsible title={props.label} open={!!param.isExpanded}>
      {() => <Fields params={param.params} binding={binding} />}
    </Collapsible>
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
      <div class={classes.row}>
        <span class={classes.label}>{props.label}</span>
        <div class={classes.control}>
          <SelectControl options={param.select.options} binding={nameBinding} />
        </div>
      </div>
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
          if (inner.type === 'group') return <Fields params={inner.params} binding={paramsBinding} nested />;
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
    <div class={props.nested ? classes.nested : undefined}>
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
