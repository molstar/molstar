/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The params the sentence does not hold (commands §5.3): a fold, closed by default, its summary what differs from the
// defaults; open, a params form whose edits go to the invocation (`set`), so the run, the preview and `check` see them.

import { derive } from '@molstar/core/reactive';
import { Color } from '@molstar/core/util/color';
import type { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import type { CommandInvocation } from '../../commands/invocation.js';
import type { ChangedOption } from '../../commands/sentence.js';
import { effect, Show, setAttr, setText, uniqueId, watch } from '../../dom/index.js';
import { ParamsFormModel } from '../../forms/params-form.js';
import { state } from '../../state.js';
import { Icon } from '../icon/icon.js';
import { ParamsForm } from '../params/params-form.js';
import { classes } from './_generated/options-fold.module.css.js';

export function OptionsFold(props: { invocation: CommandInvocation<any> }) {
  const { invocation } = props;
  // The form's params change only when the command reshapes (not at each edit): the form is made again then.
  const shape = derive(
    invocation.state,
    (s) => ({ params: s.params, formOnly: s.formOnly }),
    (a, b) => a.params === b.params && a.formOnly === b.formOnly,
  );
  return (
    <Show when={derive(shape, (x) => (x.formOnly.length ? x : undefined))}>
      {(x) => {
        const open = state(false);
        const bodyId = uniqueId('msui-command-options');
        const summary = document.createTextNode('');
        const head = (
          <button type="button" class={classes.head} aria-controls={bodyId} onClick={() => open.set(!open.value)}>
            <span class={classes.chevron}>
              <Icon name="chevron-right" />
            </span>
            <span class={classes.title}>Options</span>
            <span class={classes.summary}>{summary}</span>
          </button>
        ) as HTMLElement;
        watch([open, invocation.options], () => {
          setAttr(head, 'aria-expanded', String(open.value));
          setText(summary, summarize(invocation.options.value.changed, x.params, x.formOnly.length));
        });
        return (
          <div class={classes.fold}>
            {head}
            <Show when={open}>
              {() => {
                const own: PD.Params = {};
                for (const name of x.formOnly) own[name] = x.params[name];
                const form = new ParamsFormModel(own, invocation.options.value.values as PD.Values<PD.Params>, {
                  onChange: (values, name) => invocation.set(name, values[name]),
                });
                // Set from elsewhere (a command's own form, the assistant): the form shows it.
                effect(
                  derive(invocation.options, (o) => o.values),
                  (v) => {
                    if (v !== form.values.value) form.load(v as PD.Values<PD.Params>);
                  },
                );
                return (
                  <div id={bodyId} class={classes.body}>
                    <ParamsForm form={form} />
                  </div>
                );
              }}
            </Show>
          </div>
        );
      }}
    </Show>
  );
}

/** "3 · size 0.4, color #ff0000": how many, and what differs from the defaults (or that none does). */
export function summarize(changed: readonly ChangedOption[], params: PD.Params, count: number): string {
  if (!changed.length) return `${count} at their defaults`;
  return changed.map((c) => `${c.label.toLowerCase()} ${formatValue(params[c.name], c.value)}`).join(', ');
}

function formatValue(param: PD.Any | undefined, value: unknown): string {
  if (param?.type === 'color' && typeof value === 'number') return Color.toHexStyle(value as Color);
  if (param?.type === 'select') return String(param.options.find((o) => o[0] === value)?.[1] ?? value);
  if (typeof value === 'boolean') return value ? 'on' : 'off';
  if (typeof value === 'number') return String(Math.round(value * 1000) / 1000);
  if (typeof value === 'string') return value.length > 24 ? `${value.slice(0, 24)}…` : value;
  return 'changed';
}
