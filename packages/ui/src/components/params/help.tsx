/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { derive, type ReadonlyState, State } from '@molstar/core/reactive';
import type { Legend as MolLegend } from '@molstar/core/util/legend';
import type { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { type Child, Show, setAttr, watch } from '../../dom/index.js';
import { state } from '../../state.js';
import { colorToHex } from '../color/color-picker.js';
import { Icon } from '../icon/icon.js';
import { Legend, type LegendItem, type LegendRamp } from '../legend/legend.js';
import { classes } from './_generated/help.module.css.js';

/** What a param says about itself: a description and a legend (a color theme's table or scale). */
export interface ParamHelp {
  description?: string;
  legend?: MolLegend;
}

/**
 * A param's help for its current value: its `help(value)` when it has one (a color theme: the chosen theme's
 * description and legend), else its own description and legend; undefined when it has none.
 */
export function helpOf(param: PD.Any, value: ReadonlyState<any>): ReadonlyState<ParamHelp | undefined> | undefined {
  const has = (h: ParamHelp | undefined) => (h && (h.description || h.legend) ? h : undefined);
  if (param.help) {
    const help = param.help;
    return derive(
      value,
      (v) => {
        try {
          return has(help(v));
        } catch {
          return undefined;
        }
      },
      (a, b) => a?.description === b?.description && a?.legend === b?.legend,
    );
  }
  const own = has({ description: param.description, legend: param.legend });
  return own ? new State<ParamHelp | undefined>(own) : undefined;
}

/**
 * A row's label with its help: the label's text (the UI's own: it wraps, never cut, styles §2), and when there is help a
 * small ? after it (shown on the row's hover, or while open); the text opens the help too. The help itself goes under
 * the row (`block`), its text on the label column. `id` goes on the label's text (not the ?): the row's control is
 * labelled by it (`aria-labelledby`). `about` names the ? toggle by the label (the form's labels).
 */
export function withHelp(
  label: string,
  help: ReadonlyState<ParamHelp | undefined> | undefined,
  id: string,
  about: (label: string) => string,
): { head: Child; block: Child } {
  if (!help) {
    const head = (
      <span class={classes.text} id={id} data-tip-overflow="">
        {label}
      </span>
    );
    return { head, block: null };
  }
  const open = state(false);
  const has = derive(help, (h) => !!h);
  const text = (
    <span
      class={classes.text}
      id={id}
      data-tip-overflow=""
      onClick={(e) => {
        if (!has.value) return;
        e.stopPropagation();
        open.set(!open.value);
      }}
    >
      {label}
    </span>
  ) as HTMLElement;
  const head = (
    <>
      {text}
      <Show when={has}>
        {() => {
          const toggle = (
            <button
              type="button"
              class={classes.toggle}
              aria-label={about(label)}
              data-tip={about(label)}
              onClick={(e) => {
                e.stopPropagation();
                open.set(!open.value);
              }}
            >
              <Icon name="circle-help" />
            </button>
          ) as HTMLElement;
          watch([open], () => setAttr(toggle, 'aria-expanded', String(open.value)));
          return toggle;
        }}
      </Show>
    </>
  );
  watch([has], () => setAttr(text, 'data-help', has.value));
  return { head, block: <HelpBlock help={help} open={open} /> };
}

function HelpBlock(props: { help: ReadonlyState<ParamHelp | undefined>; open: State<boolean> }) {
  return (
    <Show when={props.open} reveal>
      {() => (
        <Show when={props.help}>
          {(h) => (
            <div class={classes.help}>
              {h.description ? <p class={classes.description}>{described(h.description)}</p> : null}
              {h.legend ? <LegendView legend={h.legend} /> : null}
            </div>
          )}
        </Show>
      )}
    </Show>
  );
}

/** A description's `code` as code; the rest as text. */
function described(text: string): Child[] {
  return text
    .split(/(`[^`]+`)/)
    .map((part) =>
      part.startsWith('`') && part.endsWith('`') && part.length > 1 ? <code>{part.slice(1, -1)}</code> : part,
    );
}

/** A Mol* legend as the UI's `Legend`: a table's entries as swatches, a scale as a ramp with its ends. */
export function LegendView(props: { legend: MolLegend }) {
  const { legend } = props;
  if (legend.kind === 'table-legend') {
    const items: LegendItem[] = legend.table.map(([label, color]) => ({ label, color: colorToHex(color) }));
    return <Legend class={classes.table} items={items} />;
  }
  return <Legend ramps={[legendRamp(legend.minLabel, legend.maxLabel, legend.colors)]} />;
}

/** A color list (plain colors spread evenly, or `[color, position]` pairs) as a ramp's stops. */
export function legendRamp(lo: string, hi: string, colors: readonly unknown[]): LegendRamp {
  const n = colors.length;
  return {
    lo,
    hi,
    stops: colors.map((c, i) =>
      Array.isArray(c)
        ? { t: c[1] as number, color: colorToHex(c[0] as number) }
        : { t: n === 1 ? 0 : i / (n - 1), color: colorToHex(c as number) },
    ),
  };
}
