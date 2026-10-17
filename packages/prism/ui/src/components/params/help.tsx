/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Legend } from '@molstar/core/util/legend';
import type { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { Show } from '../../dom/control-flow.js';
import type { Child } from '../../dom/jsx-types.js';
import { type ObservableState, observableState } from '../../reactive/observable-state.js';
import type { ReadonlyState } from '../../reactive/readonly-state.js';
import { derive } from '../../reactive/state.js';
import { colorToHex } from '../color/color-picker.js';
import { Icon } from '../icon/icon.js';
import { colorListGradient } from './extra-controls.js';
import { classes } from './_generated/help.module.css.js';

/** What a param says about itself: a description and a legend (a color theme's table or scale). */
export interface ParamHelp {
  description?: string;
  legend?: Legend;
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
  return own ? observableState<ParamHelp | undefined>(own).asReadonly : undefined;
}

/**
 * A row's label with its help: the label text, and a small ? when there is help (shown on hover, or while open); the
 * help itself opens under the row (`HelpBlock`). Returns the label's content and the block. `id` goes on the label's
 * text (not the ?), for the row's control to be labelled by (`aria-labelledby`).
 */
export function withHelp(label: string, help: ReadonlyState<ParamHelp | undefined> | undefined, id?: string) {
  if (!help)
    return {
      head: (
        <span class={classes.text} id={id} data-tip-overflow="">
          {label}
        </span>
      ) as Child,
      block: null as Child,
    };
  const open = observableState(false);
  const has = derive(help, (h) => !!h);
  const head = (
    <>
      {/* With help, the text is underlined (dotted) and opens it too, as the ? does. */}
      <span
        class={classes.text}
        id={id}
        data-tip-overflow=""
        data-help={derive(has, (h) => (h ? '' : undefined))}
        onClick={(e) => {
          if (!e.currentTarget.hasAttribute('data-help')) return;
          e.stopPropagation();
          open.set(!open.value);
        }}
      >
        {label}
      </span>
      <Show when={has}>
        {() => (
          <button
            type="button"
            class={classes.toggle}
            aria-label={`About ${label}`}
            title={`About ${label}`}
            aria-expanded={open}
            onClick={(e) => {
              e.stopPropagation();
              open.set(!open.value);
            }}
          >
            <Icon name="circle-help" size={12} />
          </button>
        )}
      </Show>
    </>
  );
  return { head, block: <HelpBlock help={help} open={open} /> };
}

function HelpBlock(props: { help: ReadonlyState<ParamHelp | undefined>; open: ObservableState<boolean> }) {
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

/** A legend: a table (a swatch and a name per entry) or a scale (the colors as a bar, its ends labelled). */
export function LegendView(props: { legend: Legend }) {
  const { legend } = props;
  if (legend.kind === 'table-legend') {
    return (
      <ul class={classes.table}>
        {legend.table.map(([name, color]) => (
          <li class={classes.entry}>
            <span class={classes.swatch} style={{ background: colorToHex(color) }} />
            <span class={classes.name}>{name}</span>
          </li>
        ))}
      </ul>
    );
  }
  return (
    <div class={classes.scale}>
      <span
        class={classes.bar}
        style={{ background: colorListGradient({ kind: 'interpolate', colors: legend.colors }) }}
      />
      <span class={classes.ends}>
        <span>{legend.minLabel}</span>
        <span>{legend.maxLabel}</span>
      </span>
    </div>
  );
}
