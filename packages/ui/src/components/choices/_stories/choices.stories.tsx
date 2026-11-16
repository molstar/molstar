/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { type Child, onCleanup } from '../../../dom/index.js';
import { story } from '../../../_stories/support.js';
import { state } from '../../../state.js';
import { Select } from '../../select/select.js';
import { TextInput } from '../../text-input/text-input.js';
import { installTooltips } from '../../tooltip/tooltip.js';
import { type Choice, Choices } from '../choices.js';

const meta: Meta = { title: 'Components/Choices' };
export default meta;

/** A surface with tooltips: what each story is in. */
function Panel(props: { width?: number; children: Child }) {
  return (
    <div
      ref={(el: HTMLElement) => onCleanup(installTooltips(el))}
      class="msui-surface"
      style={`width: ${props.width ?? 360}px; margin: var(--msui-space-4);`}
    >
      {props.children}
    </div>
  );
}

const quickStyles: Choice[] = [
  { id: 'default', label: 'Default' },
  { id: 'cartoon', label: 'Cartoon' },
  { id: 'spacefill', label: 'Spacefill' },
  { id: 'surface', label: 'Surface' },
];

/**
 * In a form, with an active choice: a single-select group (a radio group: Tab reaches the checked choice, the arrows,
 * Home and End move and pick). Its label is in the form's label column and its choices start where the fields do and
 * end where they end; each line of choices is a form row.
 */
export const InAForm: StoryObj = {
  render: () =>
    story(() => {
      const style = state<string | undefined>('cartoon');
      const name = state('1tqn');
      const color = state('chain-id');
      return (
        <Panel width={420}>
          <div class="msui-form">
            <span class="msui-label">Name</span>
            <TextInput size="sm" aria-label="Name" value={name} onInput={(v) => name.set(v)} />
            <Choices label="Quick style" choices={quickStyles} active={style} onPick={(id) => style.set(id)} />
            <span class="msui-label">Color</span>
            <Select
              size="sm"
              aria-label="Color"
              options={[
                ['chain-id', 'Chain Id'],
                ['uniform', 'Uniform'],
              ]}
              value={color}
              onChange={(v) => color.set(v)}
            />
            <span class="msui-label">Active</span>
            <span class="msui-value">{style}</span>
          </div>
        </Panel>
      );
    }),
};

/** Without an active choice: a row of one-click actions (plain buttons, each a Tab stop); the pick echoed below. */
export const Actions: StoryObj = {
  render: () =>
    story(() => {
      const picked = state('none');
      return (
        <Panel>
          <div class="msui-form">
            <Choices
              label="Wiggle"
              choices={[
                { id: 'dynamics', label: 'Dynamics' },
                { id: 'uncertainty', label: 'Uncertainty' },
                { id: 'clear', label: 'Clear', title: 'Stop the wiggle' },
              ]}
              onPick={(id) => picked.set(id)}
            />
            <span class="msui-label">Picked</span>
            <span class="msui-value">{picked}</span>
          </div>
        </Panel>
      );
    }),
};

/** Icons with labels, and icon-only choices: each named by its label, which is also its tooltip. */
export const Icons: StoryObj = {
  render: () =>
    story(() => {
      const view = state<string | undefined>('front');
      const visibility = state<string | undefined>('show');
      return (
        <Panel>
          <div class="msui-form">
            <Choices
              label="Visibility"
              choices={[
                { id: 'show', label: 'Show', icon: 'eye' },
                { id: 'hide', label: 'Hide', icon: 'eye-off' },
              ]}
              active={visibility}
              onPick={(id) => visibility.set(id)}
            />
            <Choices
              label="View"
              iconOnly
              choices={[
                { id: 'front', label: 'Front', icon: 'scan' },
                { id: 'side', label: 'Side', icon: 'flip-vertical-2' },
                { id: 'top', label: 'Top', icon: 'layers' },
                { id: 'reset', label: 'Reset', icon: 'rotate-ccw', title: 'Reset the camera' },
              ]}
              active={view}
              onPick={(id) => view.set(id)}
            />
            <span class="msui-label">View</span>
            <span class="msui-value">{view}</span>
          </div>
        </Panel>
      );
    }),
};

/** A disabled choice: greyed, its reason in its tooltip, not picked, skipped by the arrows. */
export const Disabled: StoryObj = {
  render: () =>
    story(() => {
      const style = state<string | undefined>('cartoon');
      return (
        <Panel>
          <div class="msui-form">
            <Choices
              label="Quick style"
              choices={[
                { id: 'default', label: 'Default' },
                { id: 'cartoon', label: 'Cartoon' },
                { id: 'spacefill', label: 'Spacefill', disabled: 'Too many atoms for spacefill' },
                { id: 'surface', label: 'Surface' },
              ]}
              active={style}
              onPick={(id) => style.set(id)}
            />
          </div>
        </Panel>
      );
    }),
};

const many: Choice[] = [
  'Cartoon',
  'Backbone',
  'Ball & Stick',
  'Spacefill',
  'Surface',
  'Gaussian',
  'Putty',
  'Label',
  'Ellipsoid',
].map((label) => ({ id: label.toLowerCase(), label }));

/**
 * Many choices wrap into rows of equal columns (as many as fit, or `columns`); each line is a form row, so the rows
 * after them stay on the row grid.
 */
export const Many: StoryObj = {
  render: () =>
    story(() => {
      const a = state<string | undefined>('cartoon');
      const b = state<string | undefined>('surface');
      return (
        <Panel>
          <div class="msui-form">
            <Choices label="As many as fit" choices={many} active={a} onPick={(id) => a.set(id)} />
            <Choices label="Three columns" columns={3} choices={many} active={b} onPick={(id) => b.set(id)} />
            <span class="msui-label">After</span>
            <span class="msui-value">on the row grid</span>
          </div>
        </Panel>
      );
    }),
};

/** Narrow: the choices wrap into fewer columns, each still as wide as the widest label (never cut). */
export const Narrow: StoryObj = {
  render: () =>
    story(() => {
      const active = state<string | undefined>('default');
      return (
        <Panel width={240}>
          <div class="msui-form">
            <Choices
              label="Representation"
              choices={[
                { id: 'default', label: 'Default' },
                { id: 'cartoon', label: 'Cartoon' },
                { id: 'spacefill', label: 'Spacefill' },
                { id: 'surface', label: 'Molecular surface' },
              ]}
              active={active}
              onPick={(id) => active.set(id)}
            />
          </div>
        </Panel>
      );
    }),
};

/** On its own (no form, no label): the choices take the width they are given. */
export const Standalone: StoryObj = {
  render: () =>
    story(() => (
      <Panel>
        <div style="padding-inline: var(--msui-inner);">
          <Choices aria-label="Quick style" choices={quickStyles} onPick={() => {}} />
        </div>
      </Panel>
    )),
};
