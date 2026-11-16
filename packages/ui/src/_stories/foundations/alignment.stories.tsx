/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { derive, State } from '@molstar/core/reactive';
import type { Meta, StoryObj } from '@storybook/web-components-vite';
import type { Child } from '../../dom/index.js';
import { Heading, Note, page, story } from '../support.js';

const meta: Meta = { title: 'Foundations/Alignment' };
export default meta;

/**
 * The rules (styles §5): what shares a column, and why. Each is a token, so a change moves everything that follows the
 * rule together.
 */
const Rules: readonly [rule: string, token: string, why: string][] = [
  [
    'One inset',
    '--msui-inset',
    "A header's, a row's and a form label's text start the same distance from their surface's edge: one left edge for everything in a panel.",
  ],
  [
    'Icons keep their slot',
    '--msui-icon-slot',
    "A row's icon has a fixed slot, kept when there is no icon, so the text after it starts in the same column with or without one.",
  ],
  [
    'A level is a slot',
    '--msui-indent',
    "Nesting indents by exactly an icon slot and its gap (20px): a child's icon sits under its parent's text, and a deep tree still fits a narrow panel.",
  ],
  [
    'One label column',
    '--msui-label',
    'Every form uses the same label width, so values line up across forms and across the sections of a panel.',
  ],
  [
    'Values start where controls write',
    '--msui-control-inset',
    'A plain value in a form starts where a control writes its text: a read-only row and an editable one line up.',
  ],
  [
    'A control sits in a row',
    '--msui-row',
    'A form row is a list row and a control is centered in it: a list beside a form, or a tree beside its properties, line up row by row.',
  ],
  [
    'Numbers by place',
    '.msui-num',
    'Tabular figures, right-aligned: digits stack by place, units in a column of their own.',
  ],
  [
    'On the grid',
    '--msui-unit',
    'Sizes, gaps and line heights are multiples of 4px: what the scale places lands where something else placed by it does.',
  ],
];

/** The guides: the columns, measured from the surface's edge, drawn over a demo when asked. */
const Guides: readonly [label: string, left: string, color: string][] = [
  ['inset', 'var(--msui-inset)', 'var(--msui-accent)'],
  ['text after a slot', 'calc(var(--msui-inset) + var(--msui-indent))', 'var(--msui-info)'],
  ['level 2 text', 'calc(var(--msui-inset) + 2 * var(--msui-indent))', 'var(--msui-text-subtle)'],
  ['values', 'calc(var(--msui-inset) + var(--msui-label) + var(--msui-space-3))', 'var(--msui-success)'],
  [
    'value text',
    'calc(var(--msui-inset) + var(--msui-label) + var(--msui-space-3) + var(--msui-control-inset))',
    'var(--msui-danger)',
  ],
];

function GuideLines(props: { show: State<boolean>; grid: State<boolean> }) {
  return (
    <div style="position: absolute; inset: 0; pointer-events: none;">
      <div
        style={{
          position: 'absolute',
          inset: '0',
          display: derive(props.grid, (g) => (g ? 'block' : 'none')),
          backgroundImage:
            'repeating-linear-gradient(to bottom, transparent 0 calc(var(--msui-unit) - 1px), color-mix(in oklab, var(--msui-info), transparent 80%) calc(var(--msui-unit) - 1px) var(--msui-unit))',
        }}
      />
      {Guides.map(([label, left, color]) => (
        <div
          title={label}
          style={{
            position: 'absolute',
            top: '0',
            bottom: '0',
            left,
            width: '1px',
            background: color,
            opacity: '0.8',
            display: derive(props.show, (s) => (s ? 'block' : 'none')),
          }}
        />
      ))}
    </div>
  );
}

function Icon(props: { glyph: string }) {
  return <span class="msui-slot">{props.glyph}</span>;
}

/** A panel built with the primitives: a header, a tree with and without icons, two forms, numbers. */
function AlignedPanel() {
  return (
    <div class="msui-surface" style="padding-bottom: var(--msui-space-3);">
      <div class="msui-header">
        <Icon glyph="◆" />
        Structure
      </div>
      <div class="msui-row">
        <Icon glyph="▾" />
        <span>1TQN</span>
        <span class="msui-num" style="color: var(--msui-text-subtle)">
          3,214
        </span>
      </div>
      <div class="msui-row" style="--msui-level: 1;">
        <Icon glyph="●" />
        <span>Polymer</span>
        <span class="msui-num" style="color: var(--msui-text-subtle)">
          3,012
        </span>
      </div>
      <div class="msui-row" style="--msui-level: 1;">
        <Icon glyph="" />
        <span>Ligand (no icon)</span>
        <span class="msui-num" style="color: var(--msui-text-subtle)">
          44
        </span>
      </div>
      <div class="msui-row" style="--msui-level: 1;">
        <Icon glyph="●" />
        <span>Water</span>
        <span class="msui-num" style="color: var(--msui-text-subtle)">
          158
        </span>
      </div>
      <div class="msui-header" style="height: var(--msui-row); margin-top: var(--msui-space-1);">
        Appearance
      </div>
      <div class="msui-form">
        <span class="msui-label">Representation</span>
        <span class="msui-control">Cartoon</span>
        <span class="msui-label">Color</span>
        <span class="msui-control">Chain ID</span>
        <span class="msui-label">Atoms</span>
        <span class="msui-value msui-num" style="text-align: start;">
          3,214 (read only)
        </span>
      </div>
      <div class="msui-header" style="height: var(--msui-row); margin-top: var(--msui-space-1);">
        Measurements
      </div>
      <div class="msui-form" style="grid-template-columns: var(--msui-label) minmax(0, 1fr) max-content;">
        {[
          ['Distance', '12.40', 'Å'],
          ['Angle', '104.5', '°'],
          ['RMSD', '0.83', 'Å'],
        ].map(([label, value, unit]) => [
          <span class="msui-label">{label}</span>,
          <span class="msui-value msui-num">{value}</span>,
          <span style="color: var(--msui-text-subtle)">{unit}</span>,
        ])}
      </div>
    </div>
  );
}

/** The same content laid out ad hoc: each part picked its own padding, icon width and label width. */
function AdHocPanel() {
  return (
    <div style="display: grid; align-content: start; padding: var(--msui-space-1) 0 var(--msui-space-3); border-radius: var(--msui-radius-lg); background: var(--msui-surface);">
      <div style="display: flex; align-items: center; gap: 6px; height: 40px; padding-inline: 16px; font-weight: 600;">
        <span>◆</span>Structure
      </div>
      <div style="display: flex; align-items: center; gap: 6px; height: 26px; padding-inline: 10px;">
        <span>▾</span>
        <span style="flex: 1">1TQN</span>
        <span style="color: var(--msui-text-subtle)">3,214</span>
      </div>
      {[
        ['●', 'Polymer', '3,012'],
        ['', 'Ligand (no icon)', '44'],
        ['●', 'Water', '158'],
      ].map(([glyph, label, count]) => (
        <div style="display: flex; align-items: center; gap: 6px; height: 26px; padding-inline: 10px 10px; padding-left: 26px;">
          {glyph ? <span>{glyph}</span> : null}
          <span style="flex: 1">{label}</span>
          <span style="color: var(--msui-text-subtle)">{count}</span>
        </div>
      ))}
      <div style="height: 28px; display: flex; align-items: center; margin-top: 8px; padding-inline: 14px; font-weight: 600;">
        Appearance
      </div>
      <div style="display: grid; grid-template-columns: max-content 1fr; column-gap: 10px; row-gap: 3px; align-items: center; padding-inline: 14px;">
        <span style="color: var(--msui-text-muted)">Representation</span>
        <span class="msui-control" style="height: 24px; padding-inline: 6px;">
          Cartoon
        </span>
        <span style="color: var(--msui-text-muted)">Color</span>
        <span class="msui-control" style="height: 24px; padding-inline: 6px;">
          Chain ID
        </span>
        <span style="color: var(--msui-text-muted)">Atoms</span>
        <span>3,214 (read only)</span>
      </div>
      <div style="height: 28px; display: flex; align-items: center; margin-top: 8px; padding-inline: 12px; font-weight: 600;">
        Measurements
      </div>
      <div style="display: grid; grid-template-columns: max-content max-content; column-gap: 20px; padding-inline: 12px;">
        {[
          ['Distance', '12.40 Å'],
          ['Angle', '104.5°'],
          ['RMSD', '0.83 Å'],
        ].map(([label, value]) => [<span style="color: var(--msui-text-muted)">{label}</span>, <span>{value}</span>])}
      </div>
    </div>
  );
}

function Toggles(props: { show: State<boolean>; grid: State<boolean> }) {
  return (
    <div style="display: flex; gap: var(--msui-space-5); align-items: center;">
      <label style="display: flex; gap: var(--msui-space-2); align-items: center;">
        <input type="checkbox" checked={props.show} onChange={(e) => props.show.set(e.currentTarget.checked)} />
        Column guides
      </label>
      <label style="display: flex; gap: var(--msui-space-2); align-items: center;">
        <input type="checkbox" checked={props.grid} onChange={(e) => props.grid.set(e.currentTarget.checked)} />
        4px grid
      </label>
    </div>
  );
}

function Guided(props: { show: State<boolean>; grid: State<boolean>; children: Child }) {
  return (
    <div style="position: relative;">
      {props.children}
      <GuideLines show={props.show} grid={props.grid} />
    </div>
  );
}

/**
 * Columns: a panel built with the primitives, its columns drawn over it. Header and labels on the inset; text after an
 * icon slot, with or without an icon, in one column; a child's icon under its parent's text; values and a read-only
 * value starting where a control writes; numbers by place, units in their own column.
 */
export const Columns: StoryObj = {
  render: () =>
    story(() => {
      const show = new State(true);
      const grid = new State(false);
      return (
        <div style={page}>
          <Toggles show={show} grid={grid} />
          <div style="display: grid; grid-template-columns: 360px minmax(0, 1fr); gap: var(--msui-space-5); align-items: start;">
            <Guided show={show} grid={grid}>
              <AlignedPanel />
            </Guided>
            <dl style="display: grid; gap: var(--msui-space-2); margin: 0;">
              {Guides.map(([label, , color]) => (
                <div style="display: flex; gap: var(--msui-space-2); align-items: center;">
                  <span style={{ width: '12px', height: '2px', background: color }} />
                  <span>{label}</span>
                </div>
              ))}
            </dl>
          </div>
        </div>
      );
    }),
};

/**
 * Before and after: the same content laid out ad hoc (each part with its own padding, icon width, label width and
 * row height) and with the primitives. The guides are the aligned panel's columns.
 */
export const BeforeAndAfter: StoryObj = {
  name: 'Before and after',
  render: () =>
    story(() => {
      const show = new State(true);
      const grid = new State(false);
      return (
        <div style={page}>
          <Toggles show={show} grid={grid} />
          <div style="display: grid; grid-template-columns: repeat(2, 360px); gap: var(--msui-space-5); align-items: start;">
            <div style="display: grid; gap: var(--msui-space-2);">
              <Heading>Ad hoc</Heading>
              <Guided show={show} grid={grid}>
                <AdHocPanel />
              </Guided>
            </div>
            <div style="display: grid; gap: var(--msui-space-2);">
              <Heading>Aligned</Heading>
              <Guided show={show} grid={grid}>
                <AlignedPanel />
              </Guided>
            </div>
          </div>
        </div>
      );
    }),
};

/** The rules, each with its token: what shares a column, and why. */
export const Rules_: StoryObj = {
  name: 'Rules',
  render: () =>
    story(() => (
      <div style={page}>
        <Note>
          Things that belong together visually share a column: the eye reads a panel by its edges. Each rule is a token,
          so a change to it moves everything that follows it together; ad hoc paddings are what drift.
        </Note>
        <dl style="display: grid; grid-template-columns: minmax(200px, max-content) 160px minmax(0, 1fr); gap: var(--msui-space-3) var(--msui-space-4); margin: 0; align-items: baseline;">
          {Rules.map(([rule, token, why]) => [
            <dt style="font-weight: var(--msui-weight-strong);">{rule}</dt>,
            <code style="font-size: var(--msui-text-xs);">{token}</code>,
            <dd style="margin: 0; color: var(--msui-text-muted);">{why}</dd>,
          ])}
        </dl>
      </div>
    )),
};
