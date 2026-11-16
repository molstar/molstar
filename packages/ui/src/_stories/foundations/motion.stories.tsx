/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { derive, State } from '@molstar/core/reactive';
import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { Button } from '../../components/button/button.js';
import { Icon } from '../../components/icon/icon.js';
import { Select } from '../../components/select/select.js';
import { Heading, Note, story } from '../support.js';

const meta: Meta = { title: 'Foundations/Motion' };
export default meta;

const Principles: readonly [string, string][] = [
  ['Move what the eye follows', 'A fold opening, a chevron turning, a state changing color; never decoration.'],
  [
    'Pick the move by what must not change',
    'Grow to its height for content in a flow; slide for panels; fade (opacity only) for what keeps its size.',
  ],
  [
    'Cheap to run',
    'transform, opacity and grid tracks; no layout read in a frame; nothing on the first render or for the rows of a list.',
  ],
  ['Quiet when asked', 'Every duration is 0 under reduced motion (the tokens do it).'],
];

const Tokens: readonly [string, string][] = [
  ['--msui-motion-fast', '100 ms: hover, press, a chevron'],
  ['--msui-motion-base', '200 ms: content that opens or closes'],
  ['--msui-motion-slow', '360 ms: what only fades'],
  ['--msui-ease', 'Out: state changes'],
  ['--msui-ease-in-out', 'In and out: opening and closing, so growth is even'],
];

/** One fold: the chevron turns (fast), the body grows to its height and back (base), by its grid track. */
function Fold(props: { label: string; items: readonly string[] }) {
  const open = new State(false);
  return (
    <div>
      <button
        type="button"
        class="msui-row"
        style="width: 100%; border: 0; background: none; color: inherit; font: inherit; text-align: start; cursor: pointer;"
        aria-expanded={open}
        onClick={() => open.set(!open.value)}
      >
        <span
          class="msui-slot"
          style={{
            '--msui-icon-size': 'var(--msui-icon-sm)',
            transition: 'transform var(--msui-motion-fast) var(--msui-ease)',
            transform: derive(open, (o) => (o ? 'rotate(90deg)' : 'none')),
          }}
        >
          <Icon name="chevron-right" />
        </span>
        <span>{props.label}</span>
        <span />
      </button>
      <div
        style={{
          display: 'grid',
          transition: 'grid-template-rows var(--msui-motion-base) var(--msui-ease-in-out)',
          gridTemplateRows: derive(open, (o) => (o ? '1fr' : '0fr')),
        }}
      >
        <div style="overflow: hidden;">
          {props.items.map((i) => (
            <div class="msui-row" style="--msui-level: 1;">
              <span class="msui-slot" />
              <span style="color: var(--msui-text-muted)">{i}</span>
              <span />
            </div>
          ))}
        </div>
      </div>
    </div>
  );
}

/** Motion: the principles, the tokens, and the moves (styles §7). Slow motion multiplies every duration by five. */
export const Moves: StoryObj = {
  render: () =>
    story(() => {
      const slow = new State(false);
      const scale = (ms: number) => derive(slow, (s) => (s ? `${ms * 5}ms` : ''));
      return (
        <div
          style={{
            display: 'grid',
            gap: 'var(--msui-space-4)',
            padding: 'var(--msui-space-5)',
            maxWidth: '960px',
            '--msui-motion-fast': scale(100),
            '--msui-motion-base': scale(200),
            '--msui-motion-slow': scale(360),
          }}
        >
          <Heading>Principles</Heading>
          <dl style="display: grid; grid-template-columns: minmax(200px, max-content) 1fr; gap: var(--msui-space-2) var(--msui-space-4); margin: 0;">
            {Principles.map(([rule, why]) => [
              <dt style="font-weight: var(--msui-weight-strong);">{rule}</dt>,
              <dd style="margin: 0; color: var(--msui-text-muted);">{why}</dd>,
            ])}
          </dl>
          <Heading>Tokens</Heading>
          <dl style="display: grid; grid-template-columns: 200px 1fr; gap: var(--msui-space-1) var(--msui-space-4); margin: 0;">
            {Tokens.map(([t, use]) => [
              <dt>
                <code style="font-size: var(--msui-text-xs);">{t}</code>
              </dt>,
              <dd style="margin: 0; color: var(--msui-text-muted);">{use}</dd>,
            ])}
          </dl>
          <label style="display: flex; gap: var(--msui-space-2); align-items: center;">
            <input type="checkbox" checked={slow} onChange={(e) => slow.set(e.currentTarget.checked)} />
            Slow motion (×5)
          </label>
          <div style="display: grid; grid-template-columns: repeat(auto-fill, minmax(280px, 1fr)); gap: var(--msui-space-4);">
            <section class="msui-surface">
              <div class="msui-header">Fold</div>
              <Fold label="1TQN" items={['Polymer', 'Ligand', 'Water']} />
              <Fold label="Advanced" items={['Quality', 'Smoothness']} />
            </section>
            <section class="msui-surface" style="padding-bottom: var(--msui-space-3);">
              <div class="msui-header">State</div>
              <div style="padding: 0 var(--msui-inner) var(--msui-space-2);">
                <Note>Hover and press change color (fast); an option list opens over its field.</Note>
              </div>
              <div class="msui-form">
                <span class="msui-label">Buttons</span>
                <div style="display: flex; gap: var(--msui-space-2);">
                  <Button size="sm">Hover me</Button>
                  <Button size="sm" variant="primary">
                    Press me
                  </Button>
                </div>
                <span class="msui-label">Select</span>
                <Select
                  size="sm"
                  aria-label="Color"
                  value="chain-id"
                  options={[
                    ['chain-id', 'Chain ID'],
                    ['element', 'Element'],
                    ['uniform', 'Uniform'],
                  ]}
                />
              </div>
            </section>
          </div>
        </div>
      );
    }),
};
