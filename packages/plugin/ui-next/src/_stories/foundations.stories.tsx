/** @jsxImportSource @molstar/plugin-ui-next/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { story } from './story-root.js';

const meta: Meta = { title: 'ui-next/Foundations' };
export default meta;

const Colors = [
  'bg',
  'panel',
  'raised',
  'control',
  'hover',
  'pressed',
  'border',
  'border-strong',
  'text',
  'text-muted',
  'text-subtle',
  'accent',
  'accent-strong',
  'accent-muted',
  'danger',
  'success',
  'warning',
  'info',
];

const Accents = ['cyan', 'red', 'gray', 'green', 'purple', 'blue', 'orange'];

function Swatch(props: { token: string }) {
  return (
    <div style="display: grid; gap: 4px; font-size: var(--msp-font-size-xs)">
      <div
        style={{
          height: '40px',
          borderRadius: 'var(--msp-radius-md)',
          background: `var(${props.token})`,
          border: '1px solid var(--msp-color-border)',
        }}
      />
      <code style="color: var(--msp-color-text-muted)">{props.token}</code>
    </div>
  );
}

const grid = 'display: grid; grid-template-columns: repeat(auto-fill, minmax(140px, 1fr)); gap: 12px;';
const section = 'display: grid; gap: 12px; padding: 16px;';

/** Color tokens of the current theme (switch themes in the toolbar). */
export const ColorTokens: StoryObj = {
  render: () =>
    story(() => (
      <div style={section}>
        <h2 style="font-size: var(--msp-font-size-lg)">Colors</h2>
        <div style={grid}>
          {Colors.map((c) => (
            <Swatch token={`--msp-color-${c}`} />
          ))}
        </div>
        <h2 style="font-size: var(--msp-font-size-lg)">Accents</h2>
        <div style={grid}>
          {Accents.map((c) => (
            <Swatch token={`--msp-accent-${c}`} />
          ))}
        </div>
      </div>
    )),
};

/** The same subtree in every theme: a `data-msp-theme` attribute recomputes all derived tokens below it. */
export const Themes: StoryObj = {
  render: () =>
    story(() => (
      <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(220px, 1fr));">
        {['light', 'dark', 'blue'].map((theme) => (
          <div
            data-msp-theme={theme}
            style="background: var(--msp-color-bg); color: var(--msp-color-text); padding: 16px; display: grid; gap: 8px"
          >
            <b>{theme}</b>
            <div style="background: var(--msp-color-panel); padding: 8px; border-radius: var(--msp-radius-md)">
              Panel <span style="color: var(--msp-color-text-muted)">muted</span>{' '}
              <span style="color: var(--msp-color-accent)">accent</span>
            </div>
            <div style="border: 1px solid var(--msp-color-border); padding: 8px; border-radius: var(--msp-radius-md)">
              Bordered
            </div>
          </div>
        ))}
      </div>
    )),
};

/** Type scale and spacing. */
export const TypeAndSpace: StoryObj = {
  render: () =>
    story(() => (
      <div style={section}>
        {['xs', 'sm', 'md', 'lg'].map((s) => (
          <div style={{ fontSize: `var(--msp-font-size-${s})` }}>
            --msp-font-size-{s}: The quick brown fox jumps over the lazy dog
          </div>
        ))}
        <code style="font-family: var(--msp-font-mono)">--msp-font-mono: MALYGTHSHGLFKKLGIPGPTPLP</code>
        <div style="display: flex; gap: 8px; align-items: end">
          {[1, 2, 3, 4, 5, 6].map((n) => (
            <div style="display: grid; gap: 4px; justify-items: center; font-size: var(--msp-font-size-xs)">
              <div
                style={{
                  width: `var(--msp-space-${n})`,
                  height: `var(--msp-space-${n})`,
                  background: 'var(--msp-color-accent)',
                }}
              />
              {n}
            </div>
          ))}
        </div>
      </div>
    )),
};
