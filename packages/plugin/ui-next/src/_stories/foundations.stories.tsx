/** @jsxImportSource @molstar/plugin-ui-next/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { story } from './story-root.js';

const meta: Meta = { title: 'Overview/Foundations' };
export default meta;

/** The tokens by role (`styles/theme.css`): every shade derives from a theme's base color and text color. */
const Groups: [title: string, note: string, tokens: string[]][] = [
  [
    'Surfaces',
    'From the base to the strongest: a panel, a control, hover, pressed; a well is a block inside a panel (an open action, a section), a step the other way.',
    ['bg', 'panel', 'raised', 'control', 'hover', 'pressed', 'well'],
  ],
  ['Text', 'Body, labels and secondary text, hints and placeholders.', ['text', 'text-muted', 'text-subtle']],
  [
    'States',
    'The accent only for actions and the picked item; the rest for messages.',
    ['accent', 'accent-strong', 'danger', 'success', 'warning', 'info', 'focus-ring'],
  ],
  [
    'Brand',
    "The logo's colors (the wordmark's flash): lighter on dark themes.",
    ['brand-1', 'brand-2', 'brand-3', 'brand-4'],
  ],
];

const Glass = ['--msp-glass-bg', '--msp-glass-panel-bg'];
const Accents = ['cyan', 'red', 'gray', 'green', 'purple', 'blue', 'orange'];

function Swatch(props: { token: string }) {
  return (
    <div style="display: grid; gap: 4px; font-size: var(--msp-font-size-xs)">
      <div
        style={{
          height: '36px',
          borderRadius: 'var(--msp-radius-md)',
          background: `var(${props.token})`,
          boxShadow: 'inset 0 0 0 1px var(--msp-glass-edge)',
        }}
      />
      <code style="color: var(--msp-color-text-muted)">{props.token}</code>
    </div>
  );
}

const grid = 'display: grid; grid-template-columns: repeat(auto-fill, minmax(132px, 1fr)); gap: 12px;';
const section = 'display: grid; gap: 12px; padding: 16px; max-width: 960px;';
const heading = 'margin: 8px 0 0; font-size: var(--msp-font-size-md);';
const note = 'margin: 0; color: var(--msp-color-text-muted); font-size: var(--msp-font-size-sm);';

/** The color tokens of the current theme, by role (switch themes in the toolbar). */
export const ColorTokens: StoryObj = {
  render: () =>
    story(() => (
      <div style={section}>
        {Groups.map(([title, text, tokens]) => (
          <>
            <h2 style={heading}>{title}</h2>
            <p style={note}>{text}</p>
            <div style={grid}>
              {tokens.map((t) => (
                <Swatch token={`--msp-color-${t}`} />
              ))}
            </div>
          </>
        ))}
        <h2 style={heading}>Glass</h2>
        <p style={note}>Over the canvas: translucent, blurred; the panel's is more opaque, for reading.</p>
        <div style={grid}>
          {Glass.map((t) => (
            <Swatch token={t} />
          ))}
        </div>
        <h2 style={heading}>Node kinds</h2>
        <p style={note}>For extensions' and node kinds' brands (the classic card accents).</p>
        <div style={grid}>
          {Accents.map((c) => (
            <Swatch token={`--msp-accent-${c}`} />
          ))}
        </div>
      </div>
    )),
};

/**
 * How the surfaces layer: a panel on the base, a well inside it (an open action), its fields at the control tone,
 * hover and pressed above; text in its three tones; the accent only on the action.
 */
export const Surfaces: StoryObj = {
  render: () =>
    story(() => (
      <div style="padding: 16px; display: grid; gap: 8px; max-width: 360px; font-size: var(--msp-font-size-sm)">
        <div style="padding: 12px; border-radius: var(--msp-radius-lg); background: var(--msp-color-panel); display: grid; gap: 8px">
          <b>Panel</b>
          <span style="color: var(--msp-color-text-muted)">A label, muted · a hint, subtle</span>
          <div style="padding: 8px; border-radius: var(--msp-radius-md); background: var(--msp-color-well); display: grid; gap: 4px">
            <b>Well (an open action)</b>
            <div style="height: 24px; border-radius: var(--msp-radius-md); background: var(--msp-color-control); padding: 0 8px; display: flex; align-items: center">
              Control
            </div>
            <div style="height: 24px; border-radius: var(--msp-radius-md); background: var(--msp-color-hover); padding: 0 8px; display: flex; align-items: center">
              Hover
            </div>
            <div style="height: 24px; border-radius: var(--msp-radius-md); background: var(--msp-color-pressed); padding: 0 8px; display: flex; align-items: center">
              Pressed
            </div>
            <span style="justify-self: end; color: var(--msp-color-accent); font-weight: 600">Download</span>
          </div>
          <span style="color: var(--msp-color-text-subtle)">Subtle: placeholders, counts, dates</span>
        </div>
      </div>
    )),
};

/** The same panel in every theme: `data-msp-theme` on any element recomputes every derived token below it. */
export const Themes: StoryObj = {
  render: () =>
    story(() => (
      <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(220px, 1fr));">
        {['light', 'dark', 'blue'].map((theme) => (
          <div
            data-msp-theme={theme}
            style="background: var(--msp-color-bg); color: var(--msp-color-text); padding: 16px; display: grid; gap: 8px; font-size: var(--msp-font-size-sm)"
          >
            <b style="text-transform: capitalize">{theme}</b>
            <div style="background: var(--msp-color-panel); padding: 8px; border-radius: var(--msp-radius-md); display: grid; gap: 6px">
              <span>
                Panel <span style="color: var(--msp-color-text-muted)">muted</span>{' '}
                <span style="color: var(--msp-color-text-subtle)">subtle</span>
              </span>
              <div style="background: var(--msp-color-well); padding: 6px 8px; border-radius: var(--msp-radius-md); display: flex; justify-content: space-between">
                <span style="background: var(--msp-color-control); padding: 0 8px; border-radius: var(--msp-radius-md)">
                  Field
                </span>
                <span style="color: var(--msp-color-accent); font-weight: 600">Apply</span>
              </div>
            </div>
          </div>
        ))}
      </div>
    )),
};

/** Type scale (12 px in panels; 14 px for what is typed and for panel titles) and the spacing steps. */
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
                  background: 'var(--msp-color-text-muted)',
                }}
              />
              {n}
            </div>
          ))}
        </div>
      </div>
    )),
};
