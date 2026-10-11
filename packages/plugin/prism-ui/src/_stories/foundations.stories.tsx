/** @jsxImportSource @molstar/plugin-prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { map } from 'rxjs';
import { Button } from '../components/button/button.js';
import { Icon } from '../components/icon/icon.js';
import { Select } from '../components/select/select.js';
import { Show } from '../dom/control-flow.js';
import { observableState } from '../reactive/observable-state.js';
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

/** One row that folds open: its chevron turns, its children grow in (Show's `reveal`) and collapse back. */
function FoldRow(props: { label: string; children: readonly string[] }) {
  const open = observableState(false);
  return (
    <div>
      <button
        type="button"
        style="display: flex; align-items: center; gap: 6px; width: 100%; height: 26px; padding: 0 4px; border-radius: var(--msp-radius-md); color: var(--msp-color-text); font-size: var(--msp-font-size-sm);"
        aria-expanded={open.asReadonly}
        onClick={() => open.set(!open.value)}
      >
        <span
          style={{
            display: 'inline-flex',
            transition: 'transform var(--msp-motion-fast) var(--msp-ease)',
            transform: open.asReadonly.pipe(map((o) => (o ? 'rotate(90deg)' : 'none'))),
          }}
        >
          <Icon name="chevron-right" size={14} />
        </span>
        {props.label}
      </button>
      <Show when={open.asReadonly} reveal>
        {() => (
          <div style="display: grid; padding-left: 24px;">
            {props.children.map((c) => (
              <div style="height: 26px; display: flex; align-items: center; color: var(--msp-color-text-muted); font-size: var(--msp-font-size-sm);">
                {c}
              </div>
            ))}
          </div>
        )}
      </Show>
    </div>
  );
}

/**
 * Motion: what moves, and how (`--msp-motion-fast` for state changes, `--msp-motion-base` for content that opens,
 * `--msp-ease`, and `--msp-ease-reveal` for opening and closing). Content that opens grows to its height and fades in, and closing plays it back (Show's `reveal`;
 * `dom/motion.ts` for code that hides content itself, as section lists do); panels slide in from the side; chevrons
 * turn; hover and press change color; option lists and popovers drop in. Nothing animates on the first render, on a
 * tab switch, for rows of a list, or with reduced motion. Slow motion (×5) to look at them.
 */
export const Motion: StoryObj = {
  render: () =>
    story(() => {
      const slow = observableState(false);
      const panel = observableState(false);
      const theme = observableState('chain-id');
      const card =
        'display: grid; gap: 8px; align-content: start; padding: 12px; border-radius: var(--msp-radius-lg); background: var(--msp-color-panel); box-shadow: inset 0 0 0 1px var(--msp-glass-edge);';
      return (
        <div
          style={{
            display: 'grid',
            gap: '12px',
            padding: '16px',
            maxWidth: '960px',
            '--msp-motion-fast': slow.asReadonly.pipe(map((s) => (s ? '500ms' : ''))),
            '--msp-motion-base': slow.asReadonly.pipe(map((s) => (s ? '900ms' : ''))),
          }}
        >
          <div style="display: flex; gap: 16px; align-items: center;">
            <label style="display: flex; gap: 6px; align-items: center; font-size: var(--msp-font-size-sm);">
              <input type="checkbox" checked={slow.asReadonly} onChange={(e) => slow.set(e.currentTarget.checked)} />
              Slow motion (×5)
            </label>
            <p style={note}>
              --msp-motion-fast 100 ms · --msp-motion-base 220 ms · --msp-ease (0.2, 0, 0, 1) · --msp-ease-reveal (0.4,
              0, 0.2, 1)
            </p>
          </div>
          <div style="display: grid; grid-template-columns: repeat(auto-fill, minmax(260px, 1fr)); gap: 12px;">
            <section style={card}>
              <h3 style={heading}>Fold (trees, groups, sections)</h3>
              <p style={note}>The chevron turns (fast); the content grows in (base) and collapses back.</p>
              <FoldRow label="1TQN" children={['Polymer', 'Ligand', 'Water']} />
              <FoldRow label="Advanced" children={['Quality', 'Smoothness', 'Ignore light']} />
            </section>
            <section style={card}>
              <h3 style={heading}>Panel (the right dock, cards)</h3>
              <p style={note}>Slides in from the side and back; a tab switch inside does not move.</p>
              <Button size="sm" onClick={() => panel.set(!panel.value)}>
                {panel.asReadonly.pipe(map((p) => (p ? 'Close' : 'Open')))}
              </Button>
              <div style="min-height: 120px;">
                <Show when={panel.asReadonly} reveal="side">
                  {() => (
                    <div style="height: 110px; padding: 12px; border-radius: var(--msp-radius-md); background: var(--msp-color-raised); font-size: var(--msp-font-size-sm);">
                      A panel
                    </div>
                  )}
                </Show>
              </div>
            </section>
            <section style={card}>
              <h3 style={heading}>State (hover, press, lists)</h3>
              <p style={note}>Colors change (fast); an option list drops in.</p>
              <div style="display: flex; gap: 8px;">
                <Button size="sm">Hover me</Button>
                <Button size="sm" variant="primary">
                  Press me
                </Button>
              </div>
              <Select
                size="sm"
                aria-label="Color theme"
                options={[
                  ['chain-id', 'Chain ID'],
                  ['element-symbol', 'Element symbol'],
                  ['uniform', 'Uniform'],
                ]}
                value={theme.asReadonly}
                onChange={(v) => theme.set(v)}
              />
            </section>
          </div>
        </div>
      );
    }),
};
