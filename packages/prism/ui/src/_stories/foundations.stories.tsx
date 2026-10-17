/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { Button } from '../components/button/button.js';
import { Icon } from '../components/icon/icon.js';
import { Select } from '../components/select/select.js';
import { Show } from '../dom/control-flow.js';
import { observableState } from '../reactive/observable-state.js';
import { story } from '../testing/story-root.js';

import { derive } from '../reactive/state.js';
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
    <div style="display: grid; gap: var(--msp-space-2); font-size: var(--msp-font-size-xs)">
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

const grid = 'display: grid; grid-template-columns: repeat(auto-fill, minmax(132px, 1fr)); gap: var(--msp-space-4);';
const section = 'display: grid; gap: var(--msp-space-4); padding: var(--msp-space-5); max-width: 960px;';
const heading = 'margin: var(--msp-space-3) 0 0; font-size: var(--msp-font-size-md);';
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
      <div style="padding: var(--msp-space-5); display: grid; gap: var(--msp-space-3); max-width: 360px; font-size: var(--msp-font-size-sm)">
        <div style="padding: var(--msp-space-4); border-radius: var(--msp-radius-lg); background: var(--msp-color-panel); display: grid; gap: var(--msp-space-3)">
          <b>Panel</b>
          <span style="color: var(--msp-color-text-muted)">A label, muted · a hint, subtle</span>
          <div style="padding: var(--msp-space-3); border-radius: var(--msp-radius-md); background: var(--msp-color-well); display: grid; gap: var(--msp-space-2)">
            <b>Well (an open action)</b>
            <div style="height: 24px; border-radius: var(--msp-radius-md); background: var(--msp-color-control); padding: 0 var(--msp-space-3); display: flex; align-items: center">
              Control
            </div>
            <div style="height: 24px; border-radius: var(--msp-radius-md); background: var(--msp-color-hover); padding: 0 var(--msp-space-3); display: flex; align-items: center">
              Hover
            </div>
            <div style="height: 24px; border-radius: var(--msp-radius-md); background: var(--msp-color-pressed); padding: 0 var(--msp-space-3); display: flex; align-items: center">
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
        {['light', 'dark', 'blue', 'contrast'].map((theme) => (
          <div
            data-msp-theme={theme}
            style="background: var(--msp-color-bg); color: var(--msp-color-text); padding: var(--msp-space-5); display: grid; gap: var(--msp-space-3); font-size: var(--msp-font-size-sm)"
          >
            <b style="text-transform: capitalize">{theme}</b>
            <div style="background: var(--msp-color-panel); padding: var(--msp-space-3); border-radius: var(--msp-radius-md); display: grid; gap: 6px">
              <span>
                Panel <span style="color: var(--msp-color-text-muted)">muted</span>{' '}
                <span style="color: var(--msp-color-text-subtle)">subtle</span>
              </span>
              <div style="background: var(--msp-color-well); padding: 6px var(--msp-space-3); border-radius: var(--msp-radius-md); display: flex; justify-content: space-between">
                <span style="background: var(--msp-color-control); padding: 0 var(--msp-space-3); border-radius: var(--msp-radius-md)">
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

/** Type: the sizes (11, 12, 14, 16 px; panels read at 12) and the monospace face (sequences, code). */
export const Type: StoryObj = {
  render: () =>
    story(() => (
      <div style={section}>
        {['xs', 'sm', 'md', 'lg'].map((s) => (
          <div style={{ fontSize: `var(--msp-font-size-${s})` }}>
            --msp-font-size-{s}: The quick brown fox jumps over the lazy dog
          </div>
        ))}
        <code style="font-family: var(--msp-font-mono)">--msp-font-mono: MALYGTHSHGLFKKLGIPGPTPLP</code>
      </div>
    )),
};

/** The spacing scale and what each step is for (theme.css). */
const Scale: [token: string, use: string][] = [
  ['--msp-space-1', 'Between rows of a list; between a label and its hint'],
  ['--msp-space-2', "Inside a row (icon to text); between controls in a group; a panel's inset"],
  ['--msp-space-3', "Between groups in a form; a card's inner gap; where text starts in a control"],
  ['--msp-space-4', "A card's padding; between sections"],
  ['--msp-space-5', 'Between blocks of a page'],
  ['--msp-space-6', 'Empty states; large separations'],
];

/** The rules named after the scale: what a value means, so a change moves everything that follows the rule. */
const Rules: [token: string, use: string][] = [
  ['--msp-inset', 'Where text starts inside a control or a row: fields, selects, buttons and tree rows line up'],
  ['--msp-panel-inset', "A panel's content inset at its sides and bottom (rows inside pad to --msp-inset)"],
  ['--msp-canvas-edge', 'The margin of everything floating over the canvas'],
];

const Sizes: [token: string, use: string][] = [
  ['--msp-row-height', 'A list or tree row'],
  ['--msp-control-height', 'Fields, selects, buttons'],
  ['--msp-label-width', "A form's label column"],
  ['--msp-icon-sm', 'Small icons (a chevron), for indents after them'],
  ['--msp-icon-md', "A row's kind icon"],
];

function TokenRow(props: { token: string; use: string; bar?: boolean }) {
  return (
    <div style="display: grid; grid-template-columns: 180px 48px 32px 1fr; gap: var(--msp-space-4); align-items: center; min-height: var(--msp-row-height)">
      <code style="font-size: var(--msp-font-size-xs); color: var(--msp-color-text)">{props.token}</code>
      <span style="color: var(--msp-color-text-muted); font-variant-numeric: tabular-nums" />
      <div
        style={{
          width: props.bar ? `var(${props.token})` : '0',
          height: '12px',
          borderRadius: 'var(--msp-radius-sm)',
          background: 'var(--msp-color-accent)',
        }}
      />
      <span style={note}>{props.use}</span>
    </div>
  );
}

/** A token's value as the page resolves it. */
function resolved(token: string) {
  const probe = document.createElement('div');
  probe.style.width = `var(${token})`;
  document.body.append(probe);
  const px = `${Math.round(probe.getBoundingClientRect().width * 10) / 10}px`;
  probe.remove();
  return px;
}

function TokenTable(props: { rows: [string, string][]; bars?: boolean }) {
  return (
    <div style="display: grid">
      {props.rows.map(([token, use]) => {
        const row = (<TokenRow token={token} use={use} bar={props.bars} />) as HTMLElement;
        row.children[1].textContent = resolved(token);
        return row;
      })}
    </div>
  );
}

/**
 * Spacing: every padding, margin and gap is a step of the scale or a rule named after it (styles:check rejects other
 * px there, but 1px hairlines and values marked `off-scale` with their reason). The panel below is laid out with them:
 * the panel's inset, rows' text at the inset, rows a step apart.
 */
export const Spacing: StoryObj = {
  render: () =>
    story(() => (
      <div style={section}>
        <h3 style={heading}>Scale</h3>
        <TokenTable rows={Scale} bars />
        <h3 style={heading}>Rules</h3>
        <TokenTable rows={Rules} bars />
        <h3 style={heading}>Sizes</h3>
        <TokenTable rows={Sizes} />
        <h3 style={heading}>In a panel</h3>
        <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(260px, 1fr)); gap: var(--msp-space-5); align-items: start">
          <div style="display: grid; gap: var(--msp-space-1); padding: 0 var(--msp-panel-inset) var(--msp-panel-inset); border-radius: var(--msp-radius-lg); background: var(--msp-color-panel); box-shadow: inset 0 0 0 1px var(--msp-glass-edge); font-size: var(--msp-font-size-sm)">
            <div style="display: flex; align-items: center; height: 40px; padding: 0 var(--msp-space-4); font-weight: 600">
              Panel
            </div>
            {['Polymer', 'Ligand', 'Water'].map((label, i) => (
              <div
                style={{
                  display: 'flex',
                  alignItems: 'center',
                  gap: 'var(--msp-space-2)',
                  height: 'var(--msp-row-height)',
                  padding: '0 var(--msp-inset)',
                  borderRadius: 'var(--msp-radius-md)',
                  background: i === 1 ? 'var(--msp-color-hover)' : 'transparent',
                }}
              >
                {label}
                <span style="color: var(--msp-color-text-muted)">a description</span>
              </div>
            ))}
            <div style="display: flex; align-items: center; height: var(--msp-control-height); padding: 0 var(--msp-inset); border-radius: var(--msp-radius-md); background: var(--msp-color-control)">
              A field: its text at the same inset
            </div>
          </div>
          <p style={note}>
            The panel insets its content by --msp-panel-inset (4px); rows and controls start their text at --msp-inset
            (8px), so it lines up with the header's 12px; rows are --msp-space-1 apart.
          </p>
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
        style="display: flex; align-items: center; gap: 6px; width: 100%; height: 26px; padding: 0 var(--msp-space-2); border-radius: var(--msp-radius-md); color: var(--msp-color-text); font-size: var(--msp-font-size-sm);"
        aria-expanded={open.asReadonly}
        onClick={() => open.set(!open.value)}
      >
        <span
          style={{
            display: 'inline-flex',
            transition: 'transform var(--msp-motion-fast) var(--msp-ease)',
            transform: derive(open.asReadonly, (o) => (o ? 'rotate(90deg)' : 'none')),
          }}
        >
          <Icon name="chevron-right" size={14} />
        </span>
        {props.label}
      </button>
      <Show when={open.asReadonly} reveal>
        {() => (
          <div style="display: grid; padding-left: var(--msp-space-6);">
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

/** The principles of motion, as the Motion story lists them (ui/_spec/shell.md, Motion). */
const MotionPrinciples: readonly [string, string][] = [
  [
    'Pick the move by what must not change',
    'Grow (to its height) for content in a flow: folds, sections. Slide for panels. Fade, opacity only, for what must keep its size: a band drawn on a canvas, the 3D canvas itself (growing would lay it out again every frame, and it jitters).',
  ],
  [
    'Replace by a scoped view transition',
    "One view replacing another (a right panel by another) swaps as a view transition scoped to its container where the browser has them (`Show` `swap`, `swapIn`). Never the document's: it snapshots the page and the 3D canvas stands still. Elsewhere, at once.",
  ],
  [
    'Calm before the content',
    "The canvas starts quiet and is quiet while the empty canvas's mark and splash are up (no axes, nothing half drawn); what loads comes in gradually over the canvas's own color.",
  ],
  [
    'Guard against flashes',
    'Only a settled state animates in: the mark waits for the canvas to be empty and the viewer idle (work restarts the wait), so a scene cleared to load the next never flashes it, and a load it never showed for appears at once.',
  ],
  [
    'Cheap to run',
    'transform and opacity only, never layout in a frame; nothing animates on the first render, on a tab switch inside a panel, for rows of a list; nothing at all under reduced motion.',
  ],
];

/**
 * Motion: what moves, and how (`--msp-motion-fast` for state changes, `--msp-motion-base` for content that opens,
 * `--msp-motion-fade` for what only fades, `--msp-motion-reveal` for the canvas coming in, `--msp-ease`, and
 * `--msp-ease-reveal` for opening and closing). Content that opens grows to its height and fades in, and closing plays
 * it back (Show's `reveal`; `dom/motion.ts` for code that hides content itself, as section lists do); panels slide in
 * from the side; a band only fades (`reveal="fade"`); one panel replacing another swaps by a scoped view transition
 * (`swap`); the canvas is quiet under the empty mark and comes in on a load; chevrons turn; hover and press change
 * color; option lists and popovers drop in. Nothing animates on the first render, on a tab switch inside a panel, for
 * rows of a list, or with reduced motion. Slow motion (×5) to look at them.
 */
export const Motion: StoryObj = {
  render: () =>
    story(() => {
      const slow = observableState(false);
      const panel = observableState(false);
      const theme = observableState('chain-id');
      const band = observableState(true);
      const tab = observableState<'home' | 'structure' | 'selection'>('home');
      const quiet = observableState(true);
      let swapScope: HTMLElement | undefined;
      const card =
        'display: grid; gap: var(--msp-space-3); align-content: start; padding: var(--msp-space-4); border-radius: var(--msp-radius-lg); background: var(--msp-color-panel); box-shadow: inset 0 0 0 1px var(--msp-glass-edge);';
      return (
        <div
          style={{
            display: 'grid',
            gap: 'var(--msp-space-4)',
            padding: 'var(--msp-space-5)',
            maxWidth: '960px',
            '--msp-motion-fast': derive(slow.asReadonly, (s) => (s ? '500ms' : '')),
            '--msp-motion-base': derive(slow.asReadonly, (s) => (s ? '900ms' : '')),
            '--msp-motion-fade': derive(slow.asReadonly, (s) => (s ? '1800ms' : '')),
            '--msp-motion-reveal': derive(slow.asReadonly, (s) => (s ? '4500ms' : '')),
          }}
        >
          <section style="display: grid; gap: var(--msp-space-3);">
            <h3 style={heading}>Principles</h3>
            <dl style="display: grid; grid-template-columns: minmax(180px, max-content) 1fr; gap: var(--msp-space-2) var(--msp-space-4); margin: 0;">
              {MotionPrinciples.map(([rule, why]) => [
                <dt style="font-size: var(--msp-font-size-sm); font-weight: 600;">{rule}</dt>,
                <dd style={`${note}; margin: 0;`}>{why}</dd>,
              ])}
            </dl>
          </section>
          <div style="display: flex; gap: var(--msp-space-5); align-items: center;">
            <label style="display: flex; gap: 6px; align-items: center; font-size: var(--msp-font-size-sm);">
              <input type="checkbox" checked={slow.asReadonly} onChange={(e) => slow.set(e.currentTarget.checked)} />
              Slow motion (×5)
            </label>
            <p style={note}>
              --msp-motion-fast 100 ms · --msp-motion-base 220 ms · --msp-motion-fade 360 ms · --msp-motion-reveal 900
              ms · --msp-ease (0.2, 0, 0, 1) · --msp-ease-reveal (0.4, 0, 0.2, 1)
            </p>
          </div>
          <div style="display: grid; grid-template-columns: repeat(auto-fill, minmax(260px, 1fr)); gap: var(--msp-space-4);">
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
                {derive(panel.asReadonly, (p) => (p ? 'Close' : 'Open'))}
              </Button>
              <div style="min-height: 120px;">
                <Show when={panel.asReadonly} reveal="side">
                  {() => (
                    <div style="height: 110px; padding: var(--msp-space-4); border-radius: var(--msp-radius-md); background: var(--msp-color-raised); font-size: var(--msp-font-size-sm);">
                      A panel
                    </div>
                  )}
                </Show>
              </div>
            </section>
            <section style={card}>
              <h3 style={heading}>State (hover, press, lists)</h3>
              <p style={note}>Colors change (fast); an option list drops in.</p>
              <div style="display: flex; gap: var(--msp-space-3);">
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
            <section style={card}>
              <h3 style={heading}>Fade (a band on the canvas)</h3>
              <p style={note}>Only its opacity moves: drawn on a canvas, it keeps its size (fade).</p>
              <Button size="sm" onClick={() => band.set(!band.value)}>
                {derive(band.asReadonly, (b) => (b ? 'Hide the band' : 'Show the band'))}
              </Button>
              <div style="height: 48px;">
                <Show when={band.asReadonly} reveal="fade">
                  {() => (
                    <div style="height: 48px; border-radius: var(--msp-radius-md); background: repeating-linear-gradient(90deg, var(--msp-color-raised) 0 10px, var(--msp-color-panel) 10px 14px); font: var(--msp-font-size-xs) var(--msp-font-mono); display: flex; align-items: center; padding-inline: var(--msp-space-3);">
                      MKTAYIAKQRQISFVKSHFSRQ
                    </div>
                  )}
                </Show>
              </div>
            </section>
            <section style={card}>
              <h3 style={heading}>Swap (one panel for another)</h3>
              <p style={note}>A scoped view transition where the browser has one; elsewhere at once.</p>
              <div style="display: flex; gap: var(--msp-space-2);">
                {(['home', 'structure', 'selection'] as const).map((t) => (
                  <Button size="sm" onClick={() => tab.set(t)}>
                    {t[0].toUpperCase() + t.slice(1)}
                  </Button>
                ))}
              </div>
              <div
                ref={(el) => {
                  swapScope = el;
                }}
                style="min-height: 90px;"
              >
                <Show when={tab.asReadonly} swap={() => swapScope}>
                  {(t) => (
                    <div style="padding: var(--msp-space-4); border-radius: var(--msp-radius-md); background: var(--msp-color-raised); font-size: var(--msp-font-size-sm);">
                      {`The ${t} panel`}
                    </div>
                  )}
                </Show>
              </div>
            </section>
            <section style={card}>
              <h3 style={heading}>Reveal (the canvas)</h3>
              <p style={note}>Quiet under the empty mark; what loads comes in over the canvas's own color.</p>
              <Button size="sm" onClick={() => quiet.set(!quiet.value)}>
                {derive(quiet.asReadonly, (q) => (q ? 'Load' : 'Clear'))}
              </Button>
              <div style="height: 90px; border-radius: var(--msp-radius-md); background: var(--msp-color-bg); overflow: hidden;">
                <div
                  style={{
                    height: '100%',
                    display: 'grid',
                    placeItems: 'center',
                    transition: 'opacity var(--msp-motion-reveal) var(--msp-ease)',
                    opacity: derive(quiet.asReadonly, (q) => (q ? '0' : '1')),
                  }}
                >
                  <Icon name="cartoon" />
                </div>
              </div>
            </section>
          </div>
        </div>
      );
    }),
};
