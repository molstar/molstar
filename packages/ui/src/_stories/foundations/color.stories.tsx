/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { grid, Heading, Note, page, story, Swatch } from '../support.js';

const meta: Meta = { title: 'Foundations/Color' };
export default meta;

/** The tokens by role: a theme sets the base colors; every shade below is mixed from them (styles §1). */
const Groups: readonly [title: string, note: string, tokens: readonly string[]][] = [
  [
    'Base',
    'What a theme sets: the background, the text, the accent, the shade direction.',
    ['bg', 'fg', 'accent', 'shade'],
  ],
  [
    'Surfaces',
    'Sunken, the surface, raised; a control is raised, hover and active a step each above it.',
    ['surface-sunken', 'surface', 'surface-raised', 'control', 'control-hover', 'control-active'],
  ],
  [
    'Text',
    'Body and values; labels and secondary text; hints, counts, placeholders.',
    ['text', 'text-muted', 'text-subtle'],
  ],
  ['Lines', 'Hairlines only; separation is by surface first.', ['line', 'line-strong', 'focus']],
  [
    'States',
    'The accent for actions and the picked item only; the rest for messages.',
    ['accent', 'accent-muted', 'danger', 'success', 'warning', 'info'],
  ],
];

/** The color tokens of the current theme, by role (switch themes in the toolbar). */
export const Tokens: StoryObj = {
  render: () =>
    story(() => (
      <div style={page}>
        {Groups.map(([title, text, tokens]) => (
          <section style="display: grid; gap: var(--msui-space-3);">
            <Heading>{title}</Heading>
            <Note>{text}</Note>
            <div style={grid}>
              {tokens.map((t) => (
                <Swatch token={`--msui-${t}`} />
              ))}
            </div>
          </section>
        ))}
      </div>
    )),
};

function Layers() {
  return (
    <div class="msui-surface" style="padding-bottom: var(--msui-space-3);">
      <div class="msui-header">Surface</div>
      <div class="msui-form">
        <span class="msui-label">Muted</span>
        <span class="msui-value">
          The text · <span style="color: var(--msui-text-subtle)">subtle</span>
        </span>
      </div>
      {/* A well: sunken, at the surface's padding; the form in it keeps the inset column. */}
      <div style="margin-top: var(--msui-space-2); padding-block: var(--msui-space-2); border-radius: var(--msui-radius-md); background: var(--msui-surface-sunken);">
        <div class="msui-form">
          <span class="msui-label">Sunken</span>
          <span class="msui-control">Control</span>
          <span class="msui-label">Hover</span>
          <span class="msui-control" style="background: var(--msui-control-hover)">
            Control, hovered
          </span>
          <span class="msui-label">Active</span>
          <span class="msui-control" style="background: var(--msui-control-active)">
            Control, pressed
          </span>
        </div>
      </div>
      <div style="display: flex; justify-content: flex-end; padding: var(--msui-space-3) var(--msui-inner) 0;">
        <span style="color: var(--msui-accent); font-weight: var(--msui-weight-strong);">Apply</span>
      </div>
    </div>
  );
}

/** How the surfaces layer: a well inside a surface, controls on it, hover and active above, the accent on the action. */
export const Surfaces: StoryObj = {
  render: () =>
    story(() => (
      <div style={`${page} max-width: 400px;`}>
        <Layers />
      </div>
    )),
};

/** The same surface in every theme: `data-msui-theme` on any element recomputes every derived token below it. */
export const Themes: StoryObj = {
  render: () =>
    story(() => (
      <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(260px, 1fr));">
        {['light', 'dark', 'blue', 'contrast'].map((theme) => (
          <div
            data-msui-theme={theme}
            style="display: grid; gap: var(--msui-space-3); align-content: start; padding: var(--msui-space-5); background: var(--msui-bg); color: var(--msui-text);"
          >
            <b style="text-transform: capitalize;">{theme}</b>
            <Layers />
          </div>
        ))}
      </div>
    )),
};
