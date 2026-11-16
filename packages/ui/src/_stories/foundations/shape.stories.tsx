/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { Heading, Note, page, story } from '../support.js';

const meta: Meta = { title: 'Foundations/Shape' };
export default meta;

const Radii = [
  ['--msui-radius-sm', 'Inline marks: code, chips in text'],
  ['--msui-radius-md', 'Controls and rows inside a surface'],
  ['--msui-radius-lg', 'Surfaces: panels, cards, popovers'],
] as const;

const Elevation = [
  ['--msui-edge', 'A hairline around a surface that sits on its own color'],
  ['--msui-shadow-float', 'What floats over content: a popover, a menu'],
  ['--msui-shadow-overlay', 'What covers the page: a dialog'],
] as const;

/**
 * Radii by role, nested inward: a surface's radius is larger than the rows' inside it, so their corners follow its
 * curve at the inset (styles §6). Elevation is for what floats; what sits in the page has none.
 */
export const RadiiAndElevation: StoryObj = {
  name: 'Radii and elevation',
  render: () =>
    story(() => (
      <div style={page}>
        <Heading>Radii</Heading>
        <div style="display: flex; gap: var(--msui-space-5);">
          {Radii.map(([t, use]) => (
            <div style="display: grid; gap: var(--msui-space-2); width: 160px;">
              <div style={{ height: '64px', borderRadius: `var(${t})`, background: 'var(--msui-surface-raised)' }} />
              <code style="font-size: var(--msui-text-xs);">{t}</code>
              <span style="color: var(--msui-text-muted)">{use}</span>
            </div>
          ))}
        </div>
        <Heading>Nested</Heading>
        <Note>
          A surface (lg) with rows (md) inside, a padding apart: the inner corners stay inside the outer curve.
        </Note>
        <div class="msui-surface" style="width: 280px;">
          {['Polymer', 'Ligand', 'Water'].map((l, i) => (
            <div class="msui-row msui-row-plain" style={i === 1 ? 'background: var(--msui-control-hover)' : ''}>
              <span>{l}</span>
            </div>
          ))}
        </div>
        <Heading>Elevation</Heading>
        <div style="display: flex; gap: var(--msui-space-5); padding: var(--msui-space-4) 0;">
          {Elevation.map(([t, use]) => (
            <div style="display: grid; gap: var(--msui-space-2); width: 180px;">
              <div
                style={{
                  height: '64px',
                  borderRadius: 'var(--msui-radius-lg)',
                  background: 'var(--msui-surface)',
                  boxShadow: `var(${t})`,
                }}
              />
              <code style="font-size: var(--msui-text-xs);">{t}</code>
              <span style="color: var(--msui-text-muted)">{use}</span>
            </div>
          ))}
        </div>
      </div>
    )),
};
