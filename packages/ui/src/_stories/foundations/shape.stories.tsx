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
  ['--msui-radius-md', 'Controls and rows inside a surface: the surface radius minus its padding'],
  ['--msui-radius-lg', 'Surfaces: panels, cards, popovers'],
] as const;

const Elevation = [
  ['--msui-edge', 'A hairline around a surface that sits on its own color'],
  ['--msui-shadow-float', 'What floats over content: a popover, a menu'],
  ['--msui-shadow-overlay', 'What covers the page: a dialog'],
] as const;

const squircles = typeof CSS !== 'undefined' && CSS.supports('corner-shape', 'squircle');

/** Radii by role, and elevation only for what floats (styles §6). */
export const RadiiAndElevation: StoryObj = {
  name: 'Radii and elevation',
  render: () =>
    story(() => (
      <div style={page}>
        <Note>
          {squircles
            ? 'This browser draws squircles: every corner here is one (--msui-corner-shape), with radii grown by half.'
            : 'This browser draws round corners; where corner-shape is supported, every corner becomes a squircle.'}
        </Note>
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

/** A parent of radius `outer` with a child inset by `inset`, the child's radius given. */
function NestedPair(props: {
  title: string;
  note: string;
  outer: string;
  inset: string;
  inner: string;
  good?: boolean;
}) {
  return (
    <div style="display: grid; gap: var(--msui-space-2); width: 200px; align-content: start;">
      <div
        style={{
          padding: props.inset,
          borderRadius: props.outer,
          background: 'var(--msui-surface-raised)',
          boxShadow: 'inset 0 0 0 1px var(--msui-line)',
        }}
      >
        <div style={{ height: '72px', borderRadius: props.inner, background: 'var(--msui-control-active)' }} />
      </div>
      <b style={{ color: props.good ? 'var(--msui-success)' : 'var(--msui-text)' }}>{props.title}</b>
      <span style="color: var(--msui-text-muted)">{props.note}</span>
    </div>
  );
}

/**
 * Nested corners are concentric: inset by `p` inside a radius `R`, the inner radius is `R − p`, so both curves share
 * their centers and the gap is even around the corner (styles §6). `--msui-radius-md` is derived that way from the
 * surface, so rows and controls in a surface are concentric by construction, with round corners and with squircles.
 */
export const NestedCorners: StoryObj = {
  name: 'Nested corners',
  render: () =>
    story(() => (
      <div style={page}>
        <Heading>One level: R − p</Heading>
        <div style="display: flex; flex-wrap: wrap; gap: var(--msui-space-5);">
          <NestedPair
            title="Concentric"
            note="Inner radius = outer − inset: an even gap around the corner."
            outer="20px"
            inset="8px"
            inner="12px"
            good
          />
          <NestedPair
            title="Same radius"
            note="The inner corner bulges toward the outer one: the gap pinches at the corner."
            outer="20px"
            inset="8px"
            inner="20px"
          />
          <NestedPair
            title="Too small"
            note="A wedge of extra space opens at the corner."
            outer="20px"
            inset="8px"
            inner="2px"
          />
        </div>
        <Heading>In the system</Heading>
        <Note>
          A surface (lg) pads by --msui-surface-pad; its rows and controls take md = lg − pad. Another level composes
          the same way: max(0px, calc(parent radius − inset)); an inset of the radius or more gives a square (or sm)
          corner.
        </Note>
        <div style="display: flex; flex-wrap: wrap; gap: var(--msui-space-5); align-items: start;">
          <div class="msui-surface" style="width: 240px;">
            {['Polymer', 'Ligand', 'Water'].map((l, i) => (
              <div class="msui-row msui-row-plain" style={i === 1 ? 'background: var(--msui-control-hover)' : ''}>
                <span>{l}</span>
              </div>
            ))}
          </div>
          <div class="msui-surface" style="width: 240px; padding-block: var(--msui-surface-pad);">
            <div style="padding: var(--msui-surface-pad); border-radius: var(--msui-radius-md); background: var(--msui-surface-sunken);">
              <div style="height: 48px; border-radius: max(0px, calc(var(--msui-radius-md) - var(--msui-surface-pad))); background: var(--msui-control-active);" />
            </div>
          </div>
        </div>
        <Heading>Corner shape</Heading>
        <Note>
          {squircles
            ? 'The same surface with round corners and with squircles (--msui-corner-shape on a subtree).'
            : 'This browser has no corner-shape: both look round here.'}
        </Note>
        <div style="display: flex; gap: var(--msui-space-5);">
          {(['round', 'squircle'] as const).map((shape) => (
            <div style={{ '--msui-corner-shape': shape, display: 'grid', gap: 'var(--msui-space-2)', width: '240px' }}>
              <b>{shape}</b>
              <div class="msui-surface">
                {['Polymer', 'Ligand'].map((l, i) => (
                  <div class="msui-row msui-row-plain" style={i === 0 ? 'background: var(--msui-control-hover)' : ''}>
                    <span>{l}</span>
                  </div>
                ))}
              </div>
            </div>
          ))}
        </div>
      </div>
    )),
};
