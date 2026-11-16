/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { Show } from '../../../dom/index.js';
import { Heading, Note, page, story } from '../../../_stories/support.js';
import { state } from '../../../state.js';
import { Button } from '../../button/button.js';
import { Icon } from '../../icon/icon.js';
import { Loaded, Loading, Skeleton, SkeletonRows, SkeletonText } from '../skeleton.js';

const meta: Meta = { title: 'Components/Skeleton' };
export default meta;

const Components = ['Polymer', 'Ligand HEM', 'Water', 'Ion SO4'];
const Fields = [
  ['Quality', 'Auto'],
  ['Opacity', '100%'],
  ['Color theme', 'Chain ID'],
];

/**
 * Placeholders in place: a list's rows, a form, a paragraph, a plot. Each sits on the columns and lines of what replaces
 * it (styles §5), so toggling to the loaded content moves nothing.
 */
export const InPlace: StoryObj = {
  name: 'In place',
  render: () =>
    story(() => {
      const loaded = state(false);
      return (
        <div style={page}>
          <div style="display: flex; gap: var(--msui-space-3); align-items: center;">
            <Button size="sm" onClick={() => loaded.set(!loaded.value)}>
              <Show when={loaded} fallback={() => 'Show loaded'}>
                {() => 'Show loading'}
              </Show>
            </Button>
            <Note>
              The loading state fades in after 300ms (a quicker load never shows it); loaded, it crossfades in place.
            </Note>
          </div>
          <div style="display: flex; flex-wrap: wrap; gap: var(--msui-space-5); align-items: start;">
            <div style="width: 240px;">
              <Heading>List</Heading>
              <div class="msui-surface">
                <div class="msui-header">Components</div>
                <Loaded
                  when={loaded}
                  label="Loading components"
                  placeholder={() => (
                    <>
                      <SkeletonRows count={4} />
                    </>
                  )}
                >
                  {() =>
                    Components.map((c) => (
                      <div class="msui-row">
                        <span class="msui-slot">
                          <Icon name="layers" />
                        </span>
                        <span>{c}</span>
                        <span />
                      </div>
                    ))
                  }
                </Loaded>
              </div>
            </div>
            <div style="width: 280px;">
              <Heading>Form</Heading>
              <div class="msui-surface">
                <Loaded
                  when={loaded}
                  label="Loading options"
                  placeholder={() => (
                    <>
                      <div class="msui-form">
                        {Fields.map(() => [
                          <SkeletonText width="64px" />,
                          <Skeleton height="var(--msui-control-height)" />,
                        ])}
                      </div>
                    </>
                  )}
                >
                  {() => (
                    <div class="msui-form">
                      {Fields.map(([label, value]) => [
                        <span class="msui-label">{label}</span>,
                        <span class="msui-value">{value}</span>,
                      ])}
                    </div>
                  )}
                </Loaded>
              </div>
            </div>
            <div style="width: 280px;">
              <Heading>Text and a plot</Heading>
              <div class="msui-surface" style="padding: var(--msui-surface-pad) var(--msui-inset);">
                <Loaded
                  when={loaded}
                  label="Loading the summary"
                  placeholder={() => (
                    <>
                      <SkeletonText lines={2} />
                      <div style="height: var(--msui-space-2)" />
                      <Skeleton height="64px" radius="sm" />
                    </>
                  )}
                >
                  {() => [
                    <p style="margin: 0;">
                      1TQN: cytochrome P450 3A4, 486 residues in one chain, a heme and its ligand, 2.05 Å.
                    </p>,
                    <div style="height: var(--msui-space-2)" />,
                    <div style="height: 64px; border-radius: var(--msui-radius-sm); background: var(--msui-surface-sunken);" />,
                  ]}
                </Loaded>
              </div>
            </div>
          </div>
        </div>
      );
    }),
};

/** The shapes on their own: text in the three sizes, rows with and without icons, blocks and a round one. */
export const Shapes: StoryObj = {
  render: () =>
    story(() => (
      <div style={page}>
        <Loading>
          <div style="display: grid; grid-template-columns: repeat(3, 200px); gap: var(--msui-space-5); align-items: start;">
            <div>
              <Heading>Text: xs, sm, md</Heading>
              <SkeletonText size="xs" lines={2} />
              <SkeletonText size="sm" lines={2} />
              <SkeletonText size="md" lines={2} />
            </div>
            <div>
              <Heading>Rows</Heading>
              <div class="msui-surface">
                <SkeletonRows count={2} />
                <SkeletonRows count={2} round />
                <SkeletonRows count={2} icon={false} />
              </div>
            </div>
            <div style="display: grid; gap: var(--msui-space-3);">
              <Heading>Blocks</Heading>
              <Skeleton height="var(--msui-control-height)" />
              <Skeleton height="80px" radius="lg" />
              <Skeleton width="32px" height="32px" round />
            </div>
          </div>
        </Loading>
      </div>
    )),
};
