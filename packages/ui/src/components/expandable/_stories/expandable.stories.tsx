/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { Heading, story } from '../../../_stories/support.js';
import { type Child, onCleanup } from '../../../dom/index.js';
import { ParamsFormModel } from '../../../forms/params-form.js';
import { state } from '../../../state.js';
import { Button } from '../../button/button.js';
import { ParamsForm } from '../../params/params-form.js';
import { installTooltips } from '../../tooltip/tooltip.js';
import { Expandable, ExpandableList } from '../expandable.js';

const meta: Meta = { title: 'Components/Expandable' };
export default meta;

function Panel(props: { width?: number; children: Child }) {
  return (
    <div
      ref={(el: HTMLElement) => onCleanup(installTooltips(el))}
      class="msui-surface"
      style={`width: ${props.width ?? 340}px; margin: var(--msui-space-4); padding-block: var(--msui-space-1);`}
    >
      {props.children}
    </div>
  );
}

const background = {
  color: PD.Color(0xfcfbf9 as any, { label: 'Color' }),
  transparent: PD.Boolean(false, { label: 'Transparent' }),
  environment: PD.Select('off', [
    ['off', 'Off'],
    ['studio', 'Studio'],
    ['sky', 'Sky'],
  ]),
};

const lighting = {
  ambient: PD.Numeric(0.4, { min: 0, max: 1, step: 0.01 }),
  key: PD.Numeric(0.6, { min: 0, max: 1, step: 0.01 }),
  shadows: PD.Boolean(true),
};

const camera = {
  mode: PD.Select('perspective', [
    ['perspective', 'Perspective'],
    ['orthographic', 'Orthographic'],
  ]),
  fov: PD.Numeric(45, { min: 10, max: 90, step: 1 }, { label: 'Field of view' }),
};

const form = <P extends PD.Params>(params: P) => (
  <ParamsForm form={new ParamsFormModel(params, PD.getDefaultValues(params))} />
);

/**
 * Panels stacked: each opens to a well holding its header and its form, the form's labels on the header's column.
 * Closed, a header hovers as a row. Opening grows the body in; closing plays it back, the chevron turning at once.
 */
export const Stacked: StoryObj = {
  render: () =>
    story(() => (
      <Panel>
        <ExpandableList>
          <Expandable title="Background" open>
            {() => form(background)}
          </Expandable>
          <Expandable title="Lighting">{() => form(lighting)}</Expandable>
          <Expandable title="Camera">{() => form(camera)}</Expandable>
        </ExpandableList>
      </Panel>
    )),
};

/**
 * Held by the host: the open state is a `State` the host keeps (to remember it, to open one from elsewhere); a trailing
 * count in the header; the flat variant inside a surface that already is one (an inspector).
 */
export const Variants: StoryObj = {
  render: () =>
    story(() => {
      const open = state(false);
      return (
        <Panel>
          <Heading>Held by the host</Heading>
          <div>
            <Button size="sm" onClick={() => open.set(!open.value)}>
              Toggle from outside
            </Button>
          </div>
          <ExpandableList>
            <Expandable title="Measurements" open={open} trailing="3">
              {() => form(lighting)}
            </Expandable>
          </ExpandableList>
          <Heading>Flat</Heading>
          <ExpandableList>
            <Expandable title="Background" variant="flat" open>
              {() => form(background)}
            </Expandable>
            <Expandable title="Camera" variant="flat">
              {() => form(camera)}
            </Expandable>
          </ExpandableList>
          <Heading>A long title</Heading>
          <Expandable title="Electron density around the selected ligand and its neighbours">
            {() => form(camera)}
          </Expandable>
        </Panel>
      );
    }),
};
