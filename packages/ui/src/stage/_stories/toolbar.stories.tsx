/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { Note, story } from '../../_stories/support.js';
import { onCleanup } from '../../dom/index.js';
import { state } from '../../state.js';
import { Stage } from '../stage.js';
import { StageModel } from '../stage-model.js';
import { StageToolbar } from '../toolbar.js';
import { DemoCanvas, demoGroups, Pill, StageBox } from './demo.js';

const meta: Meta = { title: 'Stage/Toolbar' };
export default meta;

function Demo(props: { collapsed: 'auto' | boolean; width?: string; height?: string; note: string }) {
  const stage = new StageModel({ empty: false });
  onCleanup(() => stage.dispose());
  const ran = state('nothing yet');
  return (
    <div>
      <div style="padding: var(--msui-space-3) var(--msui-space-4) 0">
        <Note>{props.note}</Note>
      </div>
      <StageBox
        {...(props.width ? { width: props.width } : {})}
        height={props.height ?? 'calc(100vh - 12 * var(--msui-unit))'}
      >
        <Stage
          model={stage}
          canvas={<DemoCanvas stage={stage} />}
          status={<Pill>ran: {ran}</Pill>}
          corner={
            <StageToolbar groups={demoGroups((what) => ran.set(what))} stage={stage} collapsed={props.collapsed} />
          }
        />
      </StageBox>
    </div>
  );
}

/** Expanded: menus open upward from the corner (Reset's projection, Lighting with a Style submenu); a disabled reason. */
export const Expanded: StoryObj = {
  render: () =>
    story(() => (
      <Demo
        collapsed={false}
        note="Hover Lighting: its menu opens upward, as tall as the room above; Style opens a submenu. Fullscreen is greyed with its reason."
      />
    )),
};

/** Folded: one button whose menu lists the rows (toggles checkable, menus as submenus led by their own action). */
export const Collapsed: StoryObj = {
  render: () =>
    story(() => (
      <Demo collapsed={true} note="The folded toolbar: hover or tap the wrench for every tool in one menu." />
    )),
};

/** Narrow: the row scrolls across where the stage is narrower than it (not folded). */
export const NarrowRow: StoryObj = {
  render: () =>
    story(() => (
      <Demo
        collapsed={false}
        width="200px"
        height="420px"
        note="A row wider than the stage scrolls across (the folding is off here)."
      />
    )),
};

/** Short: little room above, so the menus are short and scroll inside. */
export const ShortStage: StoryObj = {
  render: () =>
    story(() => (
      <Demo collapsed={false} width="720px" height="160px" note="A short stage: menus are bounded by the room above." />
    )),
};
