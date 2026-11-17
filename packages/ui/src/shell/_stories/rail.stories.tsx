/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { State } from '@molstar/core/reactive';
import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { Heading, Note, page, story } from '../../_stories/support.js';
import type { Child } from '../../dom/index.js';
import { domView } from '../../dom/index.js';
import { ViewEmpty, ViewStack } from '../view-frame.js';
import { type DemoOptions, DemoShell } from './demo.js';

/**
 * The rail on its own (shell §2, §4): the right group's views as a column of icon buttons at the shell's edge (its head,
 * its strip, the app's items by group, its foot), stepping down with the shell's height: normal (32px buttons) from
 * 620px, compact (28) from 520, small (24) below; what still does not fit goes under More, the highest `priority`
 * first (else the foot first, from the end).
 */
const meta: Meta = { title: 'Shell/Rail' };
export default meta;

/** A narrow desktop shell (folded: the main area and the rail), `height` tall. */
function RailBox(props: { label: string; height: number; width?: number; options?: DemoOptions }) {
  return (
    <div style="display: grid; gap: var(--msui-space-2); align-content: start;">
      <Heading>{props.label}</Heading>
      <DemoShell
        options={{
          mode: 'desktop',
          railToggle: false,
          ...props.options,
          width: `${props.width ?? 200}px`,
          height: `${props.height}px`,
        }}
      />
    </div>
  );
}

const row = 'display: flex; flex-wrap: wrap; gap: var(--msui-space-5); align-items: start;';

const stage = (content: () => Child): StoryObj => ({
  render: () => story(() => <div style={page.replace('max-width: 960px;', '')}>{content()}</div>),
});

/** More views than the demo's: a strip long enough to overflow a short rail. */
function moreViews(options: DemoOptions['setup'] = () => {}): DemoOptions['setup'] {
  return (s) => {
    s.registerViews(
      (
        [
          ['measure', 'Measurements', 'ruler', 4],
          ['camera', 'Camera', 'camera', 5],
          ['palette', 'Colors', 'palette', 6],
        ] as const
      ).map(([id, title, icon, order]) => ({
        id,
        title,
        icon,
        order,
        render: domView(() => (
          <ViewStack>
            <ViewEmpty>{title}: a strip view.</ViewEmpty>
          </ViewStack>
        )),
      })),
    );
    options?.(s);
  };
}

/**
 * The three height steps, side by side: normal (≥ 620), compact (≥ 520) and small (below). The buttons shrink, the
 * captioned one (⌘K) keeps its caption; the sections (head: Home; strip; foot: Settings, Help; tools: Theme; the corner,
 * under Hide: Commands) keep a hairline between them.
 */
export const Steps = stage(() => (
  <div style={row}>
    <RailBox label="Normal (680)" height={680} />
    <RailBox label="Compact (560)" height={560} />
    <RailBox label="Small (420)" height={420} />
  </div>
));

/**
 * Too short even at the small step: what does not fit goes under More (⋯, above the foot), the foot's and the strip's
 * last first; More's menu lists them. Press More to see it.
 */
export const Overflow = stage(() => (
  <div style={row}>
    <RailBox label="Small, everything fits (420)" height={420} options={{ setup: moreViews() }} />
    <RailBox label="More (300)" height={300} options={{ setup: moreViews() }} />
    <RailBox label="More (240)" height={240} options={{ setup: moreViews() }} />
    <RailBox label="More (180)" height={180} options={{ setup: moreViews() }} />
  </div>
));

/**
 * The buttons' states: the open view's button selected (a press on it folds the column), a pressed app item (Pin, a
 * toggle), a captioned item (⌘K), and an item of another group (a hairline before it).
 */
export const States = stage(() => (
  <>
    <div style={row}>
      <RailBox
        label="Layers open, Pin pressed"
        height={560}
        width={560}
        options={{
          setup: (s) => {
            s.registerRailItems([
              { id: 'pin', icon: 'pin', title: 'Pin the view', group: 'app', pressed: new State(true), run: () => {} },
              { id: 'search', icon: 'search', title: 'Search', group: 'find', caption: '⌘F', run: () => {} },
            ]);
            s.show('layers');
          },
        }}
      />
      <RailBox label="Folded" height={560} />
    </div>
    <Note>Right-click a rail button (or Shift+F10 on it) for its menu: Hide, Show a hidden view, Reset layout.</Note>
  </>
));
