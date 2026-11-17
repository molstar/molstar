/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { story } from '../../_stories/support.js';
import { findView } from '../layout.js';
import type { ShellModel } from '../shell-model.js';
import { DemoShell, type DemoOptions } from './demo.js';

// The split stage (shell §2.1): the main area is a region of panes around the app's stage, in the shell's one layout.
// A stage is a view (`stage: true`): beside the app's, below it, in a window; a sequence docks above or over its top; a
// panel docks beside it. Right-click a pane's header (a panel's, or a stage's tab) for its menu; the ⋯ menu of a view
// in the right column offers the stage.

const meta: Meta = { title: 'Stage/Split', parameters: { layout: 'fullscreen' } };
export default meta;

const splitStory = (options: DemoOptions = {}): StoryObj => ({
  render: () => story(() => <DemoShell options={options} />),
});

/** The group of a view's pane (to place another beside it). */
const paneOf = (s: ShellModel, ref: string) => findView(s.layout.value, ref)!.group.id;

/** Two stages side by side: the app's (Stage A) and a stage view (Stage B), each with its own toolbar. */
export const SideBySide = splitStory({ setup: (s) => s.show('stage', { key: 'b' }) });

/** Four stages placed one by one: two columns, each with its own line (drag one: the other column stays). */
export const Quad = splitStory({
  setup: (s) => {
    s.show('stage', { key: 'b' });
    s.show('stage', { key: 'c', target: { kind: 'edge', of: 'content', side: 'bottom' } });
    s.show('stage', { key: 'd', target: { kind: 'edge', of: paneOf(s, 'stage#b'), side: 'bottom' } });
  },
});

/**
 * The same four as a grid (a pane's menu: Arrange as a grid, `view.grid-panes`): equal cells, lined up; drag the
 * line across and both columns follow. Split or close a pane and it is no grid any more.
 */
export const Grid = splitStory({
  setup: (s) => {
    s.show('stage', { key: 'b' });
    s.show('stage', { key: 'c' });
    s.show('stage', { key: 'd' });
    s.arrangeGrid();
  },
});

/** The sequence docked above the stage (the shell's `top: 'split'`): a pane of its own, resizable. */
export const SequenceDocked = splitStory({ setup: (s) => s.show('sequence') });

/** The sequence over the stage's top (`top: 'overlay'`): in the stage's band slot, the canvas running under it. */
export const SequenceOverlay = splitStory({ top: 'overlay', setup: (s) => s.show('sequence') });

/** A panel beside the stage, and the right column as it was. */
export const PanelBeside = splitStory({
  setup: (s) => {
    s.show('home');
    s.show('panel');
  },
});

/** One pane maximized (its menu's Maximize, `view.maximize-pane`): the others away and paused; Escape brings them back. */
export const Maximized = splitStory({
  setup: (s) => {
    s.show('stage', { key: 'b' });
    s.show('panel');
    s.maximize(paneOf(s, 'stage#b'));
  },
});

/** A stage in a window: detached, its own chrome; docked back beside the app's stage, its state is kept. */
export const StageInWindow = splitStory({
  setup: (s) => {
    s.show('stage', { key: 'b' });
    s.float('stage#b');
  },
});

/** One toolbar shared by the panes, over the main area: it acts on the active pane (press a pane; it is ringed). */
export const ToolbarShared = splitStory({ toolbars: 'shared', setup: (s) => s.show('stage', { key: 'b' }) });

/** Mixed: each pane its camera tools, the view tools once over the main area. */
export const ToolbarMixed = splitStory({ toolbars: 'mixed', setup: (s) => s.show('stage', { key: 'b' }) });

/**
 * A toolbar only on the active stage: the layout gives each pane its status (`ctx.pane`: active, hovered, focused), and
 * each stage shows its toolbar from it (`toolbarShown`). Hover a pane, or press it; the others' toolbars fade out.
 */
export const ToolbarActiveOnly = splitStory({
  activeOnly: true,
  setup: (s) => {
    s.show('stage', { key: 'b' });
    s.show('stage', { key: 'c' });
    s.show('stage', { key: 'd' });
    s.arrangeGrid();
  },
});
