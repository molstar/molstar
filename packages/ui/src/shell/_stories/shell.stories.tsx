/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { State } from '@molstar/core/reactive';
import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { story } from '../../_stories/support.js';
import { findView, type Rect, RightGroup } from '../layout.js';
import type { ShellModel } from '../shell-model.js';
import type { ShellBusy } from '../shell-view.js';
import { DemoShell, type DemoOptions } from './demo.js';

const meta: Meta = { title: 'Shell/Shell', parameters: { layout: 'fullscreen' } };
export default meta;

const shellStory = (options: DemoOptions = {}): StoryObj => ({
  render: () => story(() => <DemoShell options={options} />),
});

/** The right group open on Home, the rail beside it: the shell as an app starts. */
export const Open = shellStory({ setup: (s) => s.show('home') });

/** Nothing loaded yet: the stage shows the glass mark (after 700 ms); drop a file on the shell to "load" it. */
export const EmptyStage = shellStory({ empty: true, dropLabel: 'Drop files to open them' });

/** Folded to the rail: the main area takes the room; a press on a rail button opens its view. */
export const Folded = shellStory();

/** The rail and the docked column hidden: the Mol* star in the main area's corner brings them back. */
export const RailHidden = shellStory({ setup: (s) => s.toggleRail(false) });

/** The right column split: Inspect moved below Layers (its ⋯ menu: Split below), Settings into the lower group, its tabs. */
export const SplitColumn = shellStory({
  setup: (s) => {
    s.show('layers');
    s.splitBelow('inspect');
    const lower = findView(s.layout.value, 'inspect')?.group.id;
    if (lower) s.move('settings', { kind: 'group', id: lower });
  },
});

/** Views popped out into windows: two instances of a details view side by side, a docked view floated, a window of two tabs. */
export const Windows = shellStory({
  setup: (s) => {
    s.show('home');
    s.show('details', { key: 'a', subject: 'Item A' });
    s.show('details', { key: 'b', subject: 'Item B' });
    s.float('inspect');
    const inspect = windowOf(s, 'inspect');
    if (inspect) s.move('help', { kind: 'group', id: inspect });
    s.select('inspect');
    place(s, 'details#a', { x: 24, y: 24, width: 320, height: 240 });
    place(s, 'details#b', { x: 360, y: 24, width: 320, height: 240 });
    place(s, 'inspect', { x: 24, y: 288, width: 340, height: 300 });
  },
});

/** The window group holding `ref`. */
function windowOf(s: ShellModel, ref: string): string | undefined {
  const at = findView(s.layout.value, ref);
  return at?.window ? at.group.id : undefined;
}

/** The window of `ref` to `rect` (in the main area). */
function place(s: ShellModel, ref: string, rect: Rect) {
  const id = windowOf(s, ref);
  if (id) s.setWindowRect(id, rect);
}

/** A card in the right column's place: Export, its sections from contributors, opened at one. */
export const Card = shellStory({ setup: (s) => s.show('export', { section: 'data' }) });

/** A page over the main area. */
export const Page = shellStory({ setup: (s) => s.show('report') });

/** A modal over the whole shell: Tab stays in it, the rest is inert. */
export const Modal = shellStory({ setup: (s) => s.show('confirm') });

/** A task holding the app: its progress and Cancel. */
export const Busy: StoryObj = {
  render: () => {
    const busy = new State<ShellBusy | undefined>({
      label: 'Loading 3 files…',
      progress: 0.4,
      cancel: () => busy.set(undefined),
    });
    return story(() => <DemoShell options={{ busy, setup: (s) => s.show('home') }} />);
  },
};

/** Compact (a narrow shell): the views in a bottom sheet, its bar of the `compactBar` views and More. */
export const Compact = shellStory({ width: '390px', height: '760px', setup: (s) => s.show('layers') });

/** Compact, a window shown: it opens as a card in the sheet. */
export const CompactWindow = shellStory({
  width: '390px',
  height: '760px',
  setup: (s) => s.show('details', { key: 'a' }),
});

/** Embedded (a viewer inside a page): the rail starts hidden; windows open smaller. */
export const Embedded = shellStory({ mode: 'embedded', width: '760px', height: '520px' });

/** A short shell: the rail's buttons step down, and what does not fit goes under More. */
export const ShortRail = shellStory({ height: '240px', setup: (s) => s.show('home') });

/** The right group widened (its splitter dragged, or arrows on it): Layers at 480px. */
export const ResizedPanel = shellStory({
  setup: (s) => {
    s.show('layers');
    s.resize(RightGroup, 480);
  },
});

/**
 * Hidden views (a rail button's menu: Hide): Inspect and Loading are out of the rail, their places kept. Their menu
 * (right-click a rail button, or Shift+F10) lists Show Inspect and Show Loading, which bring them back where they were;
 * Reset layout brings back every one. The story opens Home's menu.
 */
export const HiddenViews: StoryObj = {
  ...shellStory({
    setup: (s) => {
      s.show('home');
      s.hide('inspect');
      s.hide('loading');
    },
  }),
  play: async ({ canvasElement }) => {
    await new Promise((r) => setTimeout(r, 50));
    const button = canvasElement.querySelector<HTMLElement>('[data-rail-button][aria-label="Home"]');
    if (!button) return;
    const r = button.getBoundingClientRect();
    button.dispatchEvent(
      new MouseEvent('contextmenu', { bubbles: true, cancelable: true, clientX: r.left, clientY: r.bottom }),
    );
  },
};

/** The shell in the light, dark and high-contrast themes (`data-msui-theme` on a wrapper), side by side. */
export const Themes: StoryObj = {
  render: () =>
    story(() => (
      <div style="display: grid; grid-template-columns: repeat(3, minmax(560px, 1fr)); overflow-x: auto;">
        {['light', 'dark', 'contrast'].map((theme) => (
          <div data-msui-theme={theme} style="background: var(--msui-bg); color: var(--msui-text);">
            <DemoShell
              options={{
                mode: 'desktop',
                height: '100vh',
                setup: (s) => s.show('layers'),
              }}
            />
          </div>
        ))}
      </div>
    )),
};

/**
 * Files dropped anywhere on the shell (`onDropFiles`): drag files from your desktop over it. While a drag carries files,
 * one overlay says what a drop does (`dropLabel`); a drag of anything else (text, a link) shows nothing. Dropped, the
 * files go to the app (here, the stage fills and the status names them).
 */
export const FileDrop = shellStory({
  dropLabel: 'Drop files to open them',
  setup: (s) => s.show('home'),
});

/** The drop overlay as it shows while files are dragged over the shell (a drag started by the story). */
export const FileDropOverlay: StoryObj = {
  ...shellStory({
    dropLabel: 'Drop structures, sessions or maps to open them',
    dropFormats: 'PDB, mmCIF, BinaryCIF, MOL2, SDF, CCP4, sessions (.molx)',
  }),
  play: async ({ canvasElement }) => {
    await new Promise((r) => setTimeout(r, 50));
    const shell = canvasElement.querySelector('[role="application"]');
    if (!shell || typeof DataTransfer !== 'function') return;
    const data = new DataTransfer();
    data.items.add(new File(['data'], 'example.cif'));
    shell.dispatchEvent(new DragEvent('dragenter', { bubbles: true, cancelable: true, dataTransfer: data }));
  },
};

/**
 * Expanded (`view.expand`): the main area alone, the rail, the docked column and the windows away (paused); a card or
 * a page still opens over it. Escape, or showing a docked view, brings the rest back.
 */
export const Expanded = shellStory({
  setup: (s) => {
    s.show('home');
    s.show('details', { key: 'a', subject: 'Item A' });
    s.expand();
  },
});

/**
 * Expanded with `fill: 'page'`: a shell embedded in a page (here a box in a page of text) fills the page while expanded;
 * the page itself is not changed. The button toggles it.
 */
export const ExpandedPage: StoryObj = {
  render: () =>
    story(() => {
      let shell: ShellModel | undefined;
      return (
        <div style="display: grid; gap: var(--msui-space-3); padding: var(--msui-space-5); max-width: 720px;">
          <p style="margin: 0; color: var(--msui-text-muted);">
            A page with a viewer in it. Expand fills the page with it; Escape, or the button again, brings it back.
          </p>
          <div>
            <button type="button" onClick={() => shell?.toggleExpanded()}>
              Expand the viewer
            </button>
          </div>
          <DemoShell
            options={{
              mode: 'embedded',
              fill: 'page',
              width: '640px',
              height: '400px',
              setup: (s) => {
                shell = s;
              },
            }}
          />
        </div>
      );
    }),
};
