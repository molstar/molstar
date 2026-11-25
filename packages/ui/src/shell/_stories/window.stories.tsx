/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { Note, page, story } from '../../_stories/support.js';
import { findView, type Rect, type ViewRef } from '../layout.js';
import type { ShellModel } from '../shell-model.js';
import { type DemoOptions, DemoShell } from './demo.js';

/**
 * Windows (shell §2): views floating over the main area, never over a docked group. The header drags one, the right
 * and bottom edges and the corner resize it, a press or focus raises it; ✕ or Escape closes it, a press outside does
 * not. Each is kept inside the main area. Its ⋯ menu docks it, moves it into another window, moves and resizes it from
 * the keyboard, and resets its position.
 */
const meta: Meta = { title: 'Shell/Window', parameters: { layout: 'fullscreen' } };
export default meta;

/** The window group holding `ref`. */
function windowOf(s: ShellModel, ref: ViewRef): string | undefined {
  const place = findView(s.layout.value, ref);
  return place?.window ? place.group.id : undefined;
}

/** Puts the window of `ref` at `rect` (relative to the main area). */
function place(s: ShellModel, ref: ViewRef, rect: Rect) {
  const id = windowOf(s, ref);
  if (id) s.setWindowRect(id, rect);
}

const shellStory = (options: DemoOptions = {}): StoryObj => ({
  render: () => story(() => <DemoShell options={options} />),
});

/** One view popped out (Inspect, from its ⋯ menu's Float): a window of one view shows only its title. */
export const One = shellStory({
  setup: (s) => {
    s.float('inspect');
    place(s, 'inspect', { x: 80, y: 60, width: 340, height: 360 });
  },
});

/** Several windows stacked: the one pressed last on top (press another to raise it). */
export const Stacked = shellStory({
  setup: (s) => {
    s.float('inspect');
    s.float('layers');
    s.float('help');
    place(s, 'inspect', { x: 40, y: 40, width: 320, height: 300 });
    place(s, 'layers', { x: 120, y: 100, width: 300, height: 300 });
    place(s, 'help', { x: 200, y: 160, width: 300, height: 220 });
  },
});

/** A window of two views: their tabs in its header (a view's ⋯ menu: Move to the … window). */
export const Tabs = shellStory({
  setup: (s) => {
    s.float('inspect');
    const id = windowOf(s, 'inspect');
    if (id) s.move('help', { kind: 'group', id });
    s.select('inspect');
    place(s, 'inspect', { x: 80, y: 60, width: 360, height: 340 });
  },
});

/**
 * A small main area (480 × 360): a window put past its bottom-right corner is clamped back into it, its header
 * reachable; a main area smaller than a window shrinks it to fit. Resize the story's frame to see it follow.
 */
export const Clamped = shellStory({
  mode: 'desktop',
  width: '520px',
  height: '360px',
  setup: (s) => {
    s.float('inspect');
    place(s, 'inspect', { x: 360, y: 240, width: 300, height: 220 });
  },
});

/** Two instances of one view (`instances: 'many'`, opened by a key with a subject): two details views side by side. */
export const Instances = shellStory({
  setup: (s) => {
    s.show('details', { key: 'a', subject: 'Item A' });
    s.show('details', { key: 'b', subject: 'Item B' });
    place(s, 'details#a', { x: 24, y: 24, width: 340, height: 260 });
    place(s, 'details#b', { x: 380, y: 24, width: 340, height: 260 });
  },
});

/**
 * Move and resize from the keyboard: the window's ⋯ menu → Move with the keyboard (or Resize with the keyboard). The
 * window is outlined while it adjusts; arrows move or resize it by 4px, Shift+arrows by 16, Enter or Escape ends. The
 * mode is the window's own (not in the model), so this story opens the window and leaves the menu to you; Reset
 * position centers it again.
 */
export const KeyboardMoveResize: StoryObj = {
  name: 'Keyboard move and resize',
  render: () =>
    story(() => (
      <div style="display: grid; grid-template-rows: auto 1fr; height: 100vh;">
        <div style={page}>
          <Note>
            Open the window's ⋯ menu and pick Move with the keyboard: arrows move it (Shift: by 16px), Enter or Escape
            ends. Resize with the keyboard does the same to its size. F6 cycles the main area, the docked group and each
            window.
          </Note>
        </div>
        <DemoShell
          options={{
            height: '100%',
            setup: (s) => {
              s.float('inspect');
              place(s, 'inspect', { x: 120, y: 60, width: 340, height: 300 });
            },
          }}
        />
      </div>
    )),
};
