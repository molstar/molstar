/** @vitest-environment jsdom */
/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { afterEach, beforeAll, describe, expect, it } from 'vitest';
import { render } from '../../../dom/render.js';
import { flushFrame } from '../../../dom/scope.js';
import { domView } from '../../../dom/view.js';
import { ShellModel, type ShellViewContext } from '../../../shell/shell.js';
import { lostReleases, stubPointerCapture } from '../../../testing/gestures.js';
import { DockView } from '../dock-view.js';

// The dock's splitter is a `drag` (dom/scoped.ts): pressed, moved and left without a release again and again, it leaves
// no listener behind, and the next drag resizes the dock.

beforeAll(stubPointerCapture);
afterEach(() => {
  document.body.innerHTML = '';
});

describe('the dock splitter', () => {
  it('no listener left by lost releases, and the next drag resizes the dock', () => {
    const shell = new ShellModel();
    shell.registerTabs([
      { id: 'tree', title: 'Tree', icon: 'house', dock: 'left', render: domView<ShellViewContext>(() => <i>tree</i>) },
    ]);
    shell.show('tree');
    const dock = shell.dock('left');
    const box = document.createElement('div');
    document.body.append(box);
    render(() => <DockView shell={shell} dock={dock} />, box);
    const splitter = () => box.querySelector<HTMLElement>('[role=separator]')!;
    expect(splitter()).not.toBeNull();
    // jsdom lays nothing out: the handle's box, for whether the release was over it.
    splitter().getBoundingClientRect = () => ({ left: 0, right: 4, top: 0, bottom: 400 }) as DOMRect;
    let start = dock.state.value.size;
    const r = lostReleases({
      target: splitter,
      from: { x: 2, y: 100 },
      to: { x: 42, y: 100 },
      reset: () => {
        flushFrame();
        start = dock.state.value.size;
      },
      worked: () => {
        flushFrame();
        // Released off the handle: the light goes.
        return dock.state.value.size === start + 40 && !splitter().hasAttribute('data-active');
      },
    });
    expect(r.grown).toBe(0);
    expect(r.worked).toBe(true);
  });
});
