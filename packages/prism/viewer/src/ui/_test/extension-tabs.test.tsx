/** @vitest-environment jsdom */
/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { domView } from '@molstar/prism-ui/dom/view';
import { observableState } from '@molstar/prism-ui/reactive/observable-state';
import { ShellModel } from '@molstar/prism-ui/shell/shell';
import { describe, expect, it } from 'vitest';
import type { DockId } from '@molstar/prism-ui/managers/view-state';
import type { ViewerUIContext } from '../../spec.js';
import { extensionUIPort } from '../extensions.js';

// An extension's tab with `visible` (viewer-spec.md §3): on the rail while it holds (a kinemage's panel while one is
// loaded), off it otherwise; not one of the person's hidden tabs.

const docks: DockId[] = ['left', 'right', 'top', 'bottom'];
const onRail = (shell: ShellModel, id: string) => docks.some((d) => shell.dock(d).tabs.value.some((t) => t.id === id));

describe('extension tabs', () => {
  it('puts a tab on the rail while its visible holds, and takes it off with its extension', () => {
    const shell = new ShellModel();
    const port = extensionUIPort(shell, {} as never, {} as ViewerUIContext);
    const loaded = observableState(false);
    const unregister = port.register('x', {
      tabs: [
        { id: 'x.tab.always', title: 'Always', icon: 'box', render: domView(() => <p />) },
        { id: 'x.tab.loaded', title: 'Loaded', icon: 'box', visible: loaded.asReadonly, render: domView(() => <p />) },
      ],
    });
    expect(onRail(shell, 'x.tab.always')).toBe(true);
    expect(onRail(shell, 'x.tab.loaded')).toBe(false);
    loaded.set(true);
    expect(onRail(shell, 'x.tab.loaded')).toBe(true);
    loaded.set(false);
    expect(onRail(shell, 'x.tab.loaded')).toBe(false);
    expect(shell.hiddenTabs.value).toEqual([]);
    loaded.set(true);
    unregister();
    expect(onRail(shell, 'x.tab.loaded')).toBe(false);
    expect(onRail(shell, 'x.tab.always')).toBe(false);
  });
});
