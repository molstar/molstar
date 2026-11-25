/** @vitest-environment jsdom */
/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { State } from '@molstar/core/reactive';
import { describe, expect, it, vi } from 'vitest';
import { flush } from '../../dom/index.js';
import type { ShellBusy } from '../shell-view.js';
import { mountShell, probeView } from './fixture.js';

// Regressions of the shell's review, in the view: the compact sheet leaves the panes' bodies where they are, the busy
// overlay's Cancel follows the current task, a view registered anew is a new rail button.

describe('the shell view’s review regressions', () => {
  it('compact: a stage view’s body stays in its pane; the sheet shows only views docked beside the main area', () => {
    const stage = probeView('stage', { home: 'stage', canvas: true });
    const layers = probeView('layers', { compactBar: true });
    const m = mountShell([stage, layers], { mode: 'compact' });
    m.shell.show('stage');
    flush();
    expect(m.q('[data-pane] [data-probe="stage"]')).not.toBeNull();
    expect(m.q('[data-placement="sheet"] [data-probe="stage"]')).toBeNull();
    expect(m.probes.stage.mounts).toBe(1);
    m.shell.show('layers');
    flush();
    expect(m.q('[data-placement="sheet"] [data-probe="layers"]')).not.toBeNull();
    expect(m.q('[data-pane] [data-probe="stage"]')).not.toBeNull();
    m.dispose();
  });

  it('the busy overlay’s Cancel shows while the current task has one, and cancels that one', () => {
    const busy = new State<ShellBusy | undefined>(undefined);
    const m = mountShell([probeView('home')], {}, { busy });
    const cancel = () => m.q<HTMLButtonElement>('[role="alertdialog"] button');
    busy.set({ label: 'First' });
    flush();
    expect(cancel()).toBeNull();
    const second = vi.fn();
    busy.set({ label: 'Second', cancel: second });
    flush();
    expect(m.q('[role="alertdialog"]')!.textContent).toContain('Second');
    cancel()!.click();
    expect(second).toHaveBeenCalledTimes(1);
    busy.set({ label: 'Third' });
    flush();
    expect(cancel()).toBeNull();
    m.dispose();
  });

  it('a view registered anew under its id is a new rail button (its new title and icon)', () => {
    const m = mountShell([probeView('home')]);
    m.shell.show('home');
    flush();
    expect(m.railButton('Home')).not.toBeNull();
    m.shell.registerViews([{ ...probeView('home').def, title: 'Start', icon: 'house' }]);
    flush();
    expect(m.railButton('Home')).toBeNull();
    expect(m.railButton('Start')).not.toBeNull();
    m.dispose();
  });
});
