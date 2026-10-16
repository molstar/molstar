/** @vitest-environment jsdom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { cbs, crambin, expectOk } from '../../../_test/harness.js';
import { DefaultPrismViewerSpec } from '../../../default-spec.js';
import { createDefaultPlugin } from '../../../sources/create-plugin.js';
import { mountShell } from '../../../ui/shell.js';
import { createViewerCore, type ViewerCore } from '../../../viewer.js';

// A load the person starts from Home shows the Structure tab, so what was loaded is in view (David, 2026-10-16).

class NoResize {
  observe() {}
  unobserve() {}
  disconnect() {}
}

const cores: ViewerCore[] = [];
beforeEach(() => vi.stubGlobal('ResizeObserver', NoResize));
afterEach(async () => {
  for (const c of cores.splice(0)) await c.dispose();
  document.body.innerHTML = '';
  vi.restoreAllMocks();
  vi.unstubAllGlobals();
});

async function mounted() {
  const c = await createViewerCore({
    spec: DefaultPrismViewerSpec(),
    createPlugin: createDefaultPlugin,
    modules: ['files', 'page'],
  });
  cores.push(c);
  vi.spyOn(c, 'mountCanvas').mockResolvedValue(false);
  const target = document.createElement('div');
  document.body.append(target);
  const m = await mountShell(c, target);
  const active = () => m.shell.dock('right').active.value?.id;
  return { c, m, active, person: c.api.as({ kind: 'person' }) };
}

describe('a load from Home', () => {
  it('shows the Structure tab once something came in', async () => {
    const { m, active, person } = await mounted();
    m.shell.show('home');
    expect(active()).toBe('home');
    expectOk(await person.load.files({ files: [cbs()] }));
    expect(active()).toBe('structure');
    m.dispose();
  });

  it('stays where the person went, and ignores loads not started by the person', async () => {
    const { c, m, active, person } = await mounted();
    m.shell.show('home');
    expectOk(await c.api.load.files({ files: [cbs()] }));
    expect(active()).toBe('home');
    m.shell.show('selection');
    expectOk(await person.load.files({ files: [crambin()] }));
    expect(active()).toBe('selection');
    m.dispose();
  });
});
