/** @vitest-environment jsdom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Scope } from '@molstar/prism-ui/dom/scope';
import { ShellModel } from '@molstar/prism-ui/shell/shell';
import { afterEach, describe, expect, it, vi } from 'vitest';
import type { ViewerUIContext } from '../../spec.js';
import { extensionUIPort } from '../extensions.js';

// A lazy extension's placeholder section (viewer-spec.md §3): showing it loads the extension; a load that fails says
// so and waits for Retry, also when the host puts the placeholder back (no loop of loads).

afterEach(() => {
  document.body.innerHTML = '';
});

const flush = () => new Promise((r) => setTimeout(r, 0));

function show(shell: ShellModel) {
  const host = document.createElement('div');
  document.body.append(host);
  const section = shell.sections('export').value[0];
  section.render.mount(host, {} as never, new Scope());
  return host;
}

describe('lazy placeholders', () => {
  it('loads when shown; a failed load shows Retry, and a placeholder put back does not load again', async () => {
    const shell = new ShellModel();
    const port = extensionUIPort(shell, {} as never, {} as ViewerUIContext);
    const stub = { parts: [{ id: 'x.section.a', kind: 'section' as const, slot: 'export' as const, title: 'A' }] };
    const load = vi.fn(async () => false);

    let remove = port.registerStub('x', stub, load);
    let host = show(shell);
    expect(host.textContent).toContain('Loading…');
    await flush();
    expect(load).toHaveBeenCalledTimes(1);
    expect(host.textContent).toContain('could not be loaded');

    // The host puts the stub back after the failure: it waits.
    remove();
    remove = port.registerStub('x', stub, load);
    host = show(shell);
    await flush();
    expect(load).toHaveBeenCalledTimes(1);
    expect(host.textContent).toContain('could not be loaded');

    // Retry loads again; this time it works.
    load.mockImplementation(async () => true);
    host.querySelector('button')!.click();
    await flush();
    expect(load).toHaveBeenCalledTimes(2);
    remove();
  });

  it('loads when it is shown, not while its section is folded or its panel hidden', async () => {
    const shell = new ShellModel();
    const port = extensionUIPort(shell, {} as never, {} as ViewerUIContext);
    const stub = { parts: [{ id: 'x.section.a', kind: 'section' as const, slot: 'export' as const, title: 'A' }] };
    const load = vi.fn(async () => true);
    const remove = port.registerStub('x', stub, load);
    const host = document.createElement('div');
    document.body.append(host);
    // Rendered in a paused scope, as a folded section's body is.
    const scope = new Scope();
    scope.setActive(false);
    shell.sections('export').value[0].render.mount(host, {} as never, scope);
    await flush();
    expect(load).not.toHaveBeenCalled();
    scope.setActive(true);
    await flush();
    expect(load).toHaveBeenCalledTimes(1);
    // Paused and shown again: once.
    scope.setActive(false);
    scope.setActive(true);
    await flush();
    expect(load).toHaveBeenCalledTimes(1);
    remove();
  });
});
