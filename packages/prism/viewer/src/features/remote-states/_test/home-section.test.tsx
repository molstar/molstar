/** @vitest-environment jsdom */
/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ShellModel } from '@molstar/prism-ui/shell/shell';
import { afterEach, describe, expect, it, vi } from 'vitest';
import { headlessViewer, until } from '../../../_test/harness.js';
import { LazyRemoteStates } from '../../lazy/index.js';
import { RemoteStatesFeature } from '../index.js';
import { remoteStatesUI } from '../ui.js';

// The remote states on Home (parity.md §1, "the list on Home"): the `home.import` section lists the server's states as
// it shows, and a click opens one, over a stubbed `fetch` (no network).

const SERVER = 'https://state.test';
const LIST = [
  { id: 'a', name: 'First', description: '', timestamp: 1000, isSticky: false },
  { id: 'b', name: 'Second', description: '', timestamp: 2000, isSticky: false },
];

afterEach(() => {
  vi.unstubAllGlobals();
  document.body.innerHTML = '';
});

describe("the remote states' Home section", () => {
  it('lists the states as it shows, and a click opens one', async () => {
    const stub = vi.fn(async (input: string | URL | Request) => {
      const url = String(typeof input === 'object' && 'url' in input ? input.url : input);
      if (url === `${SERVER}/list`) return new Response(JSON.stringify(LIST));
      return new Response('not found', { status: 404 });
    });
    vi.stubGlobal('fetch', stub);
    const v = await headlessViewer({
      extensions: [RemoteStatesFeature],
      spec: { options: { data: { stateServer: SERVER }, ui: { showRemoteStates: true } } },
    });
    const ui = remoteStatesUI({ api: v.api.as({ kind: 'person' }), sources: v.sources, managers: v.managers });
    const section = ui.sections?.find((s) => s.id === 'remote-states.section.home');
    expect(section).toMatchObject({ slot: 'home.import', title: 'Remote states' });
    // The stub declares it, so the lazy feature's Home has it before it loads.
    expect(LazyRemoteStates.stub.parts?.map((p) => p.id)).toContain('remote-states.section.home');

    const host = document.createElement('div');
    document.body.append(host);
    section!.render.mount(host, { shell: new ShellModel() } as never, undefined as never);
    await until(() => host.textContent?.includes('Second') ?? false);
    expect(host.textContent).toContain('First');

    const row = [...host.querySelectorAll('button')].find((b) => b.textContent?.includes('First'))!;
    row.click();
    // Open reads the list and then the state (the stub has none: the failure is posted, the scene stays).
    await until(() => stub.mock.calls.some((c) => String(c[0]) === `${SERVER}/get/a`));
  });
});
