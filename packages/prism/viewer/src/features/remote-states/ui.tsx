/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { domView, onActive } from '@molstar/prism-ui';
import type { ViewerExtensionUI, ViewerUIContext } from '../../spec.js';
import { RemoteStatesList } from '../../views/home/home-tab.js';
import { RemoteStatesPanel } from '../../views/state/remote-states-panel.js';
import { createRemoteStatesModel } from './model.js';

// The Remote states card (parity.md §1): upload the current state, the states on the server (a click opens one) and the
// server's address, over the `session.remote.*` operations; and the same list on Home (a section of `home.import`, as
// the classic Home had it). One model for the viewer's life and both views, so the address the person typed stays; the
// list is read each time either shows (a folded section or a hidden Home reads nothing: `onActive`).

export function remoteStatesUI(ctx: ViewerUIContext): ViewerExtensionUI {
  if (!ctx.api.session.remoteList) return {};
  const model = createRemoteStatesModel({
    api: ctx.api,
    notify: (text, kind) => ctx.managers.notifications.post(text, kind),
  });
  return {
    sections: [
      {
        id: 'remote-states.section.home',
        slot: 'home.import',
        title: 'Remote states',
        icon: 'cloud',
        // After the loaders' and extensions' import forms.
        order: 10,
        keywords: ['remote states', 'state server', 'sessions'],
        render: domView(() => {
          onActive((active) => active && model.refresh());
          return <RemoteStatesList model={model} />;
        }),
      },
    ],
    overlays: [
      {
        id: 'remote-states.overlay.states',
        title: 'Remote states',
        render: domView(() => {
          onActive((active) => active && model.refresh());
          return <RemoteStatesPanel model={model} />;
        }),
      },
    ],
    toolbar: [
      {
        id: 'remote-states.toolbar.states',
        icon: 'cloud',
        title: 'Remote states',
        group: 'export',
        overlay: 'remote-states.overlay.states',
      },
    ],
  };
}
