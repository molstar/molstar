/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { bind, domView, onActive, onCleanup } from '@molstar/prism-ui';
import type { ViewerExtensionUI, ViewerUIContext } from '../../spec.js';
import { StateTreeTab } from '../../views/state-tree/state-tree-tab.js';
import { LiveStateTreeModel } from './model.js';

// The State tree tab (parity.md: the classic State tree panel): the cells of the data tree, with updates and actions for
// the selected one, over the feature's `state` area. Its forms follow the cells while the tab shows; hidden, it holds
// no subscription to the state (prism AGENTS.md, hidden UI).

export function stateTreeUI(ctx: ViewerUIContext): ViewerExtensionUI {
  const state = ctx.api.state;
  if (!state) return {};
  const model = new LiveStateTreeModel(state);
  return {
    tabs: [
      {
        id: 'state-tree.tab.tree',
        title: 'State tree',
        icon: 'list-tree',
        order: 2,
        render: domView(() => {
          bind(state.tree.changed, (refs) => model.refresh(refs));
          bind(state.tree.roots, () => model.prune());
          // Shown again: whatever changed while it was hidden.
          onActive((active) => active && model.refresh());
          onCleanup(() => model.cancel());
          return <StateTreeTab model={model} />;
        }),
      },
    ],
  };
}
