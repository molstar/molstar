/**
 * Copyright (c) 2025-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { MolstarLoadingExtension } from '@molstar/mvs/load';
import { UpdateTarget } from '@molstar/mvs/load-generic';

export const IsHiddenCustomStateExtension: MolstarLoadingExtension<{}> = {
  id: 'ww-pdb/is-hidden-custom-state',
  description: 'Allow updating initial visibility of nodes',
  local: true,
  createExtensionContext: () => ({}),
  action: (updateTarget, node) => {
    if (!node.custom || !node.custom?.is_hidden) return;
    UpdateTarget.updateState(updateTarget, { isHidden: true });
  },
};
