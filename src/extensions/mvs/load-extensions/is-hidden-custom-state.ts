/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { MolstarLoadingExtension } from '../load';

export const IsHiddenCustomStateExtension: MolstarLoadingExtension<{}> = {
    id: 'ww-pdb/is-hidden-custom-state',
    description: 'Allow updating initial visibility of nodes',
    createExtensionContext: () => ({}),
    action: (updateTarget, node) => {
        if (!node.custom || !node.custom?.is_hidden) return;
        updateTarget.update.to(updateTarget.selector).updateState({ isHidden: true });
    },
};