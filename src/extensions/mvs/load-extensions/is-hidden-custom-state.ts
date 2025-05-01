/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { MolstarLoadingExtension } from '../load';


export const IsHiddenCustomStateExtension: MolstarLoadingExtension<{}> = {
    id: 'mvs/is-hidden-custom-state',
    description: 'Allow updating initial visibility of nodes',
    createExtensionContext: () => ({}),
    action: (updateTarget, node) => {
        if (!node.custom || !node.custom?.molstar_is_hidden) return;
        updateTarget.builder.updateState({ isHidden: true });
    },
};
