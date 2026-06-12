/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { useRef } from 'react';
import { DefaultPluginUISpec, PluginUISpec } from '../../../mol-plugin-ui/spec';
import { PluginUIViewModel } from '../ui-view-model';

export function useCreatePluginUIViewModel(options?: { spec?: PluginUISpec | ((defaultSpec: PluginUISpec) => PluginUISpec) }): PluginUIViewModel {
    const model = useRef<PluginUIViewModel>();
    if (!model.current) {
        model.current = new PluginUIViewModel({
            spec: options?.spec
                ? (typeof options.spec === 'function' ? options.spec(DefaultPluginUISpec()) : options.spec)
                : options?.spec as PluginUISpec
        });
    }
    return model.current;
}