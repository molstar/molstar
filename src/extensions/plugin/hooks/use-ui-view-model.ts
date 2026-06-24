/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { useRef } from 'react';
import { DefaultPluginUISpec, PluginUISpec } from '../../../mol-plugin-ui/spec';
import { PluginUIViewModel } from '../ui-view-model';

export function useCreatePluginUIViewModel<T extends PluginUIViewModel = PluginUIViewModel>(options?: {
    spec?: PluginUISpec | ((defaultSpec: PluginUISpec) => PluginUISpec),
    model?: (spec: PluginUISpec) => T
}): T {
    const model = useRef<T>();
    if (!model.current) {
        const spec = options?.spec
            ? (typeof options.spec === 'function' ? options.spec(DefaultPluginUISpec()) : options.spec)
            : DefaultPluginUISpec();
        model.current = options?.model ? options.model(spec) : new PluginUIViewModel({ spec }) as T;
    }
    return model.current;
}