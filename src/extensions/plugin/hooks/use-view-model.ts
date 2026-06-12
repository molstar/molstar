/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { MutableRefObject, useEffect, useRef } from 'react';
import { DefaultPluginSpec, PluginSpec } from '../../../mol-plugin/spec';
import { PluginViewModel } from '../view-model';

export function useCreatePluginViewModel<T extends PluginViewModel>(options?: {
    spec?: PluginSpec | ((defaultSpec: PluginSpec) => PluginSpec),
    model?: (spec: PluginSpec) => T
}): T {
    const model = useRef<T>();
    if (!model.current) {
        const spec = options?.spec
            ? (typeof options.spec === 'function' ? options.spec(DefaultPluginSpec()) : options.spec)
            : DefaultPluginSpec();
        model.current = options?.model ? options.model(spec) : new PluginViewModel({ spec }) as T;
    }
    return model.current!;
}

export function usePluginViewModel(view: PluginViewModel, parent: MutableRefObject<HTMLElement | null>) {
    useEffect(() => {
        if (!parent.current) return;
        view.mount(parent.current);
        return () => view.unmount();
    }, [view]);
}