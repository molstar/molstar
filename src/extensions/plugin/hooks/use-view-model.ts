/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { MutableRefObject, useEffect, useRef } from 'react';
import { DefaultPluginSpec, PluginSpec } from '../../../mol-plugin/spec';
import { PluginViewModel } from '../view-model';

export function useCreatePluginViewModel(options?: { spec?: PluginSpec | ((defaultSpec?: PluginSpec) => PluginSpec) }): PluginViewModel {
    const model = useRef<PluginViewModel>();
    if (!model.current) {
        model.current = new PluginViewModel({
            spec: options?.spec
                ? (typeof options.spec === 'function' ? options.spec(DefaultPluginSpec()) : options.spec)
                : options?.spec as PluginSpec
        });
    }
    return model.current;
}

export function usePluginViewModel(view: PluginViewModel, parent: MutableRefObject<HTMLElement | null>) {
    useEffect(() => {
        if (!parent.current) return;
        view.mount(parent.current);
        return () => view.unmount();
    }, [view]);
}