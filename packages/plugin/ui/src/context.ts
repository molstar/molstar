/**
 * Copyright (c) 2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */


import { PluginContext } from '@molstar/plugin/context';
import type { PluginUISpec } from '@molstar/plugin-ui/spec';
import type { StateTransformParameters } from '@molstar/plugin-ui/state/common';

export class PluginUIContext extends PluginContext {
    readonly customParamEditors = new Map<string, StateTransformParameters.Class>();
    readonly customUIState: Record<string, any> = {};

    private initCustomParamEditors() {
        if (!this.spec.customParamEditors) return;

        for (const [t, e] of this.spec.customParamEditors) {
            this.customParamEditors.set(t.id, e);
        }
    }

    dispose(options?: { doNotForceWebGLContextLoss?: boolean }) {
        super.dispose(options);
        this.layout.dispose();
    }

    constructor(public spec: PluginUISpec) {
        super(spec);

        this.initCustomParamEditors();
    }
}