/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Resolve optional Node rendering dependencies when they are actually needed.
 */
import { createRequire } from 'node:module';
import type { ExternalModules } from './screenshot.js';

const require = createRequire(import.meta.url);

export function loadNativeModule(name: 'gl'): ExternalModules['gl'];
export function loadNativeModule(name: 'canvas'): unknown;
export function loadNativeModule(name: 'gl' | 'canvas'): unknown {
    try {
        return require(name);
    } catch (cause) {
        const setup = name === 'canvas' ? 'pnpm native:install -- --canvas' : 'pnpm native:install';
        throw new Error(`Optional native module '${name}' is unavailable. Install it in your application, or in this workspace run '${setup}' and launch the command with 'pnpm native:run -- <command>'. ${cause instanceof Error ? cause.message : String(cause)}`, { cause });
    }
}
