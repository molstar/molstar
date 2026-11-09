/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

// The types of vitest.mjs (plain JavaScript: the configs load it natively).

import type { InlineConfig } from 'vitest/node';

type Sources = {
  resolve: { conditions: string[] };
  ssr: { resolve: { conditions: string[] } };
};

export declare const GUARD_FLAGS: string[];
export declare const BROWSER_DIRS: string[];
export declare const SOURCES: Sources;
export declare const heapTests: () => Sources & { test: InlineConfig };
export declare const perfTests: (dir: string) => Sources & { test: InlineConfig };
