/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type CSSProperties, useRef } from 'react';
import type { PluginViewModel } from './view-model';
import { usePluginViewModel } from './hooks/use-view-model';

export function PluginCanvas({ model, style, className }: { model: PluginViewModel, style?: CSSProperties, className?: string }) {
    const root = useRef<HTMLDivElement>(null);
    usePluginViewModel(model, root);
    return <div ref={root} style={style} className={className} />;
}