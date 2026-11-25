/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */
import { type CSSProperties, createElement, useLayoutEffect, useRef } from 'react';

export interface UiElementHostProps {
  /** The element's tag. */
  tag: string;
  /** Set as properties on the element (models, callbacks): React 18 would turn them into attributes. */
  properties?: Readonly<Record<string, unknown>>;
  className?: string;
  style?: CSSProperties;
}

/**
 * Hosts one of the UI's elements (a `UiElement`, a `ViewElement`) in a React tree, setting its properties through a
 * ref (react §4). React 19 sets properties on custom elements itself; this works from React 18 on and with
 * `preact/compat`.
 */
export function UiElementHost(props: UiElementHostProps) {
  const ref = useRef<HTMLElement>(null);
  // The keys set last time: one no longer given is set back to `undefined`, so the element drops it (a callback with a
  // stale closure is not kept).
  const set = useRef<readonly string[]>([]);
  const properties = props.properties;
  useLayoutEffect(() => {
    const el = ref.current as (HTMLElement & Record<string, unknown>) | null;
    if (!el) return;
    const keys = properties ? Object.keys(properties) : [];
    for (const key of set.current) if (!properties || !Object.hasOwn(properties, key)) el[key] = undefined;
    // Only what changed is written: an inline object (new each render) sets nothing that is already there.
    for (const key of keys) if (el[key] !== properties![key]) el[key] = properties![key];
    set.current = keys;
  }, [properties]);
  return createElement(props.tag, { ref, className: props.className, style: props.style });
}
