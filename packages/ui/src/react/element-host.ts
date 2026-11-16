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
  const properties = props.properties;
  useLayoutEffect(() => {
    const el = ref.current as (HTMLElement & Record<string, unknown>) | null;
    if (!el || !properties) return;
    for (const key of Object.keys(properties)) el[key] = properties[key];
  }, [properties]);
  return createElement(props.tag, { ref, className: props.className, style: props.style });
}
