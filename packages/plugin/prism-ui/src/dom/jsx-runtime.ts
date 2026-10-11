/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The automatic JSX runtime (`"jsxImportSource": "@molstar/plugin-prism-ui/dom"`): JSX creates real DOM nodes.
import type { Component } from './jsx-types.js';
import { createElement, Fragment } from './render.js';

export type { JSX } from './jsx-types.js';
export { Fragment };

export function jsx(type: string | Component<any>, props: Record<string, any> | null, _key?: unknown): Node {
  return createElement(type, props);
}

export const jsxs = jsx;
