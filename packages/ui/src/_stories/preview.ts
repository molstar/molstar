/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The UI group's preview module (storybook/packages.ts): the toolbar's theme on the document root as
// `data-msui-theme`. Styles come with the stories' roots (`UiElement` adopts the registered styles).
import type { Decorator } from '@storybook/web-components-vite';

export const decorators: Decorator[] = [
  (story, context) => {
    document.documentElement.dataset.msuiTheme = context.globals.theme ?? 'light';
    return story();
  },
];
