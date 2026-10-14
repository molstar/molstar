/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { DocsContainer } from '@storybook/addon-docs/blocks';
import { type ComponentProps, createElement, useLayoutEffect, useState } from 'react';
import { GLOBALS_UPDATED } from 'storybook/internal/core-events';
import { themes } from 'storybook/theming';

type Props = ComponentProps<typeof DocsContainer>;
type Globals = Record<string, unknown>;

/**
 * The docs container: a docs page follows the toolbar's theme as a story does (`data-msp-theme` on the document root,
 * which the page's own components and Prism's tokens read; Storybook's dark theme for the page itself on dark and blue).
 */
export function ThemedDocsContainer(props: Props) {
  const docs = props.context as unknown as {
    channel: { on(event: string, fn: (e: { globals: Globals }) => void): void; off(event: string, fn: never): void };
    store: { userGlobals: { get(): Globals } };
  };
  const [theme, setTheme] = useState<unknown>(() => docs.store.userGlobals.get().theme ?? 'light');
  useLayoutEffect(() => {
    document.documentElement.dataset.mspTheme = String(theme ?? 'light');
  }, [theme]);
  useLayoutEffect(() => {
    const update = (e: { globals: Globals }) => setTheme(e.globals.theme ?? 'light');
    docs.channel.on(GLOBALS_UPDATED, update);
    return () => docs.channel.off(GLOBALS_UPDATED, update as never);
  }, [docs]);
  const dark = theme === 'dark' || theme === 'blue';
  return createElement(DocsContainer, { ...props, theme: dark ? themes.dark : themes.light });
}
