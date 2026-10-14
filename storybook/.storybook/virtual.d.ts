/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// story-previews.ts
declare module 'virtual:story-package-previews' {
  import type { Decorator } from '@storybook/web-components-vite';

  const decorators: Decorator[];
  export default decorators;
}
