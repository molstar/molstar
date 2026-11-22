/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Vite's `?raw` import: a module's source text (the Builder stories show the code they run).
declare module '*?raw' {
  const source: string;
  // biome-ignore lint/style/noDefaultExport: the shape of Vite's `?raw` modules
  export default source;
}
