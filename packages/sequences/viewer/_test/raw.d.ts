/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The stories' support files import sources as text (`?raw`, Vite): story-search.test.ts reads one of them.
declare module '*?raw' {
  const source: string;
  // biome-ignore lint/style/noDefaultExport: an ambient module declaration of a default-exporting import.
  export default source;
}
