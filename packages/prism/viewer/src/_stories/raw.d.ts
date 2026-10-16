/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Storybook's Vite: `?raw` imports give the file's text (the stories read the repository's example data this way).
declare module '*?raw' {
  const text: string;
  // biome-ignore lint/style/noDefaultExport: the shape of a `?raw` module
  export default text;
}
