/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Vite's `?raw` imports (the bundled fixtures of fixture-fs.ts).
declare module "*?raw" {
  const text: string;
  export default text;
}
