/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The skin's text, inlined by the bundler (Vite's `?inline`): the element adopts it as constructed stylesheets.
declare module "*.css?inline" {
  const css: string;
  export default css;
}

// The skin as a side-effect import (the browser tests load it to measure the legend's grid).
declare module "*.css";
