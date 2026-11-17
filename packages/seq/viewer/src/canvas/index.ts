/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The Canvas2D backend (_spec/README.md "Canvas2D backend"): draw lists on a base canvas and an emphasis overlay, and the
// text metrics port's browser adapter. Named re-exports only.

export { CanvasBackend, type CanvasBackendOptions, type CanvasLike, type Context2D } from './backend.js';
export { cssColorToHex, readCssTokens, ThemeAttribute, themeOf } from './css-tokens.js';
export { cssFont, FontCache } from './font.js';
export { browserTextMetrics, canvasTextMetrics, type MeasuringContext } from './metrics.js';
