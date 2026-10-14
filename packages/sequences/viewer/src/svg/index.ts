/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// /svg: the SVG backend and `renderSvg` (seqspec:viewer#the-render-backend). Platform-free: no DOM and no Node.
export { type ExportSvgOptions, type ExportSvgTarget, exportSvgAction } from "./action.js";
export {
  SvgBackend,
  type SvgFont,
  type SvgImage,
  type SvgOptions,
  type SvgPorts,
  type SvgShared,
  type SvgStats,
} from "./backend.js";
export {
  type EmBox,
  type FontMetrics,
  hasMetrics,
  type MetricsMeasurer,
  metricsMeasurer,
} from "./metrics.js";
export { DefaultFontMetrics } from "./metrics-tables.js";
export { encodePng, MAX_PNG_BYTES } from "./png.js";
export { type RenderSvgOptions, renderSvg } from "./render.js";
