/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The browser's pieces for the SVG export (seqspec:viewer#the-render-backend): the `deflate` port over the platform's
// `CompressionStream`, and `view.exportSvg` registered with it. The `/svg` entry takes both from the host.

import type { ActionEntry } from "../actions.js";
import { canvasEmBox } from "../canvas/shared.js";
import type { SeqViewerSpec } from "../index.js";
import { type ExportSvgOptions, exportSvgAction } from "../svg/index.js";
import { downloadSvg } from "./download.js";

/** The canvas's em boxes (`canvasMeasurer` answers `metrics` with them), here for the export's hosts. */
export { canvasEmBox };

/** `SvgPorts.deflate` for a page: a zlib stream out of `CompressionStream("deflate")`. */
export async function compressionDeflate(bytes: Uint8Array): Promise<Uint8Array> {
  const stream = new CompressionStream("deflate");
  const writer = stream.writable.getWriter();
  void writer.write(bytes as Uint8Array<ArrayBuffer>).then(() => writer.close());
  return new Uint8Array(await new Response(stream.readable).arrayBuffer());
}

/**
 * `view.exportSvg` with the browser's `deflate`: `{ actions: [browserExportAction({ offer })] }` in the registry. The em
 * boxes are the context's measurer's own (`canvasMeasurer` answers `metrics`, refreshed on a font load), unless `render` gives some.
 */
export function browserExportAction(
  o: Omit<ExportSvgOptions, "ports"> & { ports?: ExportSvgOptions["ports"] },
): ActionEntry {
  return exportSvgAction({
    ...o,
    ports: { deflate: compressionDeflate, ...o.ports },
    render: o.render,
  });
}

const EXPORT = "view.exportSvg";
/**
 * The spec of a context a UI makes: the UI's Export offers a download (`downloadSvg`) through the browser's
 * `view.exportSvg`, unless the spec registers its own. Both UIs' wrappers (`<seq-viewer>`, `SeqViewer`) use it.
 */
export function withDownloadExport(spec: SeqViewerSpec): SeqViewerSpec {
  return spec.registry?.some((e) => e.actions?.some((a) => a.name === EXPORT))
    ? spec
    : { ...spec, registry: [...(spec.registry ?? []), { actions: [browserExportAction({ offer: downloadSvg })] }] };
}
