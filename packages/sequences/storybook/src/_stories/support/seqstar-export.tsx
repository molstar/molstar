/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The export story's cards (Seq*/Export): a live viewer, buttons that run the UI's own export action, and the file it
// made shown beside the canvas. The action is `browserExportAction` (`@molstar/seq-viewer/dom`), registered on the context as
// the Export part finds it and run as the part runs it: `ctx.actions.run("view.exportSvg", { width })`, the width the
// live section's. A PNG is that SVG drawn by the page itself (an image, a canvas): the viewer has no PNG action, the
// UI's Export offers SVG only; `seqv-render --png` rasterizes the same file with Chrome.

import type { SeqViewSpecDoc } from "@molstar/seq-core/format";
import type { SeqViewerRegistryEntry } from "@molstar/seq-viewer";
import { browserExportAction, downloadSvg } from "@molstar/seq-viewer/dom";
import { SeqViewerLayout, SeqViewerProvider } from "@molstar/seq-viewer/react";
import { type CSSProperties, type ReactNode, useEffect, useMemo, useRef, useState } from "react";
import { PrismPage } from "../support/seqstar-prism.js";
import { Choices, TextAction } from "../ui/index.js";
import f from "./seqstar-format.module.css";
import { PrismCanvas, prism } from "./seqstar-prism.js";
import t from "./seqstar-tools.module.css";
import { useCaseContext, useStoryUi } from "./seqstar-viewer.js";

/** A file the export made: the SVG as the action offered it, or the PNG the page drew from it. */
interface Made {
  readonly kind: "svg" | "png";
  readonly name: string;
  readonly url: string;
  readonly bytes: number;
  /** The export's time from the click to the file: the action (render and queue), then, for a PNG, the page's raster. */
  readonly ms: number;
  readonly detail: string;
  /** The file's size in px. */
  readonly size: string;
}

const kb = (n: number) => (n >= 1e6 ? `${(n / 1e6).toFixed(2)} MB` : `${(n / 1e3).toFixed(1)} KB`);
const ms = (x: number) => `${Math.round(x).toLocaleString("en-US")} ms`;

/** The SVG drawn at its own size on a canvas and written as a PNG (transparent where the SVG is). */
async function rasterize(svg: string): Promise<{ png: Blob; w: number; h: number }> {
  const head = /<svg[^>]* width="([\d.]+)" height="([\d.]+)"/.exec(svg);
  if (!head) throw new Error("the SVG has no size");
  const [w, h] = [Math.ceil(Number(head[1])), Math.ceil(Number(head[2]))];
  const url = URL.createObjectURL(new Blob([svg], { type: "image/svg+xml" }));
  try {
    const img = new Image();
    img.src = url;
    await img.decode();
    const canvas = Object.assign(document.createElement("canvas"), { width: w, height: h });
    (canvas.getContext("2d") as CanvasRenderingContext2D).drawImage(img, 0, 0, w, h);
    const png = await new Promise<Blob | null>((ok) => canvas.toBlob(ok, "image/png"));
    if (!png) throw new Error("the canvas wrote no PNG");
    return { png, w, h };
  } finally {
    URL.revokeObjectURL(url);
  }
}

/** One document: its live viewer and the export buttons, the made files and the figure beside the canvas. */
export function ExportCard({ doc, label, note }: { doc: SeqViewSpecDoc; label: string; note: string }) {
  const [scale, setScale] = useState(1);
  const [background, setBackground] = useState(false);
  const [side, setSide] = useState(false);
  const [busy, setBusy] = useState<"svg" | "png" | null>(null);
  const [error, setError] = useState("");
  const [files, setFiles] = useState<{ svg?: Made; png?: Made }>({});
  const [shown, setShown] = useState<"svg" | "png">("svg");
  /** What the registered action's `offer` hands the file to: set by each click. */
  const offered = useRef<((svg: string, name: string) => void) | null>(null);
  /** The options the action reads on each run (its `render` is spread then): the checkbox is live. */
  const opts = useRef({ background });
  opts.current.background = background;
  const registry = useMemo<SeqViewerRegistryEntry[]>(
    () => [
      {
        actions: [
          browserExportAction({
            // A click here takes the file; the viewer's own Export menu (the same action) downloads it as the UI does.
            offer: (svg, name) => (offered.current ? offered.current(svg, name) : downloadSvg(svg, name)),
            render: {
              get background() {
                return opts.current.background;
              },
            },
          }),
        ],
      },
    ],
    [],
  );
  const options = useMemo(() => ({ registry, parts: ["Export"] as const }), [registry]);
  const m = useCaseContext(doc, options);
  const ui = useStoryUi(m);
  const made = useRef<Made[]>([]);
  useEffect(
    () => () => {
      for (const x of made.current) URL.revokeObjectURL(x.url);
    },
    [],
  );

  const keep = (x: Made) => {
    const old = made.current.filter((y) => y.kind === x.kind);
    for (const y of old) URL.revokeObjectURL(y.url);
    made.current = [...made.current.filter((y) => y.kind !== x.kind), x];
    setFiles((p) => ({ ...p, [x.kind]: x }));
    setShown(x.kind);
  };

  const run = (kind: "svg" | "png") => {
    const surface = m?.surfaces.values().next().value;
    if (!m || !surface || busy) return;
    const width = Math.round(surface.element.clientWidth || 1000);
    const t0 = performance.now();
    setBusy(kind);
    setError("");
    const watchdog = setTimeout(() => {
      offered.current = null;
      setBusy(null);
      setError("The action offered no file in 120 s (see the viewer's issues).");
    }, 120_000);
    offered.current = async (svg, name) => {
      clearTimeout(watchdog);
      offered.current = null;
      const svgMs = performance.now() - t0;
      try {
        if (kind === "svg") {
          const blob = new Blob([svg], { type: "image/svg+xml" });
          const head = /<svg[^>]* width="([\d.]+)" height="([\d.]+)"/.exec(svg);
          keep({
            kind,
            name,
            url: URL.createObjectURL(blob),
            bytes: blob.size,
            ms: svgMs,
            detail: `render ${ms(svgMs)}`,
            size: head ? `${head[1]} × ${head[2]}` : "?",
          });
        } else {
          const t1 = performance.now();
          const r = await rasterize(svg);
          const rasterMs = performance.now() - t1;
          keep({
            kind,
            name: name.replace(/\.svg$/, ".png"),
            url: URL.createObjectURL(r.png),
            bytes: r.png.size,
            ms: svgMs + rasterMs,
            detail: `render ${ms(svgMs)} + raster ${ms(rasterMs)}`,
            size: `${r.w} × ${r.h}`,
          });
        }
      } catch (e) {
        setError(String(e));
      } finally {
        setBusy(null);
      }
    };
    m.ctx.actions.run("view.exportSvg", { width, scale });
  };

  const style = { "--seq-section-max-height": "560px" } as CSSProperties;
  const figure = files[shown];
  const viewer = m ? (
    <SeqViewerProvider context={m.ctx} ui={ui}>
      <PrismCanvas>
        <SeqViewerLayout />
      </PrismCanvas>
    </SeqViewerProvider>
  ) : (
    <div className={f.status}>Loading…</div>
  );
  const exported = (
    <div className={f.stack}>
      <div className={f.subheading}>
        {figure
          ? `The exported ${figure.kind.toUpperCase()} (${figure.size} px), drawn here as an image`
          : "The exported file appears here"}
      </div>
      {figure ? (
        <img className={t.figure} src={figure.url} alt={`${label}, exported as ${figure.kind}`} />
      ) : (
        <div className={f.none}>Press Export SVG or Export PNG.</div>
      )}
    </div>
  );
  return (
    <section className={f.card} style={style}>
      <h2 className={f.heading}>{label}</h2>
      <div className={f.kicker}>{note}</div>
      <div className={f.row}>
        <TextAction tone="strong" disabled={!m || !!busy} onClick={() => run("svg")}>
          {busy === "svg" ? "Exporting…" : "Export SVG"}
        </TextAction>
        <TextAction tone="strong" disabled={!m || !!busy} onClick={() => run("png")}>
          {busy === "png" ? "Exporting…" : "Export PNG"}
        </TextAction>
        <div className={prism.narrow}>
          <Choices
            label="Scale"
            choices={[1, 2].map((x) => ({ id: String(x), label: `${x}×` }))}
            active={String(scale)}
            onPick={(id) => setScale(Number(id))}
          />
        </div>
        <label className={t.checks}>
          <input type="checkbox" checked={background} onChange={(e) => setBackground(e.target.checked)} />
          background (the surface color)
        </label>
        <label className={t.checks}>
          <input type="checkbox" checked={side} onChange={(e) => setSide(e.target.checked)} />
          side by side
        </label>
      </div>
      {error && <div className={f.status}>{error}</div>}
      <table className={t.files}>
        <thead>
          <tr>
            <th>File</th>
            <th>Size</th>
            <th>Pixels</th>
            <th>Time</th>
            <th />
            <th />
          </tr>
        </thead>
        <tbody>
          {(["svg", "png"] as const).map((k) => {
            const x = files[k];
            return (
              <tr key={k}>
                <td>{x ? x.name : `.${k}`}</td>
                <td>{x ? kb(x.bytes) : "–"}</td>
                <td>{x ? x.size : "–"}</td>
                <td>{x ? `${ms(x.ms)} (${x.detail})` : "–"}</td>
                <td>
                  {x && (
                    <a className={t.download} href={x.url} download={x.name}>
                      Download
                    </a>
                  )}
                </td>
                <td>
                  {x && shown !== k && (
                    <TextAction tone="quiet" onClick={() => setShown(k)}>
                      Show
                    </TextAction>
                  )}
                </td>
              </tr>
            );
          })}
        </tbody>
      </table>
      {side ? (
        <div className={t.pair}>
          <div className={f.stack}>
            <div className={f.subheading}>The live viewer (canvas)</div>
            {viewer}
          </div>
          {exported}
        </div>
      ) : (
        <>
          <div className={f.stack}>
            <div className={f.subheading}>The live viewer (canvas)</div>
            {viewer}
          </div>
          {exported}
        </>
      )}
    </section>
  );
}

/** The page's frame, as the other Seq* viewer stories have it. */
export function ExportFrame({ children }: { children: ReactNode }) {
  return (
    <PrismPage
      title="Exporting a figure: SVG and PNG"
      caption={
        <>
          Each viewer is live (pan, zoom). The buttons run <code>view.exportSvg</code>, the action the viewer's Export
          part runs (<code>browserExportAction</code>, <code>@molstar/seq-viewer/dom</code>), at the live canvas's
          width; the SVG is drawn beside it as an image so the two can be compared by eye. The PNG is that SVG drawn by
          the page. Time is from the click to the file. The file is transparent unless background is checked: in the
          dark theme its inks are light.
        </>
      }
    >
      {children}
    </PrismPage>
  );
}
