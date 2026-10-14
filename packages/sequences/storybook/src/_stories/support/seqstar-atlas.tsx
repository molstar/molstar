/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { SeqViewSpecDoc } from "@molstar/seq-core/format";
import { Dirty } from "@molstar/seq-viewer/draw";
import { SeqViewerLayout, SeqViewerProvider } from "@molstar/seq-viewer/react";
import { type CSSProperties, useEffect, useMemo, useRef, useState } from "react";
import { PrismPage } from "../support/seqstar-prism";
import { Choices, TextAction } from "../ui";
import { atlasDoc, atlasLoaders, ATLAS_COLUMNS as COLUMNS, ATLAS_ROWS as ROWS } from "./seqstar-docs-streamed";
import f from "./seqstar-format.module.css";
import { PrismCanvas, prism } from "./seqstar-prism";
import t from "./seqstar-tools.module.css";
import { type Mounted, useCaseContext, useStoryUi } from "./seqstar-viewer";

let DOC: SeqViewSpecDoc | null = null;

export interface Stats {
  readonly n: number;
  readonly median: number;
  readonly p95: number;
  readonly max: number;
}
function stats(xs: readonly number[]): Stats {
  const s = [...xs].sort((a, b) => a - b);
  const q = (p: number) => s[Math.min(s.length - 1, Math.floor(s.length * p))] ?? 0;
  return { n: s.length, median: q(0.5), p95: q(0.95), max: s[s.length - 1] ?? 0 };
}
const ms = (x: number) => x.toFixed(2);

type Kind = "pan" | "redraw";
const KINDS: Record<Kind, string> = { pan: "Pan", redraw: "Full redraws" };
const FRAMES = 240;
const WARMUP = 12;
/** What a crop is: device pixels of the main canvas at the first rows, read back. */
interface Crop {
  readonly image: ImageData;
  readonly dpr: number;
}

const nextFrame = () => new Promise<void>((ok) => requestAnimationFrame(() => ok()));
const wait = (n: number) => new Promise<void>((ok) => setTimeout(ok, n));

/** The first surface's pieces, once it has drawn. */
function parts(m: Mounted) {
  const [section, surface] = m.surfaces.entries().next().value ?? [];
  const vp = section ? m.ctx.viewports.get(section) : undefined;
  return surface && vp && surface.stats.frames > 0 ? { surface, vp } : null;
}
/** The view the measurement uses: cells 14 px wide (letters), rows 18 px (the row scale 1), the first column. */
function resetView(m: Mounted) {
  const p = parts(m);
  if (!p) return false;
  p.vp.setView(14, 0, "api");
  p.vp.setRowView(1, 0, "api");
  return true;
}
const signature = (d: ImageData) => {
  let h = 0;
  for (let i = 0; i < d.data.length; i += 7) h = (Math.imul(h, 31) + d.data[i]) | 0;
  return h;
};
/** 140 × 54 CSS px of the main canvas, from the middle, over the first rows (below the ruler), in device px. */
function takeCrop(m: Mounted): Crop | null {
  const p = parts(m);
  const canvas = p?.surface.element.querySelector("canvas");
  if (!p || !canvas || canvas.width === 0) return null;
  const dpr = canvas.width / Math.max(1, p.surface.element.clientWidth);
  const [w, h] = [Math.round(140 * dpr), Math.round(54 * dpr)];
  const x = Math.round(canvas.width / 2 - w / 2);
  const y = Math.round(44 * dpr);
  const g = canvas.getContext("2d", { willReadFrequently: true }) as CanvasRenderingContext2D;
  return { image: g.getImageData(x, y, w, h), dpr };
}

/** A viewer with the atlas on or off, its runs, and the crop taken once it settles. Remounted when the atlas changes. */
export function AtlasViewer({
  atlas,
  onResult,
  onCrop,
  onLive,
}: {
  atlas: boolean;
  onResult(kind: Kind, s: Stats): void;
  onCrop(c: Crop): void;
  onLive(text: string): void;
}) {
  DOC ??= atlasDoc();
  const options = useMemo(() => ({ loaders: atlasLoaders, maxHeight: 520 }), []);
  const m = useCaseContext(DOC, options);
  const surfaceOptions = useMemo(() => ({ atlas }), [atlas]);
  const ui = useStoryUi(m, undefined, surfaceOptions);
  const [running, setRunning] = useState<Kind | null>(null);
  const cancelled = useRef(false);
  const callbacks = useRef({ onResult, onCrop, onLive });
  callbacks.current = { onResult, onCrop, onLive };

  // Once the surface has drawn: the measurement's view, then a crop when the pixels stop changing (tiles arrive).
  useEffect(() => {
    if (!m) return;
    cancelled.current = false;
    let last = Number.NaN;
    let steady = 0;
    let started = false;
    let tries = 0;
    const timer = setInterval(() => {
      if (cancelled.current) return;
      if (!started) {
        started = resetView(m);
        return;
      }
      const crop = takeCrop(m);
      if (!crop || ++tries > 40) return clearInterval(timer);
      const sig = signature(crop.image);
      steady = sig === last ? steady + 1 : 0;
      last = sig;
      if (steady >= 3) {
        clearInterval(timer);
        callbacks.current.onCrop(crop);
      }
    }, 400);
    return () => {
      cancelled.current = true;
      clearInterval(timer);
    };
  }, [m]);

  const run = async (kind: Kind) => {
    if (!m || running || !resetView(m)) return;
    setRunning(kind);
    const p = parts(m);
    if (!p) return;
    const { surface, vp } = p;
    await wait(1500); // the tiles of the first view arrive
    const samples: number[] = [];
    for (let i = 0; i < WARMUP + FRAMES && !cancelled.current; i++) {
      const before = surface.stats.frames;
      if (kind === "pan") vp.panBy(i % 40 < 20 ? 37 : -37, 0, "pointer");
      else {
        surface.composer.redraw();
        surface.invalidate(Dirty.all);
      }
      for (let k = 0; k < 3 && surface.stats.frames === before; k++) await nextFrame();
      if (surface.stats.frames === before || i < WARMUP) continue;
      samples.push(surface.stats.lastMs);
      if (samples.length % 10 === 0) {
        const x = stats(samples);
        callbacks.current.onLive(
          `${atlas ? "Atlas on" : "Atlas off"}, ${KINDS[kind].toLowerCase()}: ${x.n} frames, median ${ms(x.median)} ms, p95 ${ms(x.p95)} ms, max ${ms(x.max)} ms`,
        );
      }
    }
    if (samples.length > 0) {
      const x = stats(samples);
      callbacks.current.onResult(kind, x);
      callbacks.current.onLive(
        `${atlas ? "Atlas on" : "Atlas off"}, ${KINDS[kind].toLowerCase()}, done: ${x.n} frames, median ${ms(x.median)} ms, p95 ${ms(x.p95)} ms, max ${ms(x.max)} ms`,
      );
    }
    resetView(m);
    setRunning(null);
  };

  const style = { "--seq-section-max-height": "520px" } as CSSProperties;
  return (
    <div className={f.stack} style={style}>
      <div className={f.row}>
        {(Object.keys(KINDS) as Kind[]).map((k) => (
          <TextAction key={k} tone="strong" disabled={!m || !!running} onClick={() => void run(k)}>
            {running === k ? "Running…" : `Run: ${KINDS[k].toLowerCase()} (${FRAMES} frames)`}
          </TextAction>
        ))}
        <TextAction tone="quiet" disabled={!m || !!running} onClick={() => m && resetView(m)}>
          Reset the view
        </TextAction>
      </div>
      {m ? (
        <SeqViewerProvider context={m.ctx} ui={ui}>
          <PrismCanvas>
            <SeqViewerLayout />
          </PrismCanvas>
        </SeqViewerProvider>
      ) : (
        <div className={f.status}>Loading…</div>
      )}
    </div>
  );
}

/** A crop enlarged by whole pixels (each device pixel a square), drawn by the canvas's own nearest-neighbor scaling. */
function CropCanvas({ crop }: { crop: Crop }) {
  const ref = useRef<HTMLCanvasElement>(null);
  const K = 5;
  useEffect(() => {
    const c = ref.current;
    if (!c) return;
    const { image } = crop;
    const tmp = Object.assign(document.createElement("canvas"), { width: image.width, height: image.height });
    (tmp.getContext("2d") as CanvasRenderingContext2D).putImageData(image, 0, 0);
    c.width = image.width * K;
    c.height = image.height * K;
    const g = c.getContext("2d") as CanvasRenderingContext2D;
    g.imageSmoothingEnabled = false;
    g.drawImage(tmp, 0, 0, c.width, c.height);
  }, [crop]);
  return (
    <canvas
      ref={ref}
      className={t.crop}
      style={{ width: (crop.image.width * K) / crop.dpr, height: (crop.image.height * K) / crop.dpr }}
    />
  );
}

/**
 * Pixels of two crops that differ, the way the tests compare (premultiplied, as the canvas stores them: a pixel of low
 * alpha read back unpremultiplied magnifies a one-unit difference), the largest difference, and those over 1.
 */
function compare(a: Crop, b: Crop) {
  if (a.image.width !== b.image.width || a.image.height !== b.image.height) return null;
  let differ = 0;
  let over = 0;
  let max = 0;
  const [x, y] = [a.image.data, b.image.data];
  for (let i = 0; i < x.length; i += 4) {
    let d = Math.abs(x[i + 3] - y[i + 3]);
    for (let c = 0; c < 3; c++)
      d = Math.max(d, Math.abs(Math.round((x[i + c] * x[i + 3]) / 255) - Math.round((y[i + c] * y[i + 3]) / 255)));
    if (d > 0) differ++;
    if (d > 1) over++;
    max = Math.max(max, d);
  }
  return { differ, over, max, pixels: x.length / 4 };
}

/** The story's page body: the atlas toggle, the viewer, the live line, the results and the crops. */
export function AtlasPage() {
  const [atlas, setAtlas] = useState(true);
  const [live, setLive] = useState("No run yet.");
  const [results, setResults] = useState<Record<string, Stats>>({});
  const [crops, setCrops] = useState<{ on?: Crop; off?: Crop }>({});
  const same = crops.on && crops.off ? compare(crops.on, crops.off) : null;
  const dpr = typeof window === "undefined" ? 1 : window.devicePixelRatio;
  return (
    <PrismPage
      title="The glyph atlas: lettered rows of a large lazy alignment"
      caption={
        <>
          A lazy alignment of {ROWS.toLocaleString("en-US")} rows × {COLUMNS.toLocaleString("en-US")} columns (tiles and
          rows sampled in the page, nothing kept), shown with 18 px rows and 14 px cells: every cell a letter, about
          2,500 a frame. With the atlas (the default) a letter is one <code>drawImage</code> of a cell rasterized once;
          off, every letter is a <code>fillText</code>. Run a pan or full redraws, then toggle and run again. Frame
          times are <code>SurfaceStats.lastMs</code>: the main thread's ms for sizing and the composer's draws, the
          raster after it not included. A pan copies the band and draws only the strip it uncovers, so full redraws show
          the letters' cost most.
        </>
      }
    >
      <div className={prism.controls}>
        <div className={prism.narrow}>
          <Choices
            label="Glyph atlas"
            choices={[
              { id: "on", label: "on (the default)" },
              { id: "off", label: "off (fillText)" },
            ]}
            active={atlas ? "on" : "off"}
            onPick={(id) => {
              const on = id === "on";
              setAtlas(on);
              setLive(`The viewer was remounted with the atlas ${on ? "on" : "off"}.`);
            }}
          />
        </div>
        <span className={f.code}>
          Window {typeof window === "undefined" ? "" : `${window.innerWidth} × ${window.innerHeight}`} CSS px, DPR {dpr}
          .
        </span>
      </div>
      <section className={f.card}>
        <h2 className={f.heading}>Lazy alignment, atlas {atlas ? "on" : "off"}</h2>
        <AtlasViewer
          key={String(atlas)}
          atlas={atlas}
          onLive={setLive}
          onResult={(kind, x) => setResults((r) => ({ ...r, [`${atlas ? "on" : "off"}/${kind}`]: x }))}
          onCrop={(c) => setCrops((p) => ({ ...p, [atlas ? "on" : "off"]: c }))}
        />
        <div className={f.status}>{live}</div>
      </section>
      <section className={f.card}>
        <h2 className={f.heading}>Frame times, main thread (ms)</h2>
        <table className={t.files}>
          <thead>
            <tr>
              <th>Atlas</th>
              <th>Motion</th>
              <th>Frames</th>
              <th>Median</th>
              <th>p95</th>
              <th>Max</th>
            </tr>
          </thead>
          <tbody>
            {(["on", "off"] as const).flatMap((a) =>
              (Object.keys(KINDS) as Kind[]).map((k) => {
                const x = results[`${a}/${k}`];
                return (
                  <tr key={`${a}/${k}`}>
                    <td>{a}</td>
                    <td>{KINDS[k].toLowerCase()}</td>
                    <td>{x ? x.n : "–"}</td>
                    <td>{x ? ms(x.median) : "–"}</td>
                    <td>{x ? ms(x.p95) : "–"}</td>
                    <td>{x ? ms(x.max) : "–"}</td>
                  </tr>
                );
              }),
            )}
          </tbody>
        </table>
      </section>
      <section className={f.card}>
        <h2 className={f.heading}>The same letters, magnified</h2>
        <div className={f.kicker}>
          140 × 54 CSS px of the main canvas at the first rows (cell 14, rows 18 px), read back after the view settled
          and shown at 5× with nearest-neighbor scaling, so each device pixel is a square. The crop is taken after each
          mount: toggle the atlas to fill the other one.
        </div>
        <div className={t.crops}>
          {(["on", "off"] as const).map((a) => {
            const c = crops[a];
            return (
              <div key={a} className={t.cropCell}>
                <div className={f.subheading}>Atlas {a}</div>
                {c ? (
                  <CropCanvas crop={c} />
                ) : (
                  <div className={f.none}>{atlas === (a === "on") ? "Settling…" : "Toggle the atlas to capture."}</div>
                )}
              </div>
            );
          })}
        </div>
        <div className={f.status}>
          {same
            ? `Pixels that differ: ${same.differ.toLocaleString("en-US")} of ${same.pixels.toLocaleString("en-US")}; over 1 level: ${same.over}; the largest difference ${same.max} (of 255), premultiplied as the canvas stores it. The atlas's tests pin at most 1 (_test/canvas/atlas.test.ts): a cell's 8-bit premultiplied ink blends again over the cell's translucent fill.`
            : "Both crops are needed to compare."}
        </div>
      </section>
    </PrismPage>
  );
}
