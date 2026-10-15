/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { BackendProfiles } from '../framework/backend.js';
import { CapabilityKeys, type CapabilityReport } from '../framework/capabilities.js';
import { downloadPng } from '../framework/image-io.js';
import { type RunSummary, summaryToJson, type TestResult, type TestStatus } from '../framework/runner.js';
import { testIdToFileName } from '../framework/test-case.js';

export interface GalleryOptions {
  target: string;
  reference: string;
  filter: string;
  tags: string;
  snapshots: boolean;
  zoom: number;
}

function el<K extends keyof HTMLElementTagNameMap>(
  tag: K,
  props: Partial<HTMLElementTagNameMap[K]> = {},
  children: (Node | string)[] = [],
) {
  const e = document.createElement(tag);
  Object.assign(e, props);
  for (const c of children) e.append(c);
  return e;
}

function fmt(v: number | undefined, digits = 3) {
  if (v === undefined) return '-';
  if (!isFinite(v)) return '∞';
  return v.toFixed(digits);
}

function imageCanvas(image: ImageData, zoom: number) {
  const canvas = el('canvas', { width: image.width, height: image.height });
  canvas.getContext('2d')!.putImageData(image, 0, 0);
  canvas.style.width = `${image.width * zoom}px`;
  canvas.style.height = `${image.height * zoom}px`;
  return canvas;
}

function figure(caption: string, image: ImageData, zoom: number) {
  return el('figure', {}, [imageCanvas(image, zoom), el('figcaption', { textContent: caption })]);
}

/** Alternates between two images to make small differences visible */
function flickerFigure(a: ImageData, b: ImageData, labelA: string, labelB: string, zoom: number) {
  const canvas = imageCanvas(a, zoom);
  const caption = el('figcaption', { textContent: `flicker: ${labelA}` });
  const ctx = canvas.getContext('2d')!;
  let showA = true;
  const handle = setInterval(() => {
    if (!canvas.isConnected) {
      clearInterval(handle);
      return;
    }
    showA = !showA;
    ctx.putImageData(showA ? a : b, 0, 0);
    caption.textContent = `flicker: ${showA ? labelA : labelB}`;
  }, 500);
  return el('figure', {}, [canvas, caption]);
}

export class Gallery {
  private readonly summaryEl = el('div', { className: 'summary' });
  private readonly reportsEl = el('details', { className: 'reports' });
  private readonly listEl = el('div');
  private readonly cards = new Map<string, HTMLElement>();
  private failedOnly = false;

  constructor(
    parent: HTMLElement,
    private readonly options: GalleryOptions,
    private readonly summary: RunSummary,
    private readonly rerun: (id: string) => Promise<TestResult>,
  ) {
    parent.append(this.createHeader(), el('main', {}, [this.reportsEl, this.listEl]));
    this.updateSummary();
  }

  private createHeader() {
    const profileSelect = (value: string) => {
      const s = el('select');
      for (const p of BackendProfiles)
        s.append(el('option', { value: p.id, textContent: p.id, title: p.description, selected: p.id === value }));
      return s;
    };
    const target = profileSelect(this.options.target);
    const reference = profileSelect(this.options.reference);
    const filter = el('input', {
      type: 'text',
      value: this.options.filter,
      placeholder: 'globs, e.g. transparency/**',
    });
    const tags = el('input', { type: 'text', value: this.options.tags, placeholder: 'tags, e.g. smoke' });
    const snapshots = el('input', { type: 'checkbox', checked: this.options.snapshots });
    const zoom = el('select');
    for (const z of [1, 2, 3, 4])
      zoom.append(el('option', { value: `${z}`, textContent: `${z}x`, selected: z === this.options.zoom }));

    const run = el('button', { textContent: 'Run' });
    run.onclick = () => {
      const params = new URLSearchParams({
        target: target.value,
        reference: reference.value,
        filter: filter.value,
        tags: tags.value,
        snapshots: snapshots.checked ? '1' : '0',
        zoom: zoom.value,
      });
      window.location.search = params.toString();
    };

    const failedOnly = el('input', { type: 'checkbox' });
    failedOnly.onchange = () => {
      this.failedOnly = failedOnly.checked;
      for (const r of this.summary.results) this.applyVisibility(r);
    };

    const json = el('button', { textContent: 'Download JSON' });
    json.onclick = () => {
      const blob = new Blob([JSON.stringify(summaryToJson(this.summary), undefined, 2)], { type: 'application/json' });
      const a = el('a', {
        href: URL.createObjectURL(blob),
        download: `visual-tests-${this.summary.target}-vs-${this.summary.reference}.json`,
      });
      a.click();
      setTimeout(() => URL.revokeObjectURL(a.href), 1000);
    };

    return el('header', {}, [
      el('h1', { textContent: 'Mol* Visual Tests' }),
      el('label', {}, ['Target', target]),
      el('label', {}, ['Reference', reference]),
      el('label', {}, ['Filter', filter]),
      el('label', {}, ['Tags', tags]),
      el('label', {}, [snapshots, 'Snapshots']),
      el('label', {}, ['Zoom', zoom]),
      run,
      el('label', {}, [failedOnly, 'Problems only']),
      json,
      this.summaryEl,
    ]);
  }

  updateSummary() {
    const { counts, done, results } = this.summary;
    this.summaryEl.replaceChildren(
      ...(Object.keys(counts) as TestStatus[]).map((k) =>
        el('span', {}, [el('span', { className: `status status-${k}`, textContent: k }), ` ${counts[k]}`]),
      ),
      el('span', { textContent: done ? `done (${results.length})` : `running... (${results.length})` }),
    );
  }

  setReports(reports: { target?: CapabilityReport; reference?: CapabilityReport }) {
    const { target, reference } = reports;
    const row = (label: string, a: string, b: string) =>
      el('tr', {}, [
        el('th', { textContent: label }),
        el('td', { textContent: a, className: a !== b ? 'mismatch' : '' }),
        el('td', { textContent: b, className: a !== b ? 'mismatch' : '' }),
      ]);
    const rows = [
      row('backend', target?.backend ?? '-', reference?.backend ?? '-'),
      row('version', target?.version ?? '-', reference?.version ?? '-'),
      row('renderer', target?.renderer ?? '-', reference?.renderer ?? '-'),
      row('vendor', target?.vendor ?? '-', reference?.vendor ?? '-'),
      ...CapabilityKeys.map((k) => row(k, `${target?.capabilities[k] ?? '-'}`, `${reference?.capabilities[k] ?? '-'}`)),
      ...Object.keys(target?.limits ?? reference?.limits ?? {}).map((k) =>
        row(k, `${target?.limits[k] ?? '-'}`, `${reference?.limits[k] ?? '-'}`),
      ),
    ];
    this.reportsEl.replaceChildren(
      el('summary', {
        textContent: `Capabilities: ${this.summary.target} (target) vs ${this.summary.reference} (reference) — ${target?.userAgent ?? navigator.userAgent}`,
      }),
      el('table', {}, [
        el('tr', {}, [
          el('th'),
          el('th', { textContent: `target: ${this.summary.target}` }),
          el('th', { textContent: `reference: ${this.summary.reference}` }),
        ]),
        ...rows,
      ]),
    );
  }

  private applyVisibility(result: TestResult) {
    const card = this.cards.get(result.id);
    if (!card) return;
    const problem =
      result.status === 'fail' ||
      result.status === 'error' ||
      result.status === 'known-diff' ||
      result.snapshot.status === 'mismatch';
    card.classList.toggle('hidden', this.failedOnly && !problem);
  }

  addResult(result: TestResult) {
    const zoom = this.options.zoom;
    const { target, reference } = this.summary;

    const header = el('div', { className: 'card-header' }, [
      el('span', { className: `status status-${result.status}`, textContent: result.status }),
      el('span', { className: 'id', textContent: result.id }),
      ...result.tags.map((t) => el('span', { textContent: `#${t}` })),
    ]);

    if (result.snapshot.status !== 'disabled') {
      header.append(
        el('span', {
          className: `snapshot-${result.snapshot.status}`,
          textContent: `snapshot: ${result.snapshot.status}`,
        }),
      );
    }

    const rerun = el('button', { textContent: 'Rerun' });
    rerun.onclick = async () => {
      rerun.disabled = true;
      this.replaceResult(await this.rerun(result.id));
    };
    header.append(rerun);

    if (result.reference) {
      const approve = el('button', {
        textContent: 'Approve reference',
        title: `Download as ${testIdToFileName(result.id)}; commit to src/tests/visual/references/`,
      });
      approve.onclick = () => downloadPng(result.reference!.image, testIdToFileName(result.id));
      header.append(approve);
    }

    const card = el('div', { className: 'card' }, [header]);
    if (result.description) card.append(el('div', { className: 'description', textContent: result.description }));
    if (result.reason) card.append(el('div', { className: 'reason', textContent: result.reason }));
    if (result.knownDifference) {
      card.append(
        el('div', {
          className: 'reason',
          textContent: `Known difference: ${result.knownDifference.reason}${result.knownDifference.doc ? ` (${result.knownDifference.doc})` : ''}`,
        }),
      );
    }

    const m = result.metrics;
    if (m) {
      const t = result.tolerance;
      card.append(
        el('div', {
          className: 'metrics',
          textContent:
            `diff pixels ${m.diffPixels} (${fmt(m.diffRatio * 100, 3)}%, max ${fmt(t.maxDiffRatio * 100, 3)}%) · max Δ ${fmt(m.maxDelta)} · mean Δ ${fmt(m.meanDelta, 4)} · PSNR ${fmt(m.psnr, 1)} dB · threshold ${t.threshold}` +
            ` · time ref ${fmt(result.reference?.timeMs, 0)} ms / target ${fmt(result.target?.timeMs, 0)} ms`,
        }),
      );
    }
    const sm = result.snapshot.metrics;
    if (sm) {
      card.append(
        el('div', {
          className: 'metrics',
          textContent: `snapshot diff pixels ${sm.diffPixels} (${fmt(sm.diffRatio * 100, 3)}%) · max Δ ${fmt(sm.maxDelta)} · PSNR ${fmt(sm.psnr, 1)} dB`,
        }),
      );
    }

    const images = el('div', { className: 'images' });
    if (result.reference) images.append(figure(`reference: ${reference}`, result.reference.image, zoom));
    if (result.target) images.append(figure(`target: ${target}`, result.target.image, zoom));
    if (result.diff) images.append(figure('diff (red: over threshold, yellow: below)', result.diff, zoom));
    if (result.reference && result.target)
      images.append(flickerFigure(result.reference.image, result.target.image, reference, target, zoom));
    if (result.snapshot.image) images.append(figure('committed snapshot', result.snapshot.image, zoom));
    if (result.snapshot.diff) images.append(figure('snapshot diff', result.snapshot.diff, zoom));
    card.append(images);

    if (result.errors.length) card.append(el('pre', { textContent: result.errors.join('\n') }));
    if (result.warnings.length) {
      card.append(
        el('details', {}, [
          el('summary', { textContent: `${result.warnings.length} warning(s)` }),
          el('pre', { textContent: result.warnings.join('\n') }),
        ]),
      );
    }

    const existing = this.cards.get(result.id);
    if (existing) existing.replaceWith(card);
    else this.listEl.append(card);
    this.cards.set(result.id, card);
    this.applyVisibility(result);
    this.updateSummary();
  }

  private replaceResult(result: TestResult) {
    const i = this.summary.results.findIndex((r) => r.id === result.id);
    if (i >= 0) {
      this.summary.counts[this.summary.results[i].status] -= 1;
      this.summary.results[i] = result;
      this.summary.counts[result.status] += 1;
    }
    this.addResult(result);
  }
}
