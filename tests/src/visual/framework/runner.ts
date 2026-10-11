/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { WebGPUUnsupportedError } from '@molstar/graphics/gl/webgpu/errors';
import { type BackendInstance, type BackendProfile } from './backend.js';
import { type CapabilityKey, type CapabilityReport, getMissingCapabilities } from './capabilities.js';
import {
  type CompareMetrics,
  DefaultTolerance,
  type Tolerance,
  compareImages,
  isUniform,
  isWithinTolerance,
} from './compare.js';
import { renderTest, type RenderOutput } from './harness.js';
import { loadPng } from './image-io.js';
import { findKnownDifference, type KnownDifference } from './known-differences.js';
import { testIdToFileName, type VisualTestCase } from './test-case.js';

export type TestStatus = 'pass' | 'fail' | 'known-diff' | 'skipped' | 'error';

export type SnapshotStatus = 'match' | 'mismatch' | 'missing' | 'disabled';

export interface SnapshotResult {
  status: SnapshotStatus;
  image?: ImageData;
  diff?: ImageData;
  metrics?: CompareMetrics;
}

export interface TestResult {
  id: string;
  description?: string;
  tags: readonly string[];
  status: TestStatus;
  reason?: string;
  knownDifference?: KnownDifference;
  tolerance: Tolerance;
  metrics?: CompareMetrics;
  reference?: RenderOutput;
  target?: RenderOutput;
  diff?: ImageData;
  snapshot: SnapshotResult;
  errors: string[];
  warnings: string[];
}

export interface RunOptions {
  target: BackendProfile;
  reference: BackendProfile;
  /** Compare reference renders against committed reference snapshots */
  snapshots: boolean;
  snapshotTolerance: Tolerance;
  referencesUrl: string;
}

/** Loose by default: snapshots come from a different GPU/driver/browser */
export const DefaultSnapshotTolerance: Tolerance = {
  threshold: 0.15,
  maxDiffRatio: 0.01,
  maxDiffPixels: Infinity,
};

export interface RunSummary {
  target: string;
  reference: string;
  reports: { target?: CapabilityReport; reference?: CapabilityReport };
  counts: Record<TestStatus, number>;
  results: TestResult[];
  done: boolean;
}

function createCounts(): Record<TestStatus, number> {
  return { pass: 0, fail: 0, 'known-diff': 0, skipped: 0, error: 0 };
}

export class VisualTestRunner {
  private target: BackendInstance | undefined = undefined;
  private reference: BackendInstance | undefined = undefined;

  readonly summary: RunSummary;

  private async getInstances(isolate: boolean, size: { width: number; height: number }) {
    const { target, reference } = this.options;
    if (isolate) {
      const t = await target.create(size);
      const r = reference.id === target.id ? t : await reference.create(size);
      return {
        target: t,
        reference: r,
        dispose: () => {
          t?.dispose();
          if (r !== t) r?.dispose();
        },
      };
    }
    if (!this.target) {
      this.target = await target.create(size);
      this.reference = reference.id === target.id ? this.target : await reference.create(size);
      this.summary.reports.target = this.target?.report;
      this.summary.reports.reference = this.reference?.report;
    }
    return { target: this.target, reference: this.reference, dispose: () => {} };
  }

  private async checkSnapshot(test: VisualTestCase, reference: RenderOutput | undefined): Promise<SnapshotResult> {
    if (!this.options.snapshots || !reference) return { status: 'disabled' };
    const image = await loadPng(`${this.options.referencesUrl}${testIdToFileName(test.id)}`);
    if (!image) return { status: 'missing' };
    const { metrics, diff } = compareImages(image, reference.image, this.options.snapshotTolerance.threshold);
    const status = isWithinTolerance(metrics, this.options.snapshotTolerance) ? 'match' : 'mismatch';
    return { status, image, diff, metrics };
  }

  async runTest(test: VisualTestCase): Promise<TestResult> {
    const knownDifference = findKnownDifference(this.options.target.id, test.id);
    const tolerance: Tolerance = {
      ...DefaultTolerance,
      ...test.tolerance,
      ...(knownDifference?.action.kind === 'tolerance' ? knownDifference.action.tolerance : undefined),
    };
    const result: TestResult = {
      id: test.id,
      description: test.description,
      tags: test.tags ?? [],
      status: 'skipped',
      knownDifference,
      tolerance,
      snapshot: { status: 'disabled' },
      errors: [],
      warnings: [],
    };

    if (knownDifference?.action.kind === 'unsupported') {
      result.reason = `Unsupported: ${knownDifference.reason}`;
      return result;
    }

    const instances = await this.getInstances(!!test.isolate, test.size ?? { width: 256, height: 256 });
    try {
      const { target, reference } = instances;
      if (!target || !reference) {
        result.reason = `Backend not available: ${!target ? this.options.target.id : this.options.reference.id}`;
        return result;
      }

      const requires = test.requires ?? [];
      const missing = new Set<CapabilityKey>([
        ...getMissingCapabilities(target.capabilities, requires),
        ...getMissingCapabilities(reference.capabilities, requires),
      ]);
      if (missing.size > 0) {
        result.reason = `Missing capabilities: ${Array.from(missing).join(', ')}`;
        return result;
      }

      const { fallback } = test;
      const targetLacks = !!fallback && getMissingCapabilities(target.capabilities, fallback.requires).length > 0;
      const referenceLacks = !!fallback && getMissingCapabilities(reference.capabilities, fallback.requires).length > 0;
      const withFallback = (lacks: boolean) =>
        fallback && (targetLacks || referenceLacks) && !lacks ? { ...test, ...fallback.overrides } : test;
      if (fallback && targetLacks !== referenceLacks) {
        result.reason = `Compared against the fallback for missing ${fallback.requires.join(', ')}`;
      }

      result.reference = await renderTest(reference, withFallback(referenceLacks));
      result.target = await renderTest(target, withFallback(targetLacks));
      result.errors.push(
        ...result.reference.errors.map((e) => `[reference] ${e}`),
        ...result.target.errors.map((e) => `[target] ${e}`),
      );
      result.warnings.push(
        ...result.reference.warnings.map((e) => `[reference] ${e}`),
        ...result.target.warnings.map((e) => `[target] ${e}`),
      );
      if (isUniform(result.reference.image)) result.errors.push('[reference] blank render (uniform image)');
      if (isUniform(result.target.image)) result.errors.push('[target] blank render (uniform image)');

      const { metrics, diff } = compareImages(result.reference.image, result.target.image, tolerance.threshold);
      result.metrics = metrics;
      result.diff = diff;
      result.snapshot = await this.checkSnapshot(test, result.reference);

      if (result.errors.length > 0) {
        result.status = 'error';
      } else if (isWithinTolerance(metrics, tolerance)) {
        result.status = 'pass';
      } else if (knownDifference?.action.kind === 'expected-diff') {
        result.status = 'known-diff';
        result.reason = knownDifference.reason;
      } else {
        result.status = 'fail';
      }
    } catch (e) {
      if (e instanceof WebGPUUnsupportedError) {
        result.status = 'skipped';
        result.reason = e.message;
        return result;
      }
      result.status = 'error';
      result.errors.push(e instanceof Error ? `${e.message}\n${e.stack}` : String(e));
    } finally {
      instances.dispose();
    }
    return result;
  }

  async run(tests: readonly VisualTestCase[], onResult?: (result: TestResult) => void) {
    for (const test of tests) {
      const result = await this.runTest(test);
      this.summary.results.push(result);
      this.summary.counts[result.status] += 1;
      onResult?.(result);
      // yield so the UI can update between tests
      await new Promise((r) => setTimeout(r, 0));
    }
    this.summary.done = true;
    return this.summary;
  }

  dispose() {
    if (this.reference && this.reference !== this.target) this.reference.dispose();
    this.target?.dispose();
    this.target = this.reference = undefined;
  }

  constructor(private readonly options: RunOptions) {
    this.summary = {
      target: options.target.id,
      reference: options.reference.id,
      reports: {},
      counts: createCounts(),
      results: [],
      done: false,
    };
  }
}

/** JSON-serializable summary without image data */
export function summaryToJson(summary: RunSummary) {
  return {
    target: summary.target,
    reference: summary.reference,
    reports: summary.reports,
    counts: summary.counts,
    done: summary.done,
    results: summary.results.map((r) => ({
      id: r.id,
      tags: r.tags,
      status: r.status,
      reason: r.reason,
      knownDifference: r.knownDifference && { reason: r.knownDifference.reason, doc: r.knownDifference.doc },
      tolerance: r.tolerance,
      metrics: r.metrics,
      snapshot: { status: r.snapshot.status, metrics: r.snapshot.metrics },
      timeMs: { reference: r.reference?.timeMs, target: r.target?.timeMs },
      errors: r.errors,
      warnings: r.warnings,
    })),
  };
}
