/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Sebastian Bittrich <sebastian.m.bittrich@gmail.com>
 */

import * as path from 'path';
import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { createJob, JobEntry } from '../jobs.js';
import { HttpError } from '../errors.js';
import {
  createStructureWrapperFromJobEntry,
  estimateRetainedSize,
  resolveStructures,
  StructureCache,
} from '../structure-wrapper.js';
import type { ResultWriter } from '../../utils/writer.js';

const examples = path.resolve(import.meta.dirname, '../../../../../data/examples');

function entry(file: string) {
  const e = JobEntry({
    entryId: path.join(examples, file),
    queryName: 'full',
    queryParams: {},
    copyAllCategories: false,
  });
  createJob({ entries: [e], writer: {} as ResultWriter });
  return e;
}

describe('model server structure wrapper', () => {
  beforeEach(() => {
    vi.spyOn(console, 'log').mockImplementation(() => {});
    vi.spyOn(console, 'error').mockImplementation(() => {});
  });
  afterEach(() => {
    StructureCache.expireAll();
  });

  it('reads an entry once for concurrent jobs and caches it', async () => {
    const [a, b] = await Promise.all([
      createStructureWrapperFromJobEntry(entry('1cbs_full.bcif'), undefined),
      createStructureWrapperFromJobEntry(entry('1cbs_full.bcif'), undefined),
    ]);
    expect(a).toBe(b);
    expect(StructureCache.stats.count).toBe(1);
    expect(await createStructureWrapperFromJobEntry(entry('1cbs_full.bcif'), undefined)).toBe(a);
  });

  it('builds the structure of a model once', async () => {
    const wrapper = await createStructureWrapperFromJobEntry(entry('1cbs_full.bcif'), undefined);
    const [first] = await resolveStructures(wrapper);
    const [second] = await resolveStructures(wrapper);
    expect(first.elementCount).toBeGreaterThan(0);
    expect(second).toBe(first);
  });

  it('attaches the custom properties once per model, and survives a failing one', async () => {
    const attach = vi.fn(() => [Promise.reject(new Error('unavailable')), Promise.resolve()]);
    const e = entry('1cbs_full.bcif');
    const wrapper = await createStructureWrapperFromJobEntry(e, attach);
    await resolveStructures(wrapper);
    await resolveStructures(wrapper);
    expect(attach).toHaveBeenCalledTimes(1);
  });

  it('estimates what a wrapper retains from its source size', async () => {
    const wrapper = await createStructureWrapperFromJobEntry(entry('1cbs_full.bcif'), undefined);
    expect(wrapper.isBinary).toBe(true);
    expect(wrapper.approximateSize).toBeGreaterThan(2 * 1024 * 1024);
    expect(estimateRetainedSize(100 * 1024 * 1024, false)).toBeGreaterThan(
      estimateRetainedSize(100 * 1024 * 1024, true),
    );
  });

  it('answers a missing entry with 404 and does not cache the failure', async () => {
    const missing = entry('does-not-exist.cif');
    const error = await createStructureWrapperFromJobEntry(missing, undefined).catch((e) => e);
    expect(error).toBeInstanceOf(HttpError);
    expect(error.status).toBe(404);
    expect(StructureCache.stats.count).toBe(0);
  });
});
