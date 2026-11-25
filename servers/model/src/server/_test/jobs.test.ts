/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Sebastian Bittrich <sebastian.m.bittrich@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import type { ResultWriter } from '../../utils/writer.js';
import { createJob, type Job, JobRunner, type JobRunnerParams } from '../jobs.js';

function writer(state = { cancelled: false }): ResultWriter {
  return {
    beginEntry() {},
    endEntry() {},
    writeBinary: () => true,
    writeString: () => true,
    end() {},
    isCancelled: () => state.cancelled,
  };
}

function job(state?: { cancelled: boolean }) {
  return createJob({ entries: [], writer: writer(state) });
}

/** A runner whose jobs finish when the test says so. */
function setup(params: JobRunnerParams) {
  const started: Job[] = [];
  const skipped: Job[] = [];
  const finish = new Map<Job, { resolve: () => void; reject: (e: any) => void }>();
  let running = 0;
  let maxRunning = 0;
  const runner = new JobRunner(
    (j) => {
      started.push(j);
      running++;
      maxRunning = Math.max(maxRunning, running);
      return new Promise<void>((resolve, reject) => finish.set(j, { resolve, reject })).finally(() => running--);
    },
    () => params,
    (j) => skipped.push(j),
  );
  const settle = async (j: Job, error?: any) => {
    if (error) finish.get(j)!.reject(error);
    else finish.get(j)!.resolve();
    // let the runner's `finally` and its next dispatch (two `setImmediate` turns) run
    for (let i = 0; i < 3; i++) await new Promise((r) => setImmediate(r));
  };
  return { runner, started, skipped, settle, maxRunning: () => maxRunning };
}

describe('model server job runner', () => {
  it('runs at most maxConcurrentJobs at a time, in order', async () => {
    const { runner, started, settle, maxRunning } = setup({ maxConcurrentJobs: 2, maxQueueLength: 10 });
    const jobs = [job(), job(), job(), job(), job()];
    for (const j of jobs) expect(runner.add(j)).toBe(true);

    expect(started).toEqual(jobs.slice(0, 2));
    expect(runner.stats).toEqual({ queued: 3, running: 2 });

    await settle(jobs[1]);
    expect(started).toEqual(jobs.slice(0, 3));
    await settle(jobs[0]);
    await settle(jobs[2]);
    expect(started).toEqual(jobs.slice(0, 5));
    await settle(jobs[3]);
    await settle(jobs[4]);

    expect(maxRunning()).toBe(2);
    expect(runner.stats).toEqual({ queued: 0, running: 0 });
  });

  it('refuses jobs once the queue is full', async () => {
    // molstar#1816: `maxQueueLength` was never enforced, and every request that found the queue empty started a
    // worker of its own.
    const { runner, settle } = setup({ maxConcurrentJobs: 1, maxQueueLength: 2 });
    const first = job();
    expect(runner.add(first)).toBe(true); // running
    expect(runner.add(job())).toBe(true); // queued
    expect(runner.add(job())).toBe(true); // queued
    expect(runner.add(job())).toBe(false);
    expect(runner.stats).toEqual({ queued: 2, running: 1 });

    await settle(first);
    expect(runner.add(job())).toBe(true);
  });

  it('skips jobs whose client is gone', async () => {
    const { runner, started, skipped, settle } = setup({ maxConcurrentJobs: 1, maxQueueLength: 10 });
    const gone = { cancelled: false };
    const first = job();
    const second = job(gone);
    const third = job();
    runner.add(first);
    runner.add(second);
    runner.add(third);

    gone.cancelled = true;
    await settle(first);
    expect(skipped).toEqual([second]);
    expect(started).toEqual([first, third]);
  });

  it('keeps going after a job fails', async () => {
    const { runner, started, settle } = setup({ maxConcurrentJobs: 1, maxQueueLength: 10 });
    const first = job();
    const second = job();
    runner.add(first);
    runner.add(second);

    await settle(first, new Error('failed'));
    expect(started).toEqual([first, second]);
    await settle(second);
    expect(runner.stats).toEqual({ queued: 0, running: 0 });
  });

  it('drains: refuses new jobs and resolves when the running ones are done', async () => {
    const { runner, settle } = setup({ maxConcurrentJobs: 1, maxQueueLength: 10 });
    const first = job();
    const second = job();
    runner.add(first);
    runner.add(second);

    let drained = false;
    const draining = runner.drain().then(() => (drained = true));
    expect(runner.isAccepting).toBe(false);
    expect(runner.add(job())).toBe(false);

    await settle(first);
    expect(drained).toBe(false);
    await settle(second);
    await draining;
    expect(drained).toBe(true);
  });

  it('drains at once when idle', async () => {
    const { runner } = setup({ maxConcurrentJobs: 1, maxQueueLength: 10 });
    await expect(runner.drain()).resolves.toBeUndefined();
  });
});
