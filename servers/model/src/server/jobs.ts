/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Sebastian Bittrich <sebastian.m.bittrich@gmail.com>
 */

import { UUID } from '@molstar/core/util';
import { getQueryByName, type QueryDefinition, type QueryName, type QueryParams, type Encoding } from './api.js';
import { LinkedList } from '@molstar/core/data/generic';
import type { ResultWriter } from '@molstar/model-server/utils/writer';
import { Mat4 } from '@molstar/core/math/linear-algebra';
import { HttpError } from './errors.js';

export interface ResponseFormat {
  tarball: boolean;
  encoding: Encoding;
}

export interface Job {
  id: UUID;
  datetime_utc: string;

  entries: JobEntry[];

  responseFormat: ResponseFormat;
  outputFilename?: string;

  writer: ResultWriter;
}

export interface JobDefinition {
  entries: JobEntry[];
  writer: ResultWriter;
  options?: { outputFilename?: string; binary?: boolean; tarball?: boolean; encoding?: Encoding };
}

export interface JobEntry {
  job: Job;
  sourceId: '_local_' | string;
  entryId: string;
  key: string;

  queryDefinition: QueryDefinition;
  normalizedParams: any;
  modelNums?: number[];
  copyAllCategories: boolean;
  transform?: Mat4;
}

interface JobEntryDefinition<Name extends QueryName> {
  sourceId?: string; // = '_local_',
  entryId: string;
  queryName: Name;
  queryParams: QueryParams<Name>;
  modelNums?: number[];
  copyAllCategories: boolean;
  transform?: Mat4;
}

export interface ResultWriterParams {
  encoding: Encoding;
  download: boolean;
  filename?: string;
  entryId?: string;
  queryName?: string;
}

export function JobEntry<Name extends QueryName>(definition: JobEntryDefinition<Name>): JobEntry {
  const queryDefinition = getQueryByName(definition.queryName);
  if (!queryDefinition) throw new HttpError(400, `Query '${definition.queryName}' is not supported.`);

  const normalizedParams = definition.queryParams;
  const sourceId = definition.sourceId || '_local_';

  return {
    job: void 0 as any,
    key: `${sourceId}/${definition.entryId}`,
    sourceId,
    entryId: definition.entryId,
    queryDefinition,
    normalizedParams,
    modelNums: definition.modelNums,
    copyAllCategories: !!definition.copyAllCategories,
    transform: definition.transform,
  };
}

export function createJob(definition: JobDefinition): Job {
  const job: Job = {
    id: UUID.create22(),
    datetime_utc: `${new Date().toISOString().replace(/T/, ' ').replace(/\..+/, '')}`,
    entries: definition.entries,
    writer: definition.writer,
    responseFormat: {
      tarball: !!definition?.options?.tarball,
      encoding: definition?.options?.encoding
        ? definition.options.encoding
        : !!(definition.options && definition.options.binary)
          ? 'bcif'
          : 'cif',
    },
    outputFilename: definition.options && definition.options.outputFilename,
  };
  definition.entries.forEach((e) => (e.job = job));
  return job;
}

/** True when nobody waits for the job's result anymore (the client disconnected or timed out). */
export function isJobCancelled(job: Job) {
  return !!job.writer.isCancelled?.();
}

/** The queue of the local API (`runLocal`): it sorts the jobs so that the jobs of an entry run one after another. */
class _JobQueue {
  private list: LinkedList<Job> = LinkedList();

  get size() {
    return this.list.count;
  }

  add(definition: JobDefinition) {
    const job = createJob(definition);
    this.list.addLast(job);
    return job.id;
  }

  hasNext(): boolean {
    return this.list.count > 0;
  }

  getNext(): Job {
    return this.list.removeFirst()!;
  }

  /** Sort the job list by key = sourceId/entryId */
  sort() {
    if (this.list.count === 0) return;

    const jobs: Job[] = [];
    for (let j = this.list.first; !!j; j = j.next) {
      jobs[jobs.length] = j.value;
    }

    jobs.sort((a, b) => (a.entries[0]?.key < b.entries[0]?.key ? -1 : 1));

    this.list = LinkedList();
    for (const j of jobs) {
      this.list.addLast(j);
    }
  }
}

export const JobManager = new _JobQueue();

export interface JobRunnerParams {
  /** Jobs processed at the same time. */
  maxConcurrentJobs: number;
  /** Jobs waiting for a slot; `add` refuses more. */
  maxQueueLength: number;
}

export interface JobRunnerStats {
  queued: number;
  running: number;
}

/**
 * Runs the web API's jobs: at most `maxConcurrentJobs` at a time, in the order they came, with at most
 * `maxQueueLength` waiting. A job whose client is gone when its turn comes is skipped.
 */
export class JobRunner {
  private queue: LinkedList<Job> = LinkedList();
  private running = 0;
  private accepting = true;
  private idle: (() => void)[] = [];

  get stats(): JobRunnerStats {
    return { queued: this.queue.count, running: this.running };
  }

  get isAccepting() {
    return this.accepting;
  }

  /** Queues the job; false when the queue is full or the runner no longer accepts jobs (`drain`). */
  add(job: Job) {
    if (!this.accepting || this.queue.count >= this.params().maxQueueLength) return false;
    this.queue.addLast(job);
    this.next();
    return true;
  }

  /** Stops accepting jobs; resolves when the queued and running ones are done. */
  drain(): Promise<void> {
    this.accepting = false;
    if (this.running === 0 && this.queue.count === 0) return Promise.resolve();
    return new Promise((resolve) => this.idle.push(resolve));
  }

  private next() {
    const { maxConcurrentJobs } = this.params();
    while (this.running < maxConcurrentJobs && this.queue.count > 0) {
      const job = this.queue.removeFirst()!;
      if (isJobCancelled(job)) {
        this.skip?.(job);
        continue;
      }
      this.running++;
      void this.run(job);
    }
    if (this.running === 0 && this.queue.count === 0) {
      for (const resolve of this.idle.splice(0)) resolve();
    }
  }

  private async run(job: Job) {
    try {
      await this.resolve(job);
    } catch {
      // `resolve` reports its own errors, the runner only must not stop
    } finally {
      this.running--;
      // A long job blocks the event loop: the hang-ups of clients that left meanwhile are only read in the next I/O
      // (poll) phase. Picked right away, or in this loop turn's `setImmediate`, their jobs would run for nobody; the
      // second `setImmediate` runs after that poll phase.
      setImmediate(() => setImmediate(() => this.next()));
    }
  }

  constructor(
    private resolve: (job: Job) => Promise<void>,
    private params: () => JobRunnerParams,
    private skip?: (job: Job) => void,
  ) {}
}
