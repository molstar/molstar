/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Json } from '@molstar/ai-agent';
import type {
  Clock,
  Execution,
  HarnessClient,
  HarnessId,
  HarnessObservable,
  LaunchKey,
  LaunchReport,
  Runner,
  Supervisor,
  SupervisorPolicy,
} from '@molstar/ai-harness';

export interface InProcessSupervisorOptions extends SupervisorPolicy {
  readonly client: HarnessClient;
  startRunner(launch: LaunchKey): Runner;
  readonly clock?: Clock;
  /** How long after a launch window or a lease ends to `check` (default 1 s). */
  readonly marginMs?: number;
}

const systemClock: Clock = {
  now: () => Date.now(),
  sleep: (ms, signal) =>
    new Promise<void>((resolve, reject) => {
      if (signal?.aborted) return reject(signal.reason);
      const timer = setTimeout(resolve, ms);
      signal?.addEventListener(
        'abort',
        () => {
          clearTimeout(timer);
          reject(signal.reason);
        },
        { once: true },
      );
    }),
};

/** One launch's report: replayed to late subscribers, completed after `picked-up` or `failed`. */
class Report implements HarnessObservable<LaunchReport> {
  readonly history: LaunchReport[] = [];
  done = false;
  /** Events after this seq may hold its claim. */
  since = 0;
  private readonly subscribers = new Set<{ next(r: LaunchReport): void; complete?(): void }>();

  subscribe(onNext: (r: LaunchReport) => void, _onError?: (e: unknown) => void, onComplete?: () => void) {
    for (const report of this.history) onNext(report);
    if (this.done) {
      onComplete?.();
      return () => undefined;
    }
    const subscriber = { next: onNext, ...(onComplete !== undefined ? { complete: onComplete } : {}) };
    this.subscribers.add(subscriber);
    return () => {
      this.subscribers.delete(subscriber);
    };
  }

  push(report: LaunchReport, last = false): void {
    if (this.done) return;
    this.history.push(report);
    this.done = last;
    for (const subscriber of [...this.subscribers]) {
      try {
        subscriber.next(report);
        if (last) subscriber.complete?.();
      } catch {
        // A subscriber cannot affect the launch.
      }
    }
    if (last) this.subscribers.clear();
  }
}

/**
 * Liveness in one process (harness.md §4), for a tab or a Node worker: starts a runner per launch key,
 * starts any launch it sees for its harness, checks when a launch window or lease ends, and reports
 * stalls. One supervisor per harness client.
 */
export function createInProcessSupervisor(options: InProcessSupervisorOptions): Supervisor & { close(): void } {
  const { client } = options;
  const clock = options.clock ?? systemClock;
  const margin = options.marginMs ?? 1_000;
  const reports = new Map<LaunchKey, Report>();
  const closed = new AbortController();
  let wake: AbortController | undefined;
  let held: LaunchKey | undefined;
  let stalledSince: number | undefined;

  const start = (launch: LaunchKey, report: Report, execution: Execution) => {
    const waiting = execution.inbox.map((item) => ({
      key: item.key,
      messages: item.messages,
      ...(item.options !== undefined ? { options: item.options } : {}),
      submittedAt: item.submittedAt,
      interrupt: item.interrupt !== undefined,
    }));
    if (execution.turn === undefined && waiting.length > 0 && options.shouldStart?.(waiting) === false) {
      held = launch;
      return report.push({ phase: 'failed', error: { reason: 'held' }, relaunch: false }, true);
    }
    let runner: Runner;
    try {
      runner = options.startRunner(launch);
    } catch (error) {
      return report.push({ phase: 'failed', error: String(error), relaunch: true }, true);
    }
    runner.run(launch).then(
      (exit) => {
        if (exit.reason === 'stale-launch' || exit.reason === 'stale-schema')
          report.push({ phase: 'failed', error: { reason: exit.reason }, relaunch: false }, true);
        wake?.abort();
      },
      (error: unknown) => {
        report.push({ phase: 'failed', error: String(error) as Json, relaunch: true }, true);
        wake?.abort();
      },
    );
  };

  /** Reads the Execution: starts unseen launches, reports pick-ups, and says when to look again. */
  async function survey(): Promise<number | undefined> {
    const execution = await client.execution();
    const runner = execution.runner;
    const now = clock.now();
    await pickUps(execution);
    if (runner.phase === 'launching') {
      if (!reports.has(runner.key) && held !== runner.key) ensure(client.id, runner.key);
      return runner.until + margin;
    }
    if (runner.phase === 'none') return undefined;
    let next = runner.leaseUntil + margin;
    const after = runner.progress !== undefined ? options.stalledAfterMs?.(runner.progress.phase) : undefined;
    if (runner.progress !== undefined && after !== undefined) {
      const at = runner.progress.since + after;
      if (now >= at && stalledSince !== runner.progress.since) {
        stalledSince = runner.progress.since;
        options.onStalled?.(runner);
      } else if (now < at) next = Math.min(next, at);
    }
    return next;
  }

  /** Reports the launches whose `claimed` event has appeared. */
  async function pickUps(execution: Execution): Promise<void> {
    const pending = [...reports.values()].filter((report) => !report.done);
    if (pending.length === 0) return;
    const after = Math.min(...pending.map((report) => report.since));
    const events = await client.events({ after, limit: 1_000 });
    for (const event of events) {
      if (event.kind !== 'control' || event.control !== 'claimed') continue;
      const detail = event.detail as { launch?: string; runnerId?: string };
      reports.get(detail.launch ?? '')?.push(
        {
          phase: 'picked-up',
          runnerId: detail.runnerId ?? '',
          ...(execution.turn !== undefined ? { turnId: execution.turn.id } : {}),
        },
        true,
      );
    }
    const read = events.at(-1)?.seq ?? after;
    for (const report of pending) report.since = Math.max(report.since, read);
  }

  async function watch(): Promise<void> {
    const unsubscribe = client.watch().subscribe(
      () => wake?.abort(),
      () => wake?.abort(),
    );
    closed.signal.addEventListener('abort', unsubscribe, { once: true });
    while (!closed.signal.aborted) {
      wake = new AbortController();
      const signal = AbortSignal.any([wake.signal, closed.signal]);
      try {
        const next = await survey();
        if (next !== undefined && clock.now() >= next) {
          const { launch } = await client.check();
          if (launch !== undefined) {
            held = undefined;
            ensure(client.id, launch);
          }
          continue;
        }
        if (next !== undefined) await clock.sleep(next - clock.now(), signal);
        // A wake during the read has already aborted the signal: its listener would never fire.
        else
          await new Promise((_, reject) => {
            if (signal.aborted) return reject(signal.reason);
            signal.addEventListener('abort', reject, { once: true });
          });
      } catch {
        // Woken (a change, a runner's end) or a failed read: look again.
        if (!wake.signal.aborted && !closed.signal.aborted)
          await clock.sleep(margin, closed.signal).catch(() => undefined);
      }
    }
  }

  function ensure(harnessId: HarnessId, launch: LaunchKey): HarnessObservable<LaunchReport> {
    if (harnessId !== client.id) throw new Error(`this supervisor serves ${client.id}, not ${harnessId}`);
    let report = reports.get(launch);
    if (report !== undefined) return report;
    const made = new Report();
    report = made;
    reports.set(launch, made);
    made.push({ phase: 'launching', key: launch });
    client.execution().then(
      (execution) => {
        made.since = Math.max(0, execution.seq - 1);
        start(launch, made, execution);
      },
      (error: unknown) => made.push({ phase: 'failed', error: String(error), relaunch: true }, true),
    );
    wake?.abort();
    return made;
  }

  // Watches from the start: a lease that runs out (an owner that died) is checked and relaunched
  // even when nobody asked this supervisor to launch anything.
  void watch();
  return { ensure, close: () => closed.abort() };
}
