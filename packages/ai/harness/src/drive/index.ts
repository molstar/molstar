/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type {
  Answer,
  CallView,
  Decision,
  DecisionKind,
  Drive,
  DriveContext,
  DriveHandlers,
  SkipAnswer,
} from '../contracts/drive.js';

export interface DefaultDriveOptions {
  readonly retry?: { readonly max: number; backoffMs(attempt: number): number };
  readonly maxConcurrency?: number;
}

const isObject = (value: unknown): value is Record<string, unknown> => typeof value === 'object' && value !== null;
const has = (value: unknown, key: string) => isObject(value) && key in value;

/** Whether `answer` is legal for `decision` (harness.md §7); `stop` and `fail` always are. */
export function legal(decision: Decision, answer: unknown): boolean {
  if (answer === 'stop' || has(answer, 'fail')) return true;
  switch (decision.kind) {
    case 'permission':
      return answer === 'allow' || has(answer, 'deny') || has(answer, 'park');
    case 'schedule': {
      if (!isObject(answer) || !Array.isArray(answer.stages)) return false;
      const calls = (decision as Decision<'schedule'>).payload.calls.map((call) => call.callId).sort();
      const listed = (answer.stages as unknown[][]).flat().map(String).sort();
      return listed.length === calls.length && listed.every((id, i) => id === calls[i]);
    }
    case 'budget':
      return has(answer, 'extend');
    case 'failure': {
      const payload = (decision as Decision<'failure'>).payload;
      if (has(answer, 'retry')) return true;
      if (has(answer, 'replace')) return payload.call === undefined || payload.call.toolId === '';
      if (has(answer, 'report'))
        return (
          payload.call !== undefined &&
          payload.call.toolId !== '' &&
          skipsUnstarted((answer as { report: { skip?: SkipAnswer } }).report.skip, payload.unstarted ?? [])
        );
      if (has(answer, 'compact')) return payload.attempt === undefined;
      return false;
    }
    case 'model-result':
      return answer === 'accept' || has(answer, 'replace');
    case 'repair':
      return (
        has(answer, 'report') &&
        skipsUnstarted(
          (answer as { report: { skip?: SkipAnswer } }).report.skip,
          (decision as Decision<'repair'>).payload.unstarted,
        )
      );
    case 'boundary': {
      const completing = (decision as Decision<'boundary'>).payload.completing;
      return answer === 'proceed' || has(answer, 'steer') || has(answer, 'compact') || (completing && answer === 'end');
    }
    case 'interrupt':
      return (
        answer === 'abandon' || (answer === 'skip' && (decision as Decision<'interrupt'>).payload.stepKind === 'tools')
      );
    case 'next-turn':
      return answer === 'release' || has(answer, 'begin');
  }
}

/** A skip names only calls that have not started (or none: all of them). */
function skipsUnstarted(skip: SkipAnswer | undefined, unstarted: readonly CallView[]): boolean {
  if (skip?.calls === undefined) return true;
  const ids = new Set(unstarted.map((call) => call.callId));
  return skip.calls.every((id) => ids.has(id));
}

/** A drive that answers some kinds and defers the rest, and its handlers' `pass`, to `fallback`. */
export function drive(handlers: DriveHandlers, fallback?: Drive): Drive {
  return {
    async decide<K extends DecisionKind>(decision: Decision<K>, context: DriveContext) {
      const handler = handlers[decision.kind] as
        | ((d: Decision<K>, c: DriveContext) => Answer<K> | 'pass' | Promise<Answer<K> | 'pass'>)
        | undefined;
      const answer = handler !== undefined ? await handler(decision, context) : 'pass';
      return answer === 'pass' && fallback !== undefined ? await fallback.decide(decision, context) : answer;
    },
  };
}

const backoff = (attempt: number) => Math.min(30_000, 500 * 2 ** (attempt - 1));

/** The defaults for every kind (harness.md §7). */
export function defaultDrive(options: DefaultDriveOptions = {}): Drive {
  const retry = options.retry ?? { max: 3, backoffMs: backoff };
  return drive({
    permission: () => 'allow',
    schedule: ({ payload }) => ({
      stages: payload.calls.map((call) => [call.callId]),
      ...(options.maxConcurrency !== undefined ? { maxConcurrency: options.maxConcurrency } : {}),
    }),
    budget: () => ({ fail: { reason: 'budget' } }),
    failure: ({ payload }) => {
      const tool = payload.call !== undefined && payload.call.toolId !== '';
      const attempt = payload.attempt ?? 0;
      // Not before the provider said it would take the call again (`retryAfterMs`).
      if (payload.attempt !== undefined && payload.failure.transient && attempt < retry.max)
        return { retry: { afterMs: Math.max(retry.backoffMs(attempt), payload.failure.retryAfterMs ?? 0) } };
      return tool ? { report: {} } : { fail: { reason: payload.failure.code, message: payload.failure.message } };
    },
    'model-result': () => 'accept',
    repair: () => ({ report: {} }),
    boundary: ({ payload }, context) =>
      payload.completing && context.inbox.length > 0 ? { steer: context.inbox.map((item) => item.key) } : 'proceed',
    interrupt: ({ payload }) => (payload.stepKind === 'tools' ? 'skip' : 'abandon'),
    'next-turn': (_decision, context) =>
      context.inbox.length > 0 ? { begin: { keys: context.inbox.map((item) => item.key) } } : 'release',
  });
}
