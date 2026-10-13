/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { expect } from 'vitest';
import { hashMessage, type Message } from '@molstar/ai-agent';
import {
  buildMove,
  type ControlEvent,
  type Drive,
  type DriveHandlers,
  defaultDrive,
  drive,
  type HarnessEvent,
  type HarnessReader,
  type HeadId,
  type Move,
  type RunnerPort,
  type ScopeId,
  type StepEvent,
  type TurnConfiguration,
} from '@molstar/ai-harness';

/** A turn configuration for protocol-level cases (no tools). */
export const TURN: TurnConfiguration = {
  model: { target: { provider: 'scripted', model: 'scripted-1' }, parameters: {} },
  tools: { enabled: [], preauthorized: [] },
  limits: {},
};

export const hashes = (messages: readonly Message[]): Promise<string[]> => Promise.all(messages.map(hashMessage));

/** A move of `scope` from head `from` to that head's entries plus `messages`, for protocol-level cases. */
export async function appendMove(
  reader: HarnessReader,
  scope: ScopeId,
  from: HeadId,
  messages: readonly Message[],
): Promise<Move> {
  const [head] = await reader.heads([from]);
  const entries = [...(head?.entries ?? []), ...(await hashes(messages))];
  const move = await buildMove({ pointer: scope, edit: { entries }, read: () => [] });
  return { ...move, from };
}

/** A drive that answers the given kinds and defers the rest to the defaults. */
export const with_ = (handlers: DriveHandlers): Drive => drive(handlers, defaultDrive());

export const stepEvents = (events: readonly HarnessEvent[]): StepEvent[] =>
  events.filter((event): event is StepEvent => event.kind === 'step');

export const controls = (events: readonly HarnessEvent[], control?: ControlEvent['control']): ControlEvent[] =>
  events.filter(
    (event): event is ControlEvent => event.kind === 'control' && (control === undefined || event.control === control),
  );

/** The last step event, or a failure if there is none. */
export function lastStep(events: readonly HarnessEvent[]): StepEvent {
  const last = stepEvents(events).at(-1);
  if (last === undefined) throw new Error('no step event');
  return last;
}

/** A port that forwards to `port`, with some operations replaced. */
export const delegate = (port: RunnerPort, over: Partial<Omit<RunnerPort, 'id'>> = {}): RunnerPort => ({
  get id() {
    return port.id;
  },
  execution: () => port.execution(),
  headRecords: (ids) => port.headRecords(ids),
  segments: (ids) => port.segments(ids),
  lastCompaction: (scope) => port.lastCompaction(scope),
  heads: (ids) => port.heads(ids),
  pointers: (scopes) => port.pointers(scopes),
  entries: (ids) => port.entries(ids),
  events: (input) => port.events(input),
  steps: (ids) => port.steps(ids),
  turns: (input) => port.turns(input),
  blobs: (refs) => port.blobs(refs),
  watch: () => port.watch(),
  claim: (input) => port.claim(input),
  putBlobs: (lease, blobs) => port.putBlobs(lease, blobs),
  begin: (lease, input) => port.begin(lease, input),
  commit: (lease, batch) => port.commit(lease, batch),
  renew: (lease, input) => port.renew(lease, input),
  release: (lease, input) => port.release(lease, input),
  ...over,
});

type Expectation<T> = ReturnType<typeof expect<T>>;

/**
 * `expect(actual, message)`: Jest's `expect` takes no message, so a failure of any matcher on the
 * returned expectation (`.not`, `.resolves` and `.rejects` included) is rethrown with `message` in front.
 */
export function expectWith<T>(actual: T, message: string): Expectation<T> {
  const annotate = (error: unknown): unknown => {
    if (error instanceof Error) error.message = `${message}\n${error.message}`;
    return error;
  };
  const wrap = (target: object): any =>
    new Proxy(target, {
      get(object, key) {
        const value = Reflect.get(object, key);
        if (typeof value === 'function') {
          return (...args: unknown[]) => {
            try {
              const result = value.apply(object, args);
              return result instanceof Promise
                ? result.catch((error: unknown) => {
                    throw annotate(error);
                  })
                : result;
            } catch (error) {
              throw annotate(error);
            }
          };
        }
        return typeof value === 'object' && value !== null ? wrap(value) : value;
      },
    });
  return wrap(expect(actual));
}
