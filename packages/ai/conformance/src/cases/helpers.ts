/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { expect } from 'vitest';
import {
  type ControlEvent,
  type Drive,
  type DriveHandlers,
  defaultDrive,
  drive,
  type HarnessEvent,
  type StepEvent,
  type TurnConfiguration,
} from '@molstar/ai-harness';
import { settle } from '../testing/clock.js';

/** Polls (letting the loop's promises settle) until the condition holds. */
export async function until(condition: () => boolean | Promise<boolean>, rounds = 500): Promise<void> {
  for (let i = 0; i < rounds; i += 1) {
    if (await condition()) return;
    await settle();
  }
  throw new Error('until: the condition never held');
}

/** A turn configuration for cases that need one of their own (no tools). */
export const TURN: TurnConfiguration = {
  model: { target: { provider: 'scripted', model: 'scripted-1' }, parameters: {} },
  tools: { enabled: [], preauthorized: [] },
  limits: {},
};

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
