/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { HarnessObservable, LaunchKey, LaunchReport, Supervisor } from '@molstar/ai-harness';

/**
 * A simulated supervisor for the liveness model (S1): `ensure` starts `start(launch)` once per key,
 * unless `failLaunch` says this launch is lost (as a dead launcher would lose it). The report is a
 * plain replay of what happened, enough for the property model.
 */
export function simulatedSupervisor(options: {
  start(launch: LaunchKey): Promise<unknown>;
  failLaunch?(launch: LaunchKey, attempt: number): boolean;
}): Supervisor & { readonly started: LaunchKey[]; readonly lost: LaunchKey[]; settled(): Promise<void> } {
  const started: LaunchKey[] = [];
  const lost: LaunchKey[] = [];
  const running: Promise<unknown>[] = [];
  const seen = new Set<LaunchKey>();
  return {
    started,
    lost,
    async settled() {
      await Promise.allSettled(running);
    },
    ensure(_harnessId, launch): HarnessObservable<LaunchReport> {
      if (!seen.has(launch)) {
        seen.add(launch);
        if (options.failLaunch?.(launch, seen.size) === true) lost.push(launch);
        else {
          started.push(launch);
          running.push(options.start(launch).catch(() => undefined));
        }
      }
      return {
        subscribe(onNext, _onError, onComplete) {
          onNext({ phase: 'launching', key: launch });
          onComplete?.();
          return () => undefined;
        },
      };
    },
  };
}
