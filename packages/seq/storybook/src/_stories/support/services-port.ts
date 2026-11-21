/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { DefaultAssetManager } from '@molstar/core/util/assets';
import {
  type AssetPolicy,
  AssetRefused,
  BrowserAssetPolicy,
  policyAssets,
  type ReaderAssetsPort,
} from '@molstar/seqx-core/context';
import type { AssetReader } from '@molstar/seqx-core/data';
import { Recorded, serviceOf, ServiceHosts } from '@molstar/seqx-examples/services/endpoints';

// The assets port of an application composing services (examples.md, Showcases): a spec names each service's own location
// and this port answers it, in one of two modes the application switches at run time. `fixtures`: the recorded file of that
// location (a location with none is `NotRecorded`, said so), after a delay on the port's clock (a service's latency: the
// stories and their twins share it, so a twin runs on a fake clock). `live`: the service itself, through the asset policy
// (`https:` and the services' hosts only, no credentials, redirects refused, the browser's size and time bounds): a host
// the policy does not list is refused (`AssetRefused`) whatever the spec says. A location that is no service's (a relative
// one, a dropped file's) is answered by the fixtures reader in either mode and refused by the policy otherwise.

/** The mode of the port: recorded answers, or the services themselves. */
export type ServiceMode = 'fixtures' | 'live';

/** A location of a service with no recorded answer: the fixtures cover the showcase's entries only. */
export class NotRecorded extends Error {
  constructor(readonly uri: string) {
    super(`No recorded answer for ${uri}: switch to the live services to read it`);
    this.name = 'NotRecorded';
  }
}

/** The timers the port's delay runs on. */
export interface PortClock {
  setTimeout(fn: () => void, ms: number): unknown;
  clearTimeout(handle: unknown): void;
}

/** One read: what the application's log shows. */
export interface ServiceRead {
  readonly uri: string;
  readonly mode: ServiceMode;
  /** The service's name, or none for a location that is no service's. */
  readonly service?: string;
  /** The recorded file that answered (fixtures mode). */
  readonly file?: string;
}

export interface ServicesPortOptions {
  /** The mode now (read at every read: the switch takes effect on the next one). */
  readonly mode: () => ServiceMode;
  /** Reads a recorded file by its location relative to the examples (the page fetches it, a twin reads it from disk). */
  readonly fixture: (file: string) => Promise<Uint8Array>;
  /** Milliseconds a recorded answer waits (a service's latency); 0. */
  readonly delayMs?: number;
  /** The clock of the delay: the page's timers by default, a twin's fake. */
  readonly clock?: PortClock;
  /** The live reader (a twin's refuses the network); by default the manager behind `LivePolicy`. */
  readonly live?: AssetReader;
  /** Told of each read as it starts. */
  readonly onRead?: (read: ServiceRead) => void;
}

/** What the live mode may reach: `https:` and the services' hosts, no relative location. */
export const LivePolicy: AssetPolicy = {
  ...BrowserAssetPolicy,
  relative: false,
  hosts: ServiceHosts,
};

const realClock: PortClock = {
  setTimeout: (fn, ms) => setTimeout(fn, ms),
  clearTimeout: (h) => clearTimeout(h as ReturnType<typeof setTimeout>),
};

/** The port: a reader answering by the mode, and the reads it has made (newest last). */
export function servicesPort(
  options: ServicesPortOptions,
): ReaderAssetsPort & { readonly reads: readonly ServiceRead[] } {
  const reads: ServiceRead[] = [];
  const clock = options.clock ?? realClock;
  let live = options.live;
  let manager: DefaultAssetManager | undefined;
  const liveReader = (): AssetReader => {
    if (!live) {
      manager = new DefaultAssetManager();
      live = policyAssets(manager, LivePolicy);
    }
    return live;
  };
  const wait = (ms: number, signal: AbortSignal) =>
    new Promise<void>((resolve, reject) => {
      if (ms <= 0) return resolve();
      if (signal.aborted) return reject(signal.reason ?? new Error('aborted'));
      const handle = clock.setTimeout(() => {
        signal.removeEventListener('abort', onAbort);
        resolve();
      }, ms);
      const onAbort = () => {
        clock.clearTimeout(handle);
        reject(signal.reason ?? new Error('aborted'));
      };
      signal.addEventListener('abort', onAbort, { once: true });
    });

  /** What answers the location now: a recorded file, or the live reader. */
  async function answer<T>(
    uri: string,
    signal: AbortSignal,
    fromFile: (file: string) => Promise<T>,
    fromLive: (r: AssetReader) => Promise<T>,
  ): Promise<T> {
    const mode = options.mode();
    const service = serviceOf(uri);
    const file = Recorded.get(uri);
    const read: ServiceRead = {
      uri,
      mode,
      ...(service ? { service: service.name } : {}),
      ...(mode === 'fixtures' && file ? { file } : {}),
    };
    reads.push(read);
    options.onRead?.(read);
    if (mode === 'live') return fromLive(liveReader());
    if (service) {
      if (!file) throw new NotRecorded(uri);
      await wait(options.delayMs ?? 0, signal);
      return fromFile(file);
    }
    // not a service's: a relative location of the examples (a recorded file by name), nothing absolute
    if (/^[a-z][a-z0-9+.-]*:/i.test(uri)) throw new AssetRefused(uri, 'not a service of the showcase');
    return fromFile(uri);
  }

  const text = new TextDecoder();
  return {
    reads,
    reader: {
      readText: (uri, signal) =>
        answer(
          uri,
          signal,
          async (f) => text.decode(await options.fixture(f)),
          (r) => r.readText(uri, signal),
        ),
      readBytes: (uri, signal) =>
        answer(
          uri,
          signal,
          (f) => options.fixture(f),
          (r) => r.readBytes(uri, signal),
        ),
    },
    dispose() {
      (manager as { dispose?: () => void } | undefined)?.dispose?.();
    },
  };
}
