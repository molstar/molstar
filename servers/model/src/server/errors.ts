/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Sebastian Bittrich <sebastian.m.bittrich@gmail.com>
 */

/** An error the web API answers with `status`, any other error is a 500. */
export class HttpError extends Error {
  constructor(
    readonly status: number,
    message: string,
  ) {
    super(message);
  }
}

export function errorStatus(e: unknown) {
  return e instanceof HttpError ? e.status : 500;
}
