import type { StringLike } from './string-like.js';

/** Text ranges represented by start/end pairs. */
export interface Tokens {
    data: StringLike,
    count: number,
    indices: ArrayLike<number>
}
