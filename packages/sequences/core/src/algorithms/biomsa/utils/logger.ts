/**
 * @file logger
 * @author Paul Pillot <paul.pillot@libmol.org>
 * @license MIT
 * @copyright 2020
 */

// STORIES-UNTOLD(B4): the logger is a no-op: no console output, no timing (the
// original timed events with `performance.now()` and printed them with
// `console.table`; its calls are behind `DEBUG`, which is false).

function start() {}

function add(_name: string) {}

function summary() {}

export default {
  start,
  add,
  summary,
};
