# Drift: how a range wrapping a cyclic origin is written

Where: building-blocks.md and EC-9, EC-10 (a wrapped range is two runs; `start > stop` wraps only on a cyclic target).
Lane B's check reported such ranges as errors.

Resolved by the orchestrator (following EC-10): in a spec, `start > end` on a `cyclic` target means the range wraps the
origin; on a non-cyclic target it is an error. Internally it becomes two runs.
