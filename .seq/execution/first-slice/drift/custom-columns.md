# Drift: typed columns outside declared column params

Where: format.md (columns declared in the description). A column in `custom`, in an `any` param or under an unknown
kind came back from JSON as a plain array, breaking the codec's round trip (lane A's review, blocker 2).

Resolved by the orchestrator (2026-11-17), to confirm with the author at the stage 1 gate: undeclared columns are tagged
in JSON (`{"$column": "<type>", "values": […]}`) and in the S-expression printing; declared column params stay plain
arrays. format.md is edited at the merge.
