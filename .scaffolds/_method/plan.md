# Plan: <work id>

<!-- The template for a work folder's plan, `execution/<work-id>/plan.md` (§plan-and-checklist): what and why, about a
     page or two, written once when the author agrees it. A re-plan adds a dated line at the end and never rewrites the
     plan. No status and no ticks here: the items and their ticks are in `checklist.md` beside it ([checklist.md](checklist.md)),
     where the current state is `state.md`'s. -->

<**Draft, <date>** | **Agreed, <date>** (the author: "<words>")>. The loop is [loop.md](loop.md); the items are
[checklist.md](checklist.md).

Work id `<work-id>`. <The milestone in one paragraph: what is built, from which designs, the author's words that set
it.>

## Where this is heading

<!-- One table: the earlier works (closed), this one, the next; so a reader sees the whole road. -->

| Work | Milestone | In one line |
|---|---|---|
| <work id> | <milestone> | <line> |

## The designs

<!-- Each agreed design this work builds and its lanes. The features check reads this table with the checklist's
     ticks: a design's row in the features table fails once every lane named here is ticked and the row still says
     planned. -->

| Design | Lanes |
|---|---|
| [<design>](../../design/<design>.md) | <LANE>, <LANE> |

## How this work runs

<!-- The loop's rules this plan leans on, and the fold-back it applies: one line each. -->

## Done when

<!-- What is true at the close, by feature, each checkable: a test, a view, a measurement. -->

## Gate 0

<!-- The decisions asked before stage 1, each a question with a recommendation. The answers, with the author's words,
     go into the designs or this plan's dated lines, never into the checklist. -->

## Testing strategy

<!-- What every lane tests, what each stage adds (properties, conformance suites, drawing against the data), and each
     design's performance representatives with their paths. -->

## Stages and lanes

<!-- Each stage in a few lines: its lanes (brief, model, what), the gate that ends it, the fix point after it. A stage
     with a new visible part starts with its mock lane beside the draft; the gate plays the mock. The checklist has one
     line for each lane, gate, fix point, the novice pass, the measurements and Finalize. -->

### Stage 0: designs

<The drafts, their mocks, their reviews; the gate where the author agrees them.>

### Stage <n>: <name>

<**<LANE>** (`briefs/<lane>.md`, <model>): what. The gate. The fix point after it, if any.>

### Last stage: use, measure, finalize

<The novice pass after the last lane that draws; the measurements at full size against their budgets; Finalize (the
loop's).>

## Measurements

<!-- Budgets, not wishes: each bundle's expected growth a stage, as a number from its designs ("not grown" only where
     nothing is added); each time and allocation target with its path and its test, measured at the stage that builds
     the path. -->

| What | Path and test | Budget | Stage |
|---|---|---|---|
| <bundle or frame> | <entry or scenario> | <+n kB, n ms, 0 B warm> | <n> |

## Order and parallelism

<!-- The tracks, what may run beside what, and what one lane owns at a time. -->

## Hot files (the orchestrator's only)

<!-- Files no lane stages: the lockfile, manifests, shared configs, the changelog folder, the work folder. -->

## Known risks

## Out of scope

## Briefs to write

<!-- The briefs the orchestrator writes before each lane starts. -->

<!-- Re-plans, one dated line each, appended here: "Re-planned <date>: <what changed> (the author: "<words>")." -->
