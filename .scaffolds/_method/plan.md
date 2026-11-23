# Plan: <work id>

<!-- The template for a work folder's plan, `execution/<work-id>/plan.md` (§design-plan-implement). Its stages are the
     loop's items, ticked as they merge. Draft until the author agrees it; the agreement goes in its first line. -->

<**Draft, <date>** | **Agreed, <date>** (the author: "<words>")>. The loop is [loop.md](loop.md).

Work id `<work-id>`. <The milestone in one paragraph: what is built, from which designs, the author's words that set
it.>

## Where this is heading

<!-- One table: the earlier works (closed), this one, the next; so a reader sees the whole road. -->

| Work | Milestone | In one line |
|---|---|---|
| <work id> | <milestone> | <line> |

## The designs

<!-- Each agreed design this work builds and its lanes. The features check reads this table: a design's row in the
     features table fails once every lane named here is merged and the row still says planned. -->

| Design | Lanes |
|---|---|
| [<design>](../../design/<design>.md) | <LANE>, <LANE> |

## How this work runs

<!-- The loop's rules this plan leans on, and the fold-back it applies: one line each. -->

## Done when

<!-- What is true at the close, by feature, each checkable: a test, a view, a measurement. -->

## Gate 0

<!-- The decisions asked before stage 1, each a yes or no with a recommendation; the answers with the author's words. -->

## Testing strategy

<!-- What every lane tests, what each stage adds (properties, conformance suites, drawing against the data), and each
     design's performance representatives with their paths. -->

## Stages and lanes

<!-- One checkbox a lane or a gate: `- [ ] **<LANE>** (<brief>, <model>) <what>`; ticked with the merge commit. A stage
     with a new visible part starts with its mock lane beside the draft; the gate plays the mock. -->

### Stage 0: designs

- [ ] <drafts, their mocks, their reviews>
- [ ] GATE: the drafts agreed, with their mocks

### Fix points

- [ ] **FX<n>** after <stage>: `fixes.md`'s open items.

### Stage <n>: <name>

- [ ] **<LANE>** (`briefs/<lane>.md`, <model>) <what>.
- [ ] GATE

### Last stage: use, measure, finalize

- [ ] **NOVICE PASS** after the last lane that draws.
- [ ] **MEASURE**: <the representatives at full size, the bundles against their budgets>.
- [ ] **FINALIZE** (loop, Finalize).

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
