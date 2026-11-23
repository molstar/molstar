# State: <work id>

<!-- The template for a work folder's state, `execution/<work-id>/state.md` (§state-file): written for a cold pickup,
     updated after every merge, before every new lane and before stopping. If a fresh agent cannot continue from it in a
     few minutes it is too thin. Its last line, at the close, is "Closed <date>: <record file>"; the closed-work check
     reads it. -->

## Decided

<!-- Decisions taken, newest first: the decision, the author's words and date, where it now lives (a design, a brief, the plan). -->

- <decision> (the author, <date>: "<words>"); <where>.

## Done

<!-- What is merged or finished, newest first, with its commit: enough to know what not to redo. -->

- <date>: <lane or stage> (`<commit>`).

## Next

<!-- The next steps in order, enough to resume: which brief, which model, what it waits on. -->

1. <step>

## In progress

<!-- Every agent that runs, with its start time, before it starts: lanes, reviewers, measurements. -->

- <date time>: <agent>, <model>, <what>.

## Blockers and open questions

<!-- Each with its drift note or gap and what would settle it. -->

## Verify

<!-- The commands that show where it stands, and when they were last green. -->

## Findings by lane

One row a merged lane: the model, the brief's items ticked, the review's and the look's findings by kind, who looked
and at what.

| Lane | Merged | Model | Items | Bug | Missing design | Look | Who looked | Notes |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
