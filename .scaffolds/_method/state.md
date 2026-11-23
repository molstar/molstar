# State: <work id>   (updated <date time>, by <the orchestrator>)

<!-- The template for a work folder's state, `execution/<work-id>/state.md` (§state-file): the current state only,
     about a screen, written for a cold pickup. Rewritten in place after every merge, before every new lane and before
     stopping; never appended. What happened is in git, `checklist.md`'s ticks, `findings.md` and the record; decisions
     are in their design or the plan's dated lines. If a fresh agent cannot continue from it in a few minutes it is too
     thin; past 600 words the closed-work check notes it. Its last line, at the close, is "Closed <date>: <record
     file>"; the checks read it. Delete the comments once filled. -->

Branch: `<integration branch>`   Head: `<commit>`

## Where it stands

<!-- A few lines: what is merged of the checklist, what is half done and where. -->

## In progress

<!-- Every agent that runs now, with its start time, written before it starts: lanes, reviewers, measurements. A
     finished one leaves this list. -->

- <date time>: <agent>, <model>, <what>, worktree <path from the repository root>.

## Next

<!-- The next steps in order, enough to resume: which checklist item, which brief, which model, what it waits on. -->

1. <step>

## Blockers and open questions

<!-- Each with its drift note or gap and what would settle it. -->

## Verify

<!-- The commands that show where it stands, and when they were last green. -->
