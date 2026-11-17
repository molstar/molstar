# Drift: positional arguments in query forms

Where: locations-and-queries.md, "The first forms" (`(extent <elements>)`, `(and …)`). repr's address grammar has keyword
params only; lane C built `(extent :of …)`, `(and :of [ … ])`, `(or :of [ … ])`, `(not :of …)`.

Options: (1) keep `:of` and update the design (recommended); (2) positional data forms in the address reader.
Also: add `seq.query.*` to the issue namespaces in language-description.md.
