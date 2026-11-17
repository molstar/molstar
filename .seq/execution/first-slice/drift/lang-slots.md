# Drift: the layer's `color` and `data` slots before the theme and query languages exist

Where: language-description.md ("a theme form", "an address or a query"). Lane B made `layer.color` an `any` param and
`layer.data` a `ref` only.

Resolution (recommended): switch `color` to a node slot of the theme sublanguage and `data` to a union of a ref and a
query slot when those land (lane C's forms; the theme forms with stage 4). Not blocking.
