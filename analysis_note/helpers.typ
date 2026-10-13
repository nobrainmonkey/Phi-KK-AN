// Layout helpers shared by the section files (#import "../helpers.typ": ...).

// Two figures side by side, floated together so they stay next to each other.
// Each keeps its own number, caption and label:
//   #figure-pair(
//     [#figure(image(...), caption: [...]) <fig:a>],
//     [#figure(image(...), caption: [...]) <fig:b>],
//   )
#let figure-pair(left, right) = place(auto, float: true, grid(
  columns: (1fr, 1fr),
  column-gutter: 1em,
  align: top,
  left, right,
))
