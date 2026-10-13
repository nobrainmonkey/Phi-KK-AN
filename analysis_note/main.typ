// UPC phi meson analysis note -- entry point.
//
//   typst compile main.typ     // -> main.pdf
//   typst watch main.typ       // rebuild on every save
//
// Each \section of the original LaTeX note lives in its own file under
// sections/. To add a section, create a file there and #include it below.
// Shared layout helpers (e.g. figure-pair) live in helpers.typ.
//
// Math symbol note: LaTeX's \phi (straight phi) is `phi.alt` in Typst; plain
// Typst `phi` is the looped variant, i.e. LaTeX's \varphi.

#set document(title: [UPC $phi.alt$ Meson Analysis Note], author: "Xihe Han")

// ---- Page and text (close to the LaTeX `article` look) ---------------------
#set page(paper: "us-letter", margin: 1.25in, numbering: "1")
#set text(font: "New Computer Modern", size: 11pt, lang: "en")
#set par(justify: true, first-line-indent: 1.5em, spacing: 0.65em)

// ---- Math, headings, lists, tables, links ----------------------------------
// Typst defaults to the heavier "Book" math font, which looks bolder than the
// body text; pin the Regular weight so inline math blends in.
#show math.equation: set text(weight: 400)
#set math.equation(numbering: "(1)") // number display equations, as in LaTeX
#set heading(numbering: "1.1")
#show list: set block(above: 1em, below: 1em)
#show figure: set block(above: 1.4em, below: 1.4em)
#show figure.where(kind: table): set figure.caption(position: top)
#show link: set text(fill: rgb("#1a5fb4"))
#show ref: set text(fill: rgb("#1a5fb4"))
#set par.line(numbering: "1")

// ---- Title block (plays the role of \maketitle) ----------------------------
#align(center)[
  #title()
  #v(0.2em)
  #text(size: 1.2em)[Xihe Han]
  #v(0.2em)
  #datetime.today().display("[month repr:long] [day], [year]")
]
#v(1em)

// ---- Table of contents (plays the role of \tableofcontents) ----------------
// Built from the numbered headings; delete these two lines to remove it.
#show outline: set par.line(numbering: none)
#outline(depth: 2)
#pagebreak()

// Stray line from the original .tex: it sat above \section{Introduction},
// after a lone "S". Kept here as a comment (not typeset) in case it is
// meant for a later section.
//
//   uls_data->Draw("uls_phi_mass>>h(20,1,1.2)","uls_phi_pt<0.1 && abs(uls_kaon1_nsigma_kaon-1.6)<3 && abs(uls_kaon2_nsigma_kaon-1.6)<3 && uls_pair_BBC_Veto_East==0 && uls_pair_BBC_Veto_West==0 && abs(uls_phi_rapidity)<0.25")

// ---- Sections --------------------------------------------------------------
#include "sections/01_introduction.typ"
#include "sections/02_dataset_and_trigger.typ"
#include "sections/03_event_selection.typ"
#include "sections/04_Signal_extraction.typ"
#include "sections/05_efficiencies.typ"
