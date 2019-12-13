(TeX-add-style-hook
 "imac2020_intred_extabs"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("article" "10pt" "print")))
   (TeX-add-to-alist 'LaTeX-provided-package-options
                     '(("geometry" "top=0.75in" "left=0.75in" "right=0.75in" "bottom=1in") ("natbib" "numbers" "sort&compress")))
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperref")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperimage")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperbaseurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "nolinkurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "url")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "path")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "path")
   (TeX-run-style-hooks
    "latex2e"
    "article"
    "art10"
    "amsmath"
    "amssymb"
    "authblk"
    "geometry"
    "times"
    "graphicx"
    "subcaption"
    "hyperref"
    "cleveref"
    "xcolor"
    "natbib")
   (TeX-add-symbols
    '("utilde" 1)
    '("ned" 1))
   (LaTeX-add-labels
    "sec:intr-motiv"
    "sec:descr-appr"
    "fig:iraref"
    "sec:whole-joint-appr"
    "fig:wja"
    "sec:interf-remesh-appr"
    "eq:interp"
    "eq:dynsys"
    "eq:ROMsys"
    "eq:hyperred"
    "fig:ira"
    "sec:results"
    "fig:res"
    "sec:disc-concl")
   (LaTeX-add-bibliographies
    "Refs"))
 :latex)

