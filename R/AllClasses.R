
# GGPA2 class definition

setClass( Class="GGPA2",
  representation=representation(
    fit="list",
    summary="list",
	  setting="list",
    gwasPval="matrix",
    pgraph="matrix"
  )
)
