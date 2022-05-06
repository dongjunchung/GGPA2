
# generic methods for "GGPA2" class

setGeneric( "fdr",
    function( object, ... )
    standardGeneric("fdr")
)

setGeneric( "assoc",
    function( object, ... )
    standardGeneric("assoc")
)

setGeneric( "estimates",
    function( object, ... )
    standardGeneric("estimates")
)

setGeneric(
  "get_fit",
  function(x) standardGeneric("get_fit")
)

setGeneric(
  "get_summary",
  function(x) standardGeneric("get_summary")
)

setGeneric(
  "get_setting",
  function(x) standardGeneric("get_setting")
)

setGeneric(
  "get_gwasPval",
  function(x) standardGeneric("get_gwasPval")
)

setGeneric(
  "get_pgraph",
  function(x) standardGeneric("get_pgraph")
)

# get the annotation matrix
setGeneric(
  "get_annotMat",
  function(x) standardGeneric("get_annotMat")
)

