DEFAULT_LAS_DATA = c(
  "X",
  "Y",
  "Z",
  "gpstime",
  "ReturnNumber",
  "NumberOfReturns",
  "ScanDirectionFlag",
  "EdgeOfFlightlin",
  "Classificati",
  "Synthetic_flag",
  "Withheld_flag",
  "Keypoint_flag",
  "ScanAngle",
  "UserData",
  "PointSourceID",
  "R" ,
  "G" ,
  "B",
  "NIR"
)
#' Merge neighbors variables by a function
#'
#' @return A numeric vector with the neighbors aggregated value for a given variable
#'
#' @param las an object of class LAS to be the reference
#' @param las_from an object of class LAS from which neighbors will be aggregated
#' @param k interger. The number of k-nearest neighbors.
#' @param func character.A function name to aggregate the neighbors values, accepted values
#' ("mean", "median", "q1", "q3", "var", "sd", "max", "min") or an R function (slower), default "mean"
#' @param var string. The name of a point variable to aggregate data from.
#' @param name string. The name variable to store the result.
#' @param name string. Description for the variable.
#'
#'@export
lasNN = function(las, las_from, k = 1, func = "mean", var = "Intensity", name="V", desc="Description for V")
{
  stopifnot(
    nrow(las@data) > k,
    var %in% names(las@data)
  )

  # Check if func is an R function or a string
  if (is.function(func))
  {
    stopifnot(
      is.numeric(func(c(1,2,3))[1])
    )
    out <- C_lasmerge_neighbors_Rfunc(las, las_from, k, func, var)
  } else {
    out <- C_lasmerge_neighbors(las, las_from, k, func, var)
  }


  # If is a default name then use default
  if (name %in% DEFAULT_LAS_DATA) {
    las@data[[name]] = out
  } else {
    las = lidR::lasadddata(las, out, name)
    las = lidR::lasaddextrabytes(las, name = name, desc=desc)
  }

  return(las)
}
