#' Run LDAK executable for Linux, Mac, and Windows systems.
#'
#' This function automatically detects the operating system and runs
#' the corresponding LDAK binary for the input arguments.
#'
#' @param args Character vector of arguments to pass to LDAK
#' @return The output of the system call
#' @export
ldak <- function(args = character()) {
  os <- Sys.info()[["sysname"]]

  exe_name <- switch(
    os,
    "Linux"   = "ldak6.1.linux",
    "Darwin"  = "ldak6.1.mac",
    "Windows" = "ldak6.1.exe",
    stop("Unsupported OS: ", os)
  )

  exe_path <- system.file("bin", exe_name, package = "LDAK")

  if (exe_path == "") {
    stop("Could not find the LDAK binary for your operating system.")
  }

  # Ensure executable permissions on Unix-like systems
  if (os != "Windows") {
    Sys.chmod(exe_path, mode = "755")
  }

  # Run the executable
  exit_code <- system2(exe_path, args = args, stdout = "", stderr = "")

  if (!exit_code %in% c(0,1)) {
    stop("LDAK execution failed with exit code ", exit_code)
  }

  invisible(exit_code)
}
