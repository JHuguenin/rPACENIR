#' rPACENIR : read NIR data from PACE
#'
#' import, pretraitement and tools for analyze data from NIR of UMR CEFE 5175
#'
#' @docType package
#' @name rPACENIR
#'
#' @import magrittr
#' @import stringr
#' @import rnirs
#' @import Rcpm
NULL

#' import dx files
#'
#' import et fait pleins d'autres choses
#'
#' @param file file csv.
#' @param wd working directory
#' @param export logical
#'
#' @return data.frame with data and or csv file
#' @export
#'
#' @examples
#' # blabla exemple
np.split.Multiblock.DX <- function (file = NULL, wd = NULL, export = FALSE){
  # intro ####
  if(is.null(wd) == TRUE) wd <- getwd()
  if((("samples" %in% dir(wd))==FALSE)&(export == TRUE)) dir.create("samples")
  if (is.null(file))  stop("You must provide a file name")

  # start import ####
  lines <- readLines(file)
  block_pat <- "##BLOCKS=(.*)"
  bc <- grep(block_pat, lines)
  bc <- sub(block_pat, "\\1", lines[bc])
  bc <- as.integer(bc)
  st_pat <- "##TITLE="
  st <- grep(st_pat, lines)
  st <- st[-1]
  end_pat <- "##END="
  end <- grep(end_pat, lines)
  end <- end[-length(end)]
  if (length(st) != length(end)) stop("Block starts and stops did not match")
  nb <- length(st)
  if (nb != bc) stop("Block count in link block did not match the number of blocks found")

  # block decomposition ####
  blocks <- vector("list", nb)
  for (i in 1:nb) blocks[[i]] <- (lines[st[i]:end[i]])

  fnames <- sapply(st,function(X,vec,sub_ch) gsub(sub_ch, "\\1", vec[X]), vec = lines, sub_ch = st_pat)
  fnames <- str_squish(fnames) %>% str_replace_all(" ","_") %>% str_replace_all("\\+","p") %>% str_replace_all("\\-","m")

  if (anyDuplicated(fnames))warning("Duplicated sample names found.\n\t\tLater samples will overwrite earlier samples\n\t\tunless you edit the original multiblock file.")
  names(blocks) <- paste0("S_",fnames)
  if(export == TRUE) for (i in 1:nb) writeLines(blocks[[i]], paste0(wd,"/samples/",fnames[i],".dx"))
  invisible(blocks)
}

#' read block
#'
#' lit les block du fichier et les transforment au bon format
#'
#' @param jdx file jdx
#'
#' @return a data.frame
#' @export
#'
#' @examples
#' # blabla exemple
np.read.block <- function (jdx = lnir[[1]]){
  sstt <- grep( "^\\s*##XYDATA\\s*=\\s*\\(X\\+\\+\\(Y\\.\\.Y\\)\\)$", jdx)
  send <- grep("^\\s*##END\\s*=", jdx)

  metadata <- jdx[1:(sstt[1]-1)]

  Format <- c("metadata", "XYY")
  FirstLine <- c(1, sstt)
  LastLine <- c(sstt[1] - 1, send)

  DF <- data.frame(Format, FirstLine, LastLine, stringsAsFactors = FALSE)
  keep_lines <- (1+DF[2,2]):(DF[2,3]-2)
  fmr <- sapply(jdx[keep_lines],str_split,pattern = " ", simplify = TRUE) %>% sapply(as.numeric) %>% t()
  rownames(fmr) <- NULL

  VL <- list()
  VL$Dataguide <- DF
  VL$Metadata <- metadata
  VL$spectra <- rowMeans(fmr[,-1]) %>% cbind(fmr[,1],.) %>% as.data.frame()
  colnames(VL$spectra) <- c("wl","int")

  fmr <- metadata[grep("YFACTOR",metadata)] %>% str_split(" ")
  VL$spectra$int <- as.numeric(fmr[[1]][2])*VL$spectra$int

  fmr <- (grep("CONCENTRATIONS",metadata)+1):(grep("DELTAX",metadata)-1)
  VL$concentration <- str_remove_all(metadata[fmr],"\\(") %>%
                        str_remove("\\)") %>% str_split(",",simplify = TRUE)
  fmr <- paste0(VL$concentration[,1],"_",VL$concentration[,3])
  VL$concentration <- as.numeric(VL$concentration[,2])
  names(VL$concentration) <- fmr
  return(VL)
}

#' Merge
#'
#' recupere les metadata et les data NIRS puis fait la correspondence
#'
#' @param spl csv file
#'
#' @return data.frame
#' @export
#'
#' @examples
#' # blabla exemple
np.export.date <- function(spl = lnir$S_115){
  ii <- c(grep("##TITLE",spl$Metadata), grep("##DATE",spl$Metadata), grep("##TIME",spl$Metadata))
  return(str_split(spl$Metadata[ii],"= ",simplify = TRUE)[,2])
}
# importer les metada du NIRS et matcher les deux =)

#' logical merge
#'
#' @param a logical
#' @param b  logical
#'
#' @return c logical
#' @export
#'
#' @examples
#' np.match(TRUE,FALSE)
np.match <- function(a,b){
  c <- merge(a,b)
  print(c)
  return(c)
}
