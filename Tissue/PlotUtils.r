gseaScores <- getFromNamespace("gseaScores", "DOSE")

gsInfo <- function(object, geneSetID) {
    geneList <- object@geneList

    if (is.numeric(geneSetID))
        geneSetID <- object@result[geneSetID, "ID"]

    geneSet <- object@geneSets[[geneSetID]]
    exponent <- object@params[["exponent"]]
    df <- gseaScores(geneList, geneSet, exponent, fortify=TRUE)
    df$ymin <- 0
    df$ymax <- 0
    pos <- df$position == 1
    h <- diff(range(df$runningScore))/20
    df$ymin[pos] <- -h
    df$ymax[pos] <- h
    df$geneList <- geneList
    if (length(object@gene2Symbol) == 0) {
        df$gene <- names(geneList)
    } else {
        df$gene <- object@gene2Symbol[names(geneList)]
    }

    df$Description <- object@result[geneSetID, "Description"]
    return(df)
}

get_gsdata <- function(x, geneSetID) {
    if (length(geneSetID) == 1) {
        gsdata <- gsInfo(x, geneSetID)
        return(gsdata)
    } 
    
    lapply(geneSetID, gsInfo, object = x) |>
        yulab.utils::rbindlist()
}

abs.size <- function (page, orient) {
	size.dict <- list(
		A4 = c(210, 297),
		A3 = c(297, 420),
		A5 = c(148, 210),
		B3 = c(353, 500),
		B4 = c(250, 353),
		B5 = c(176, 250)
	)
	if (page %in% names(size.dict) == FALSE) { stop("Invalid page size.") }
	if (orient %in% c("h", "v") == FALSE) { stop("Invalid orientation.") }
	if (orient == "v") {
		return(size.dict[[page]])
	} else {
		return(rev(size.dict[[page]]))
	}
}

sized.ggsave <- function (plot,
						  filename,
						  page = "A4",
						  orient = "h",
						  rel.size = c(1, 1),
						  ratio = c(1, 1),
						  ppi = 300) {
	if (all(rel.size == 0)) { stop("Invalid relative size.") }
	if (any(rel.size < 0)) { stop("Invalid relative size.") }
	img.size <- abs.size(page, orient) * rel.size * ppi / 25.4
	if (img.size[1] == 0) { img.size[1] <- img.size[2] / ratio[2] * ratio[1] }
	if (img.size[2] == 0) { img.size[2] <- img.size[1] / ratio[1] * ratio[2] }
	ggsave(filename, plot, width = img.size[1], height = img.size[2], dpi = ppi, unit = "px")
}
