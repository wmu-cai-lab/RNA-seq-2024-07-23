my.ridgeplot <- function(x, show.id, description = show.id, fill = "p.adjust", pos.lab.shift = NA) {
	if (is.list(x)) {
		for (y in 1:length(x)) {
			if (x[[y]]@readable && length(x[[y]]@gene2Symbol) > 0) {
				id <- match(names(x[[y]]@geneList), names(x[[y]]@gene2Symbol))
				names(x[[y]]@geneList) <- x[[y]]@gene2Symbol[id]
			} 
		}
		if (!is.list(show.id) || !is.list(description)) { stop("Inconsistent length of x, show.id and description.") }
		len <- length(x)
		if (len != length(show.id) || len != length(description)) { stop("Inconsistent length of x, show.id and description.") }
		gs2val <- lapply(1:len, function(y) {
			gs2id <- geneInCategory(x[[y]])[show.id[[y]]]
			lapply(gs2id, function(id) {
				res <- x[[y]]@geneList[id]
				res <- res[!is.na(res)]
			}) %>% return()
		})
		gs2val <- do.call(append, gs2val)
		description <- unlist(description)
		color <- lapply(1:len, function(y) {
			return(x[[y]][show.id[[y]], fill])
		})
		color <- do.call(c, color)
		NES <- lapply(1:len, function(y) {
			return(x[[y]][show.id[[y]], "NES"])
		})
		NES <- do.call(c, NES)
	} else {
		if (x@readable && length(x@gene2Symbol) > 0) {
			id <- match(names(x@geneList), names(x@gene2Symbol))
			names(x@geneList) <- x@gene2Symbol[id]
		} 
    	gs2id <- geneInCategory(x)[show.id]
		gs2val <- lapply(gs2id, function(id) {
			res <- x@geneList[id]
			res <- res[!is.na(res)]
		})
		color <- x[show.id, fill]
		NES <- x[show.id, "NES"]
	}

    len <- sapply(gs2val, length)
    gs2val.df <- data.frame(category = rep(description, times = len),
                            color = rep(color, times = len),
                            value = unlist(gs2val))
    gs2val.df$category <- factor(gs2val.df$category, levels = description)
	NES.df <- data.frame(category = description,
						 NES = paste0("NES: ", sprintf("%.2f", NES)),
						 NES.val = NES)
		
	if (!is.na(pos.lab.shift)) {
		ggplot(gs2val.df, aes(x = value, y = category, fill = color)) +
			ggridges::geom_density_ridges() +
			geom_text(data = NES.df,
				mapping = aes(x = ifelse(NES.val < 0, pos.lab.shift, floor(min(gs2val.df$value))), y = category, label = category),
				nudge_y = 0.6, hjust = 0)
	} else {
		ggplot(gs2val.df, aes(x = value, y = category, fill = color)) +
			ggridges::geom_density_ridges() +
			geom_text(data = NES.df, mapping = aes(x = floor(min(gs2val.df$value)), y = category, label = category), nudge_y = 0.6, hjust = 0)
	}
}
