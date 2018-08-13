

#' Get a biostrig object from different souces
#'
#' @param seqs loction of a fasta file or a character object in memory
#' @param ... additonal arguments to be passed to other functions
#'
#' @return an XStringSet object
#' @export
#'
#' @examples
get.seqs <- function(seqs, ...) {
    UseMethod("get.seqs")
}


#' @describeIn get.seqs
get.seqs.XStringSet <- function(seqs, ...) {
    return(seqs)
}

#' @describeIn get.seqs
get.seqs.character <- function(seqs, ...) {
    if (all(file.exists(seqs))) return(Biostrings::readBStringSet(seqs, ...))
    return(Biostrings::BStringSet(seqs, ...))
}


#' Cluster sequences by calling Clustal Omega
#'
#' @param seqs either
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
cluster.sequences <- function(seqs, clustalo.exe = 'clustalo',
                              out.distfile = tempfile(), out.fastafle = tempfile(),
                              in.fastafile = tempfile(), ...) {
    UseMethod("cluster.sequences")
}

#' @describeIn cluster.sequences
cluster.sequences.XStringSet <- function(seqs, clustalo.exe = 'clustalo',
                                         out.distfile = tempfile(), out.fastafle = tempfile(),
                                         in.fastafile = tempfile(), ...) {
    print(glue::glue("Generated temp fasta in {in.fastafile}"))
    Biostrings::writeXStringSet(seqs, in.fastafile)
    cluster.sequences.character(in.fastafile, ...)
}

#' @describeIn cluster.sequences
cluster.sequences.character <- function(seqs, clustalo.exe = 'clustalo',
                                      out.distfile = tempfile(), out.fastafle = tempfile(),
                                      ...) {

    if (!all(file.exists(seqs))) return(cluster.sequences.XStringSet(Biostrings::BStringSet(seqs, ...)))
    return(cluster.sequences.default(seqs))

}


#' @describeIn cluster.sequences
cluster.sequences.default <- function(seqs, clustalo.exe = 'clustalo',
                                        out.distfile = tempfile(), out.fastafle = tempfile(),
                                        ...) {

    stopifnot(file.exists(seqs))

    if (127 == suppressWarnings({system2(clustalo.exe, args = "--help", stdout = NULL)})) {
        stop("Clustal Omega executable no found")
    }

    system2(clustalo.exe,
            glue::glue("--full -i {seqs}",
                       " --distmat-out {out.distfile}",
                       " -o {out.fastafle} --force" ))

    print(glue::glue("Clustal run finished, \n",
                     "Distance matrix stored in ",
                     "{out.distfile}\n",
                     "Output fasta locatd in ",
                     "{out.fastafle}"))

    return(
        list(out.fastafle = out.fastafle,
             out.distfile = out.distfile)
    )
}

# Function to parse the CO output
#
parse.clustal.dist.matrix <- function(file) {
    num_values <- as.numeric(scan(file, n = 1))
    distmat <- data.table::fread(file, skip = 1)
    if (is.character(distmat[[1]])) {
        rownames(distmat) <- distmat[[1]]
        colnames(distmat) <- rownames(distmat)
    }
    distmat <- as.dist(data.matrix(distmat))
    return(distmat)
}


#' Cluster and draw sequences
#'
#' Uses clustal omega to custer the sequeces prvided in a fata file and
#' Makes sequence logos of the cropped branches
#'
#' @param fastafile file, location of the fasta file
#' @param cutoff double, height to use as a cutoff for the hierarchial clustering
#' @param distances dist, input ditance matrix, if not providd the fasta file will be clustered and will be given internally
#' @param clustalo.exe sting, location of the custal omega executable or name of the command that calls it in your path
#'
#' @return ggplot object
#'
#' @examples
#'
#'
#' @importFrom data.table fread
#' @importFrom Biostrings readAAStringSet writeXStringSet AAStringSet
#' @importFrom ggdendro dendro_data ggdendrogram
#' @importFrom ggseqlogo ggseqlogo theme_logo
#' @importFrom glue glue
#' @importFrom ggplot2 theme ggplot aes element_blank
#' @importFrom magrittr "%>%"
#' @import patchwork
#' @export
draw.msa.dendrogam <- function(
    sequences,
    input.type = c("Fasta", "XStingSet"),
    distances = NULL,
    cutoff = 0.25,
    clustalo.exe = 'clustalo') {

    if (is.null(distances)) {
        out <- cluster.sequences(sequences, clustalo.exe = clustalo.exe)

        sequences <- Biostrings::readAAStringSet(out$out.fastafle)
        distances <- parse.clustal.dist.matrix(out$out.distfile)
    }

    sequences <- get.seqs(sequences)

    attr(distances, "Labels") <- as.character(sequences)

    mydendro <- as.dendrogram(hclust(distances))
    cutdendro <- cut(mydendro, h = cutoff)
    dendrodata <- ggdendro::dendro_data(cutdendro$upper)
    mygdendro.gg <- ggdendro::ggdendrogram(dendrodata, rotate = TRUE) +
        theme(axis.title.x = element_blank(),
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank())

    # TODO add a way to lbel and ordr the web logos
    if (length(cutdendro$lower) > 20) {warning("More than 20 tree branches, review if you want that many")}
    myplots <- ggseqlogo(
        lapply(cutdendro$lower, labels),
        ncol = 1) | mygdendro.gg

    return(myplots)
}




