#' @import Biostrings
#' @import data.table
#' @import GenomeInfoDb
#' @import GenomicRanges
#' @import ggplot2
#' @import gTrack
#' @import gGnome
#' @import gUtils
#' @import IRanges
#' @import magrittr
#' @importFrom VariantAnnotation readVcf info
#' @importFrom MatrixGenerics rowRanges
#' @importFrom BSgenome.Hsapiens.UCSC.hg19 Hsapiens
#' @importFrom igraph graph.adjacency shortest_paths
#' @importFrom khtools .gc setcols gr2df df2gr
#' @importFrom rtracklayer import import.bw
#' @importFrom skitools rel2abs ppng

#' @name easy.cut
#' @title easy.cut
#'
#' @description
#' wrapper around cut to save literally one line of code
#'
#' @param x (numeric) vector of values
#' @param start (numeric) lowest number to start at, default -1e3
#' @param end (numeric) highest number to end at, default - start
#' @param step (numeric) step size, default 100
#'
#' @return numeric vector with bin midpoints
#' @export
easy.cut = function(x, start = 1e3, end = -start, step = 1e2) {
    y = cut(x, breaks = seq(start, end, by = step), labels = seq(start + step / 2, end - step / 2, by = step))
    return(as.numeric(as.character(y)))
}

#' @name grab.allelic.cn.gtrack
#' @title grab.allelic.cn.gtrack
#'
#' @description
#' grab jittered gTrack showing allelic CN
#'
#' @param gg (character) path to gGraph or a GRanges
#' @param melted (logical) default TRUE
#' @param afield (character) allele field (if melted is TRUE)
#' @param cfield (cn field) CN field (if melted is FALSE)
#' @param colormap (character) named character vector
#' @param jitter (numeric) jitter major nodes up this amount and minor nodes down this amount
#' @param ywid (numeric)
#' @param ... other inputs to gTrack
#'
#' @return gTrack with desired properties
#' @export
grab.allelic.cn.gtrack = function(gg, melted = TRUE,
                                  afield = "allele",
                                  cfield = "cn", ##c("acn", "bcn"),
                                  colormap = c(major = alpha('red', 0.5), minor = alpha('blue', 0.5)),
                                  jitter = 0.18,
                                  ywid = 0.3,
                                  ...)
{

    if (inherits(gg, "GRanges")) {
        if (melted) {
            gr = gg
            values(gr)[, "allele"] = values(gr)[, afield]
            values(gr)[, "cn"] = values(gr)[, cfield]
        } else {
            dc.dt = as.data.table(gg)
            dt = rbind(dc.dt[, .(seqnames, start, end, allele = "major", cn = get(cfield[1]))],
                       dc.dt[, .(seqnames, start, end, allele = "minor", cn = get(cfield[2]))])
            gr = dt2gr(dt)
        }
    } else {
        jab = gG(jabba = gg)
        if (melted) {
            dt = jab$nodes$dt[, .(seqnames, start, end, allele = get(afield), cn = get(cfield))]
        } else {
            dt = rbind(jab$nodes$dt[, .(seqnames, start, end, allele = "major", cn = get(cfield[1]))],
                       jab$nodes$dt[, .(seqnames, start, end, allele = "minor", cn = get(cfield[2]))])
        }
        gr = dt2gr(dt)
    }
    values(gr)[, "cn.jitter"] = ifelse(values(gr)[, "allele"] == "major",
                                       values(gr)[, "cn"] + jitter,
                                       values(gr)[, "cn"] - jitter)
    gt = gTrack(gr, y.field = "cn.jitter", ywid = ywid, gr.colorfield = "allele", colormap = colormap, ...)
    return(gt)
}
                                  
                                  

#' @name starcode
#' @title starcode
#'
#' @description
#' convert p values to stars
#' 
#' @param x (numeric) vector of p values
#' @return character vector of stars corresponding with p value magnitude
starcode = function(x) {
    return( ifelse(x < 1e-4, "****", ifelse(x < 1e-3, "***", ifelse(x < 1e-2, "**", ifelse(x < 5e-2, "*", ifelse(x < 1e-1, ".", ""))))))
}

#' @name hets2roh
#' @title hets2roh
#'
#' @description
#'
#' from sites.txt get a GRanges with heterozygous sites
#'
#' @param hets (character) path to sites.txt
#' @param min.frac (numeric) minimum normal allele frequency default 0.2
#' @param exclude.centromeres (logical) default FALSE
#' @param return.type (character) GRanges or data.table
#'
#' @return GRanges of all ROH
hets2roh = function(sites,
                    min.frac = 0.2,
                    exclude.centromeres = FALSE,
                    return.type = "GRanges")
{
    if (!check_file(sites)) {
        stop("supplied file invalid")
    }

    hets.gr = gr.nochr(gr.stripstrand(grab.hets(sites, min.frac = min.frac)))

    ## just take the major alleles, these should be duplicated
    hets.gr = hets.gr %Q% (allele == "major")

    gaps.gr = gaps(hets.gr)
    sel = which(seqnames(gaps.gr) %in% c(as.character(1:22), "X", "Y"))
    gaps.gr = gaps.gr[sel] %Q% (strand == "*")

    ## get rid of centromeres
    if (exclude.centromeres) {
        centromeres.gr = grab_centromeres()
        gaps.gr = gaps.gr %Q% (!gaps.gr %^% centromeres.gr)
    }

    if (return.type == "data.table") {
        return(as.data.table(gaps.gr))
    }
    return(gaps.gr)
}

#' @name all_breakends
#' @title all_breakends
#'
#' @description
#'
#' Get all breakends from JaBbA and say whether they are loose or junctions
#'
#' @param jabba_rds (character) path to JaBbA graph
#' @param verbose (logical)
#' @param return.type (character) default GRanges
#'
#' @return GRanges or data.table
#' @export
all_breakends = function(jabba_rds, return.type = "GRanges", verbose = FALSE) {

    gg = gG(jabba = jabba_rds)

    loose = gg$loose %Q% (!terminal)
    if (length(loose)) {
        loose = loose[, c("node.id", "node.cn", "cn")]
    }

    if (length(gg$junctions) > 0 && length(gg$junctions[cn > 0 & class != "REF"]) > 0) {
        junctions = stack(gg$junctions[cn > 0 & class != "REF"]$grl)
        cols = intersect(c("node.id", "edge.id", "linkedsv"), names(values(junctions)))
        junctions = junctions[, cols]
        values(junctions)[, "node.cn"] = gg$nodes$dt[values(junctions)[, "node.id"], cn]
        values(junctions)[, "cn"] = gg$edges$dt[values(junctions)[, "edge.id"], cn]
        values(junctions)[, "jstring"] = rep(grl.string(gg$junctions[cn > 0 & class != "REF"]$grl),
                                             each = 2)
        out = grbind(loose, junctions)
        values(out)[, "loose.or.junction"] = c(rep("loose", times = length(loose)),
                                               rep("junction", times = length(junctions)))
    } else {
        out = loose
        values(out)[, "loose.or.junction"] = rep("loose", times = length(loose))
    }

    if (return.type == "data.table") {
        return(as.data.table(out))
    }

    return(out)
}

#' @name telomeric_junctions
#' @title telomeric_junctions
#'
#' @description
#'
#' Given a set of junctions
#' identify whether bp1/bp2 contains a telomere trimerf
#' on the (+) strand adjacent to the junction breakend
#'
#' @param junctions (gGnome Junction object or GRangesList)
#' @param pad (numeric) number of base pairs padding, default 101
#' @param verbose (logical) default FALSE
#'
#' @return Junctions object with logical metadata fields bp1.telomere and bp2.telomere
#' @export
telomeric_junctions = function(junctions, pad = 101, verbose = FALSE)
{
    if (!inherits(junctions, "Junction")) {
        stop("Must supply Junction object to junctions")
    }

    if (!length(junctions)) {
        if (verbose) { message("Length = 0 junctions") }
        return(junctions)
    }

    grl = junctions$grl
    bp1 = gUtils::grl.pivot(grl)[[1]]
    bp2 = gUtils::grl.pivot(grl)[[2]]

    ## helper function for getting mapqs
    grab_telos = function(bps) {

        ## add padding
        bp1.pad = GenomicRanges::resize(bps, width = pad, fix = "start")

        ## unify seqlevels with reference
        bp1.fixed = trim(gr.fix(gr.chr(bp1.pad), BSgenome.Hsapiens.UCSC.hg19::Hsapiens,drop = TRUE))
        strand(bp1.fixed) = "+"

        ## get DNAStringSet
        bp1.strings = Biostrings::getSeq(BSgenome.Hsapiens.UCSC.hg19::Hsapiens, bp1.fixed)

        g.telomeric = find_telomeres(as.character(bp1.strings), gorc = "g")
        c.telomeric = find_telomeres(as.character(bp1.strings), gorc = "c")

        return(data.table(g.telomeric = g.telomeric, c.telomeric = c.telomeric))
    }

    bp1.telos = grab_telos(bp1)
    bp2.telos = grab_telos(bp2)

    junctions = copy(junctions)
    junctions$set(bp1.g = bp1.telos[, g.telomeric],
                  bp1.c = bp1.telos[, c.telomeric],
                  bp2.g = bp2.telos[, g.telomeric],
                  bp2.c = bp2.telos[, c.telomeric])
    
    return(junctions)
}   

#' @name junction_mappability
#' @title junction_mappability
#'
#' @description
#'
#' Given a set of junctions and a specified pad
#' adds fields bp1.mapq and bp2.mapq
#' referring to the MAPQ of each breakend
#'
#' @param junctions (gGnome Junction object or GRangesList)
#' @param ref (BWA object) reference genome
#' @param pad (numeric) number of base pairs padding, default 101
#' @param verbose (logical) default FALSE
#' 
#' @return Junctions object with numeric metadata fields bp1.mapq and bp2.mapq
junction_mappability = function(junctions, ref, pad = 101,
                                build = "hg19",
                                verbose = FALSE)
{
    if (!inherits(junctions, "Junction")) {
        stop("Must supply Junction object to junctions")
    }

    if (!inherits(ref, "BWA")) {
        stop("Must supply BWA object to ref")
    }

    if (!length(junctions)) {
        if (verbose) { message("Length = 0 junctions") }
        return(junctions)
    }

    grl = junctions$grl
    bp1 = gUtils::grl.pivot(grl)[[1]]
    bp2 = gUtils::grl.pivot(grl)[[2]]

    if (build == "hg19") {
        hg = BSgenome.Hsapiens.UCSC.hg19::Hsapiens
    } else {
        hg = BSgenome.Hsapiens.UCSC.hg38::Hsapiens
    }

    ## helper function for getting mapqs
    grab_mapqs = function(bps) {

        ## add padding
        bp1.pad = trim(GenomicRanges::resize(bps, width = pad, fix = "start"))

        ## unify seqlevels with reference
        bp1.fixed = trim(gr.fix(gr.chr(bp1.pad), hg, drop = TRUE))

        ## get DNAStringSet
        bp1.strings = Biostrings::getSeq(hg, bp1.fixed)

        ## get mapq by qname of primary alignment
        bp1.alns = as.data.table(ref[bp1.strings])
        if (bp1.alns[, .N]) {
            bp1.alns[, primary := bamUtils::bamflag(flag)[, "isNotPrimaryRead"] == 0]
            bp1.mapq.dt = bp1.alns[(primary), .(mapq = mapq[1]), by = qname]
            mapqs = bp1.mapq.dt[, as.numeric(as.character(mapq))][match(1:length(bps), bp1.mapq.dt[, qname])]
        } else {
            mapqs = c()
        }
        
        return(mapqs)
    }

    bp1.mapqs = grab_mapqs(bp1)
    bp2.mapqs = grab_mapqs(bp2)

    junctions = copy(junctions)
    junctions$set(bp1.mapq = bp1.mapqs)
    junctions$set(bp2.mapq = bp2.mapqs)

    return(junctions)
}

#' @name svbnd
#' @title svbnd
#'
#' @description
#'
#' Helper function I stole from Marcin's blog file for converting ranges to breakends
#'
#' @param svs (GRanges) with annotation field SVTYPE (or character vec specified in field)
#' @param field (character) default SVTYPE
#' @param return.type (character) one of jJ, grl, default jJ
svbnd = function(svs, field = "SVTYPE", return.type = "grl") {
    out = jJ(seqlengths = seqlengths(svs))  
    if (!length(svs)) {
        if (return.type == "grl") {
            return(out$grl)
        }
        return(out)
    }

    ## check if proper annotation is present
    if (!field %in% names(mcols(svs))) {
        stop("Metadata of svs missing field: ", field)
    }
    mcols(svs)[, "SVTYPE"] = mcols(svs)[, field]
    
    dels = svs %Q% (SVTYPE == 'DEL')
    if (length(dels))
    {
        gr1 = flank(dels, start = TRUE, 1)[, c()]
        strand(gr1) = '-'
        gr2 = flank(dels, start = FALSE,  1)[, c()]
        strand(gr2) = '+'
        grl = grl.pivot(GRangesList(gr1, gr2))
        values(grl) = values(dels)
        dels = jJ(grl)
        out = c(dels, out)
    }  
    dups = svs %Q% (SVTYPE == 'DUP')
    if (length(dups))
    {
        gr1 = flank(dups, start = TRUE, 1)[, c()]
        strand(gr1) = '+'
        gr2 = flank(dups, start = FALSE,  1)[, c()]
        strand(gr2) = '-'
        grl = grl.pivot(GRangesList(gr1, gr2))
        values(grl) = values(dups)
        dups = jJ(grl)##bp1 = gr1, bp2 = gr2, meta = values(dels))
        ##dups = jJ(bp1 = gr1, bp2 = gr2, meta = values(dups))
        out = c(dups, out)
    } 
    inv = svs %Q% (SVTYPE == 'INV')
    if (length(inv))
    {
        gr1 = flank(inv, start = TRUE, 1)[, c()]
        strand(gr1) = '-'
        gr2 = gr.end(inv)
        strand(gr2) = '-'
        grl = grl.pivot(GRangesList(gr1, gr2))
        values(grl) = values(inv)
        inv1 = jJ(grl)##bp1 = gr1, bp2 = gr2, meta = values(dels))
        ##inv1 = jJ(bp1 = gr1, bp2 = gr2, meta = values(inv))
        gr1 = flank(inv, start = FALSE, 1)[, c()]
        strand(gr1) = '+'
        gr2 = gr.start(inv)    
        strand(gr2) = '+'
        inv2 = jJ(bp1 = gr1, bp2 = gr2, meta = values(inv))
        out = c(out, inv1, inv2)
        grl = grl.pivot(GRangesList(gr1, gr2))
        values(grl) = values(inv)
        inv2 = jJ(grl)##bp1 = gr1, bp2 = gr2, meta = values(dels))
    }

    if (return.type == "grl") {
        return(out$grl)
    }
    return(out)
}


#' @name dgv2grl
#' @title dgv2grl
#'
#' @description
#'
#' convert table from DGV to GRangesList
#'
#' @param dgv.fn (character) path to text file
#' @param colnums (numeric) column numbers to extract
#' @param colnames (column names)
#' @param svtypes (character) default c("deletion", "duplication", "gain", "loss") (ignores case)
#'
#' @return GRangesList with values field
dgv2grl = function(dgv.fn, colnums = numeric(), colnames = character(),
                   mcols = c(),
                   field = "varianttype",
                   svtypes = c("deletion", "gain", "loss", "duplication")) {

    dgv.dt = fread(dgv.fn)

    if (length(colnums) & length(colnames)) {
        if (length(colnums) != length(colnames)) {
            stop("Number of columns should be equal to colnames")
        }
        dgv.dt = dgv.dt[, names(dgv.dt)[colnums], with = FALSE]
        setnames(dgv.dt, colnames)
    }

    setnames(dgv.dt, c("chr", "start", "end", field), c("seqnames", "start", "end", "SVTYPE"))

    dgv.dt = dgv.dt[tolower(SVTYPE) %in% svtypes,]
    dgv.gr = GRanges(seqnames = dgv.dt[, seqnames],
                     ranges = IRanges(start = dgv.dt[, start], end = dgv.dt[, end]),
                     SVTYPE = ifelse(dgv.dt[, SVTYPE] %in% c("gain", "duplication"), "DUP", "DEL"))
    if (length(mcols)) {
        values(dgv.gr) = cbind(values(dgv.gr), dgv.dt[, ..mcols])
    }

    res = svbnd(dgv.gr, field = "SVTYPE", return.type = "jJ")
    return(res)
}

#' @name sample_granges
#' @title sample_granges
#'
#' @param gr (GRanges)
#' @param k (numeric) number of samples
#' @param width (numeric) width of each sample
#' @param verbose
sample_granges = function(gr, k, width = 1, verbose = FALSE, seqlengths = NULL, mc.cores = 8) {

    ## pick a GRanges in proportion to the width of the range
    grs = sample(1:length(gr), size = k, prob = width(gr) / sum(width(gr)), replace = TRUE)

    ## optimize a bit (e.g. order by seqnames, and count)
    dt = data.table(grix = grs)[, .N, by = grix]

    starts = mclapply(1:dt[, .N],
                      function(i, tbl) {
                          ix = tbl[i, grix]
                          n = tbl[i, N]
                          sm = sample(GenomicRanges::start(gr[ix]):pmax(GenomicRanges::end(gr[ix]) - width + 1, GenomicRanges::start(gr[ix])),
                                      size = n,
                                      replace = TRUE)
                          return(sm)
                      },
                      dt,
                      mc.cores = mc.cores)

    flat.starts = purrr::flatten(starts) %>% as.numeric
    snames = rep(as.character(seqnames(gr[dt$grix])), times = dt$N)

    if (is.null(seqlengths)) {
        seqlengths = seqlengths(gr)
    }
    
    return(GRanges(seqnames = snames, ranges = IRanges(start = flat.starts, width = width), seqlengths = seqlengths))
            
    ## pick a starting point
    ##starts = sapply(grs, function(ix) { return(sample(1:pmax((width(gr)[ix] - width + 1), 1), size = 1, replace = TRUE)) })
    ##return(GRanges(seqnames = seqnames(gr)[grs], ranges = IRanges(start = starts, width = width), seqlengths = seqlengths(gr)))
}

#' @name create_mappability_fasta
#' @title create_mappability_fasta
#'
#' @param ref (DNAStringSet or fasta)
#' @param width (numeric) bp, default 101
#' @param gr (GRanges) limit output to bases starting within these ranges
#' @param output.seqnames (character) limit output to these particular seqnames
#' @param std.chrs (logical) default TRUE
#' @param output.dir (character) path to output directory
#' @param format (character) either fasta or fastq
#' @param verbose (logical) print stuff
#'
#' @return (invisibly) returns path to output FASTA/FASTQ file
create_mappability_fasta = function(ref,
                                    width = 101,
                                    gr = NULL,
                                    output.seqnames = c(),
                                    std.chrs = TRUE,
                                    tile = FALSE,
                                    output.dir = "~/testing_tmp/mappability",
                                    format = "fastq",
                                    verbose = TRUE) {

    if (!dir.exists(output.dir)) {
        dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)
    }

    ## grab eligible seqnames from reference
    sl = seqlengths(ref)
    if (std.chrs) {
        sl = sl[which(grepl('^(chr)*[0-9XY]+$', names(sl)))]
    }

    ## intersect eligible seqnames with reference
    if (length(output.seqnames)) {
        output.seqnames = intersect(output.seqnames, names(sl))
        sl = sl[names(sl) %in% output.seqnames]
        if (!length(output.seqnames)) {
            stop("None of the desired output seqnames intersect with seqnames in the reference")
        }
    } else {
        output.seqnames = names(sl)
    }

    ## check if we should be subsetting for ranges
    if ((!is.null(gr)) && inherits(gr, "GRanges")) {
        ## subset for sequences with width > input width
        gr = gr %Q% ((as.character(seqnames(gr)) %in% output.seqnames) & (GenomicRanges::width(gr) >= width))
        sl = sl[unique(as.character(seqnames(gr)))]
    } else {
        gr = si2gr(si = sl)
    }

    if (!length(gr)) {
        stop("No valid ranges supplied")
    }

    ## gr = reduce(gr) maybe don't lol

    if (tile) {
        ## write a separate FASTA file per chromosome
        lapply(seq_along(sl),
               function(ix) {
                   message(ix)
                   ## get all the GRanges associated with target chromosome
                   current.ranges = gr[which(seqnames(gr) == names(sl)[ix])]
                   if (verbose) {message("Making GRanges for chromosome: ", names(sl)[ix])}
                   ## get all the starting points for the fasta
                   grix = lapply(seq_along(current.ranges),
                                 function(jx) {
                                     jwidth = width(current.ranges)[jx]
                                     jstart = GenomicRanges::start(current.ranges)[jx]
                                     jend = GenomicRanges::end(current.ranges)[jx]
                                     grjx = GRanges(seqnames = names(sl)[ix],
                                                    ranges = IRanges(start = jstart:(jend - width + 1),
                                                                     end = (jstart + width - 1):jend),
                                                    strand = "+",
                                                    seqlengths = sl)
                                     return(grjx)
                                 })
                   grix = do.call("grbind", grix)
                   names(grix) = as.character(GenomicRanges::start(grix))
                   if (verbose) {message("Grabbing sequences...")}
                   bs = ref[grix]
                   if (verbose) { message("Writing FASTA for chromosome: ", names(sl)[ix]) }
                   fn = file.path(normalizePath(output.dir), paste0(names(sl)[ix], ".", format))
                   if (verbose) {message("Output file name : ", fn)}
                   Biostrings::writeXStringSet(bs,
                                               filepath = fn,
                                               append = FALSE, ## separate file for each chromosome??
                                               format = format)
                   rm(bs)
                   gc()
               })
    } else {
        if (verbose) {message("Grabbing sequences...")}
        grix = gr[which(as.character(seqnames(gr)) %in% output.seqnames)]
        bs = ref[grix]
        names(bs) = start(grix)
        fn = file.path(normalizePath(output.dir), paste0(output.seqnames, ".", format))
        if (verbose) {message("Output file name : ", fn)}
        Biostrings::writeXStringSet(bs,
                                    filepath = fn,
                                    append = FALSE, ## separate file for each chromosome??
                                    format = format)
    }

           
    invisible(normalizePath(output.dir))
}

#' @name count_jabba_sv
#' @title count_jabba_sv
#'
#' @description
#'
#' Given a JaBbA graph, count how many loose ends and junctions
#'
#' @param jabba_rds (character)
#' @param id (character) default ""
#' @param verbose (logical) default FALSE
#'
#' @return data.table with columns
#' - pair
#' - n.jun
#' - n.loose
count_jabba_sv = function(jabba_rds = NA_character_,
                          id = "",
                          verbose = FALSE) {

    if (!check_file(jabba_rds)) {
        stop("Supplied jabba_rds does not exist or is empty: ", jabba_rds)
    }

    gg = gGnome::gG(jabba = jabba_rds)

    n.jun = length(gg$junctions[type == "ALT"])
    n.loose = length(gg_loose_end(jabba_rds = jabba_rds))
    purity = gg$meta$purity
    ploidy = gg$meta$ploidy

    return(data.table(pair = id,
                      n.jun = n.jun,
                      n.loose = n.loose,
                      purity = purity,
                      ploidy = ploidy))
}

#' @name check_loose_recovery
#' @title check_loose_recovery
#'
#' @description
#'
#' Given a JaBbA graph and set of true positive breakpoints:
#' check how many TP loose ends are recovered by either JaBbA junctions or loose ends
#'
#' @param jabba_rds (path to JaBbA gGraph)
#' @param le_dt (path to data table)
#' @param id (character) id of sample
#' @param pad (numeric) pad for overlaps
#' @param return.type character, c("GRanges", "data.table")
#' @param verbose (logical) default FALSE
#' @return loose ends with metadata of ov.junction and ov.loose and ov.mask


check_loose_recovery = function(jabba_rds = "/dev/null",
                                le_dt = "~/projects/lambda/db/le.rds",
                                id = "",
                                pad = 1e4,
                                mask = "~/projects/gGnome/files/zc_stash/maskA_re.rds",
                                return.type = "data.table",
                                verbose = FALSE) {

    if (!check_file(jabba_rds)) {stop ("Invalid jabba_rds file: ", jabba_rds)}
    if (!check_file(le_dt)) {stop ("Invalid le_dt file: ", le_dt)}

    if (check_file(mask)) {
        mask.gr = readRDS(mask)
    } else {
        mask.gr = GRanges()
    }

    gg = gG(jabba = jabba_rds)
    if (length(gg$junctions[type == "ALT"])) {
        bps = unlist(gg$junctions[type == "ALT"]$grl)
    } else {
        bps = GRanges()
    }

    le = gg_loose_end(jabba_rds = jabba_rds, return.type = "GRanges")

    le.dt = readRDS(le_dt)
    if (!is.null(le.dt$sample)) {
        le.dt = le.dt[(sample == id),]
    }

    if (le.dt[, .N]) {
        le.gr = dt2gr(le.dt)
        le.dt[, ov.junction := le.gr %^^% (bps + pad)]
        le.dt[, ov.loose := le.gr %^^% (le + pad)]
        le.dt[, ov.mask := le.gr %^^% (mask.gr)]
        le.dt[, recovered := ov.junction | ov.loose | ov.mask]
    }
    return(le.dt)
}

#' @name anchorlift_coverage
#' @title anchorlift_coverage
#'
#' @param bps (GRanges or data.table coercible to GRanges)
#' @param cov_rds (character) path to coverage file
#' @param normal_cov_rds (character) path to normal coverage file
#' @param field (character) coverage field, default foreground
#' @param purity (numeric) estimated purity
#' @param ploidy (numeric) estimated ploidy
#' @param return.type (character) default data.table
#' @param mask (character) path to coverage mask
#' @param window (numeric)  anchorlift window, default 1e5
#' @param ignore.strand (logical) ignore strand for anchorlift
#' @param round.digits (round anchorlift start for binning)
#' @param verbose (logical) default FALSE
#' @param ... other inputs to anchorlift
#'
#' @return anchorlifted coverage relative to breakpoints
anchorlift_coverage = function(bps = NULL,
                               cov_rds = NA_character_,
                               normal_cov_rds = NA_character_,
                               field = "foreground",
                               purity = NA,
                               ploidy = NA,
                               return.type = "data.table",
                               mask = "~/projects/gGnome/files/zc_stash/maskC.rds",
                               window = 1e5,
                               ignore.strand = FALSE,
                               round.digits = -3,
                               verbose = FALSE,
                               ...) {

    if (is.null(bps)) {
        stop("bps must be GRanges or data.table")
    }

    if (inherits(bps, "data.table")) {
        bp.gr = dt2gr(bps)
    } else {
        bp.gr = copy(bps)
    }

    if (!check_file(cov_rds)) {
        stop("invalid file for cov_rds: ", cov_rds)
    }

    if (check_file(mask)) {
        if (verbose) {
            message("Reading mask")
        }
        mask.gr = readRDS(mask)
    } else {
        mask.gr = GRanges()
    }

    cov.gr = readRDS(cov_rds)

    if (!(field %in% names(values(cov.gr)))) {
        stop("coverage missing field: ", field)
    }

    cov.gr = cov.gr[, field] %Q% (!(cov.gr %^% mask.gr))

    ## add values from normal coverage
    if (check_file(normal_cov_rds)) {
        if (verbose) { message("Reading normal coverage") }
        normal.cov.gr = readRDS(normal_cov_rds)
        names(values(normal.cov.gr))[which(names(values(normal.cov.gr)) == field)] = "normal"
        ## change column name in normal coverage
        values(cov.gr)[, "normal"] = values(cov.gr %$% normal.cov.gr[, "normal"])[, "normal"]
        values(cov.gr)[, "ratio"] = values(cov.gr)[, field] / values(cov.gr)[, "normal"]
    } else {
        message("Normal coverage not provided")
        ## don't use tumor/normal ratio
        values(cov.gr)[, "ratio"] = values(cov.gr)[, field]
    }

    ## set field to be the t/n ratio
    field = "ratio"

    ## if purity and ploidy values are supplied, add cn field
    if (!(is.na(purity) | is.na(ploidy))) {
        values(cov.gr)[, "cn"] = rel2abs(cov.gr, field = field, purity = purity, ploidy = ploidy)
    } else {
        values(cov.gr)[, "cn"] = NA
    }

    al = gUtils::anchorlift(query = cov.gr, subject = bp.gr,
                            window = window,
                            ignore.strand = ignore.strand,
                            ...)

    values(al)[, "fused"] = ifelse(start(al) < 0, "unfused", "fused")
    values(al)[, "round.start"] = round(start(al), digits = round.digits)

    ## add back information from the original breakpoint
    values(al)[, "breakpoint"] = gr.string(bp.gr[values(al)[, "subject.id"]])

    if (return.type == "data.table") {
        return(as.data.table(al))
    }
    return(al)
}
                               

#' @name compare_linked_junctions
#' @title compare_linked_junctions
#'
#' @param jj (Junction) Junction object
#' @param exclude (Junction) Junction object with junctions to exclude
#' @param jabba_rds (character) path to jabba file
#' @param ascat_seg (character) path to ascat segmentation
#' @param sequenza_seg (character) path to sequenza segmentation
#' @param germline_vcf (character) path to germline junctions file
#' @param mask_gr (character) path to mask
#' @param pad (numeric) pad around breakends, default 1e4
#' @param min_size (numeric) min span of indels, default 1e4
#' @param std_chrs (logical) default TRUE
#' @param return.type (character) one of "data.table", "GRanges"
#' @param verbose (logical) default FALSE
#'
#' @return data.table of breakpoints with columns
#' - seqnames, start, end, strand (location and orientation)
#' - id (character) id of junction in original dataset
compare_linked_junctions = function(jj = jJ(),
                                    exclude = jJ(),
                                    jabba_rds = NA_character_,
                                    ascat_seg = NA_character_,
                                    sequenza_seg = NA_character_,
                                    germline_vcf = NA_character_,
                                    mask_gr = "~/projects/gGnome/files/zc_stash/maskA_re.rds",
                                    pad = 1e4,
                                    min_size = 1e4,
                                    std_chrs = TRUE,
                                    return.type = "data.table",
                                    keep.inv = FALSE,
                                    verbose = FALSE) {

    res = data.table()

    if (!length(jj)) {
        return(res)
    }

    if (min_size > 0) {
        filtered.jj = remove_small_junctions(jj, thresh = min_size, keep.inv = FALSE)
    } else {
        filtered.jj = jJ(jj$grl)
    }

    if (std_chrs) {
        filtered.jj = jJ(standard_grl(filtered.jj$grl))
    }

    if (length(exclude)) {
        filtered.jj = setdiff(filtered.jj, exclude)
    }

    ## read supplied inputs
    if (check_file(jabba_rds)) {
        if (verbose) {
            message("Found JaBbA! Grabbing JaBbA junctions and CN change points...")
        }
        this.jab = gG(jabba = jabba_rds)
        cncp.gr = grab_cncp(this.jab$nodes$gr[, "cn"])
        jabba.jj = jabba_alt(gg = this.jab, min_size = 0)
    } else {
        if (verbose) {
            message("JaBbA not supplied - skipping")
        }
        cncp.gr = GRanges(seqlengths = seqlengths(jj))
        jabba.jj = jJ()
    }

    if (check_file(germline_vcf)) {
        if (verbose) {
            message("Found germline junctions file! Reading...")
        }
        germline.jj = jJ(germline_vcf)
    } else {
        if (verbose) {
            message("germline junctions not supplied - skipping")
        }
        germline.jj = jJ()
    }

    ascat.sequenza.segs = get_consensus(ascat = ascat_seg,
                                        sequenza = sequenza_seg,
                                        breaks = "all",
                                        pad = pad,
                                        min.width = 0,
                                        max.gap = 1e5,
                                        verbose = TRUE)
    if (check_file(mask_gr)) {
        if (verbose) {
            message("Found mask! Reading...")
        }
        mask.gr = readRDS(mask_gr)
    } else {
        if (verbose) {
            message("Mask not found! Skipping...")
        }
        mask.gr = GRanges(seqlengths = seqlengths(jj))
    }

    ## start overlapping junctions!
    if (verbose) {
        message("Computing junction overlap")
    }

    if (length(germline.jj)) {
        mg.jj = merge.Junction(x = filtered.jj,
                               y = jabba.jj,
                               z = germline.jj,
                               pad = pad,
                               all = TRUE)
    } else {
        mg.jj = merge.Junction(x = filtered.jj,
                               y = jabba.jj,
                               pad = pad,
                               all = TRUE)
    }

    if (verbose) {
        message("Preparing and annotating breakpoints")
    }

    ## collate breakpoints
    mg.grl = mg.jj[seen.by.x == TRUE]$grl
    values(mg.grl)[, "junction"] = grl.string(mg.grl)
    

    ## if no germline junctions...
    if (!("seen.by.y" %in% names(values(mg.grl)))) {
        values(mg.grl)[, "seen.by.y"] = FALSE
    }
    
    if (!("seen.by.z" %in% names(values(mg.grl)))) {
        values(mg.grl)[, "seen.by.z"] = FALSE
    }

    mg.bps = grbind(unlist(grl.pivot(mg.grl)))

    ### make sure to add back any metadata corresponding to the original lrs junctions
    ## this would be any column ending with .x
    if (verbose) { message("Grabbing original junction metadata") }
    keep.cols = names(values(mg.grl))[which(grepl(".+\\.x$", names(values(mg.grl))))]

    ## store breakpoint metadata
    if (verbose) { message("Renaming columns") }
    mj.cols = c("x", "junction", "seen.by.y", "seen.by.z")
    cols = c(mj.cols, keep.cols)
    mg.dt = setnames(as.data.table(values(mg.grl)[, cols]),
                     old = mj.cols,
                     new = c("id", "junction", "jabba.junction", "germline.junction"))
    ## strip .x suffix from column names
    new.keep.cols = sapply(keep.cols, function(x) {substring(x, first = 1, last = nchar(x) - 2)})
    setnames(mg.dt, old = keep.cols, new = new.keep.cols)
    values(mg.bps) = rbind(mg.dt, mg.dt)

    ## get junction bends
    jab.jj = this.jab$junctions[type == "ALT" & cn > 0]
    if (length(jab.jj)) {
        jab.grl = jab.jj$grl
    } else {
        jab.grl = GRangesList()
    }

    ## overlap breakpoints
    values(mg.bps)[, "mask"] = mg.bps %^% (mask.gr + pad)
    values(mg.bps)[, "junction.bp"] = mg.bps %^^% (unlist(jab.grl) + pad)
    values(mg.bps)[, "jabba.cncp"] = mg.bps %^% (cncp.gr + pad)
    values(mg.bps)[, "germline.cncp"] = mg.bps %^^% (unlist(germline.jj$grl) + pad)
    values(mg.bps)[, "ascat.sequenza"] = mg.bps %^^% (ascat.sequenza.segs + pad)
    values(mg.bps)[, "short.read.somatic"] = apply(values(mg.bps)[, c("jabba.junction",
                                                                      "jabba.cncp",
                                                                      "ascat.sequenza")], 1, any)
    values(mg.bps)[, "short.read.all"] = apply(values(mg.bps)[, c("jabba.junction",
                                                                  "germline.junction",
                                                                  "jabba.cncp",
                                                                  "germline.cncp",
                                                                  "ascat.sequenza")], 1, any)
    
    values(mg.bps)[, "private.bp"] = (!values(mg.bps)[, "short.read.all"]) & (!values(mg.bps)[, "mask"])

    mg.dt = as.data.table(mg.bps)
    ## private junction refers to junctions private to linked read candidate set
    ## shared junction refers to SVs overlapping short read SVs
    mg.dt[, private.junction := any(private.bp) & (!any(short.read.somatic)) & (!all(mask)), by = id]
    mg.dt[, shared.junction := ((!any(private.bp)) | (all(short.read.somatic))) & (!all(mask)), by = id]
    values(mg.bps)[, "private.junction"] = mg.dt[, private.junction]
    values(mg.bps)[, "shared.junction"] = mg.dt[, shared.junction]

    if (return.type == "data.table") {
        return(mg.dt)
    } 
    return(mg.bps)
}

#' @name jabba_alt
#' @title jabba_alt
#'
#' @param gg (JaBbA gGraph)
#' @param min_size (numeric) min indel size, default 0
#'
#' @return Junction object with ALT junctions
jabba_alt = function(gg = NULL, min_size = 0) {
    if (is.null(gg)) {
        return(jJ())
    }

    if (length(gg$junctions[type == "ALT"])) {
        spans = gg$junctions[type == "ALT"]$span
        classes = gg$junctions[type == "ALT"]$dt$class
        keep = (spans > min_size | classes == "INV-like" | classes == "TRA-like")
        return(jJ(gg$junctions[keep]$grl))
    } else {
        return(jJ(seqlengths = seqlengths(gg)))
    }
}


#' @name remove_small_junctions
#' @title remove_small_junctions
#'
#' @param jj (Junction) Junction object
#' @param thresh (numeric) min size (default 1e4)
#' @param keep.inv (logical) keep inversions? default TRUE
#' @return Junction object
remove_small_junctions = function(jj = jJ(), thresh = 1e4, keep.inv = TRUE) {

    if (!length(jj)) {
        return(jj)
    }

    diff.orientation = sapply(strand(jj$grl), function(x) (length(unique(x)) == 2))
    small.span = sapply(width(range(jj$grl, ignore.strand = TRUE)), function(x) (length(x) == 1 && x < thresh))
    if (keep.inv) {
        sel = which(!(diff.orientation & small.span))
    } else {
        sel = which(!(small.span))
    }

    ## return(jJ(jj$grl[sel]))
    return(jj[sel])
}

#' @name proccov
#' @title proccov
proccov = function(p, field = "tumor_cov_raw", proc.field = "^count$|^reads$", save = TRUE, template, mock = FALSE, force = FALSE, pairs) {
    if (missing(pairs)) pairs = dg(pairs)
    thisenv =environment()
    message("processing: ", thisenv$p)
    covpath = pairs[thisenv$p][[field]]
    message("covpath: ", covpath)
    bp200path = paste0(dirname(covpath), "/", p, "_200_raw_mnorm.rds")
    bp1000path = paste0(dirname(covpath), "/", p, "_1000_raw_mnorm.rds")
    browser()
    if (force || (!mock && (!check_file(bp200path) || !check_file(bp1000path)))) {
        cov = readRDS(covpath) %>% gr.nochr
        seqlevels(cov, pruning.mode = "coarse") = seqlevels(template)
        cov = sort(cov)
        cov = khtools:::.gc(cov, proc.field);
        mcols(cov) = khtools:::setcols(mcols(cov), proc.field, "reads")
        covv.new.out = khtools:::gr2df(cov)[seqnames %in% c(as.character(1:22), "X")][, -1]
        covv.new.out[, nend := pmin((ceiling((end) / 1000) * 1000), max(end)), by = seqnames]
        covv.new.out[, nstart := pmin((floor((start) / 1000) * 1000) + 1, max(start)), by = seqnames]

        ## Step2: Mean normalization
        covv.new.out[, mean.r := mean(reads, na.rm = T)]
        covv.new.out[, reads.corrected := reads/mean.r]

        if (save) saveRDS(khtools::df2gr(covv.new.out), bp200path, compress = FALSE)

        bp_1000 = covv.new.out[, .(reads = sum(reads)), by = .(seqnames, start = nstart, end = nend)]
        bp_1000[, mean.r := mean(reads, na.rm = T)]
        bp_1000[, reads.corrected := reads/mean.r]
        if (save) saveRDS(khtools::df2gr(bp_1000), bp1000path)
    }

    message("finished: ", p)
    return(data.table(pair = p,
                      tumor_cov_mnorm_200 = bp200path,
                      tumor_cov_mnorm_1k = bp1000path))

}

#' @name mask_gg_cn
#' @title mask_gg_cn
#'
#' @description
#' returns only the segments with less than mask_thresh overlapping the mask
#' 
#' @param jabba_rds (character) path to jabba rds
#' @param mask (character) path to mask
#' @param mask_thresh (numeric) between zero and one, fraction of width masked
#' @param base_thresh (numeric) min number of CN-supporting bases
#' @param return.type (character)
#' @param verbose
mask_gg_cn = function(jabba_rds,
                      mask = "~/projects/gGnome/files/lowmap.rds",
                      frac_thresh = 0.5,
                      base_thresh = 1e6,
                      return.type = "GRanges",
                      verbose = FALSE) {

    if (!check_file(jabba_rds)) {
        stop("jabba_rds not valid")
    }
    if (!check_file(mask)) {
        stop("mask not valid")
    }

    gg = gG(jabba = jabba_rds)
    mask.gr = readRDS(mask)

    if (!check_class(mask.gr, "GRanges")) {
        stop("mask must be GRanges")
    }
    
    nodes = gr.stripstrand(gg$nodes$gr)
    nodes$ov.frac = nodes %O% mask.gr
    nodes$good.bases  = (width(nodes) - nodes %o% mask.gr)
    nodes = nodes %Q% (ov.frac < frac_thresh | good.bases > base_thresh)
    if (return.type == "data.table") {
        return(as.data.table(nodes))
    }
    return(nodes)
}

#' @name mnorm_cov
#' @title mnorm_cov
#'
#' @description
#'
#' mean normalize reads and coerce to template for dryclean
#' 
#' @param cov_rds (character) path to coverage file
#' @param template.gr (GRanges) template as GRanges
#' @param template (character) path to template
#' @param field (character) default counts
#' @param chrsub (logical) get rid of chr prefix
#' @param verbose (logical)
#' 
#' @return GRanges of cov_rds with counts coerced to template and mean normalized in $reads.corrected
mnorm_cov = function(cov_rds = "/dev/null",
                     cov.gr = GRanges(),
                     template.gr = GRanges(),
                     template = "~/projects/dryclean/MPON_raw_inputs/WCM-4_1kb.rds",
                     field = "auto",
                     clean = TRUE, ## NA out GC NAs and MAP0
                     chrsub = TRUE,
                     verbose = FALSE) {

    ## if (!check_file(template)) {
    ##     stop("Invalid template")
    ## }

    ## if (!check_file(cov_rds)) {
    ##     stop("Invalid coverage file")
    ## }

    if (check_file(cov_rds)) {
        if (verbose) {
            message("Reading coverage from: ", cov_rds)
        }
        cov.gr = readRDS(cov_rds)
    }

    if (!inherits(cov.gr, "GRanges")) {
        stop("cov must be GRanges")
    }

    if (!length(cov.gr)) {
        stop("cov is empty!")
    }

    if (field == "auto") {
        mnames = names(values(cov.gr))
        field = mnames[mnames %like% "count"][1]
        if (verbose) {
            message("Using field: ", field)
        }
    }

    if (inherits(template.gr, "GRanges") && length(template.gr)) {
        message("Using supplied template")
    } else {
        if (check_file(template)) {
            if (verbose) {
                message("Reading template from: ", template)
            }
            template.gr = readRDS(template)[, c()] ## remove metadata
        } else {
            stop("No template!!")
        }
    }

    if (!inherits(template.gr, "GRanges")) {
        stop("template must be GRanges!")
    }

    if (!length(template.gr)) {
        stop("No GRanges supplied in template!")
    }

    if (chrsub) {
        if (verbose) {
            message("Replacing chr prefix")
        }
        cov.gr = gr.nochr(cov.gr)
    }

    if (clean) {
        gc.map = which(is.na(cov.gr$gc) | is.na(cov.gr$map) | cov.gr$map == 0)
        values(cov.gr)[gc.map, field] = NA
    }

    if (verbose) {
        message("Starting analysis")
    }

    ov.dt = gr.findoverlaps(template.gr, cov.gr, return.type = "data.table")
    ov.dt[, count := values(cov.gr)[[field]][subject.id]]
    template.dt = ov.dt[, .(count = sum(count, na.rm = TRUE)), by = query.id]
    template.gr$reads = template.dt$count[match(1:length(template.gr), template.dt$query.id)]
    mean.reads = mean(template.gr$reads, na.rm = TRUE)
    template.gr$reads.corrected = template.gr$reads / mean.reads

    return(template.gr)
}

#' @name grab_cnloh
#' @title grab_cnloh
#'
#' @param gr (GRanges)
#' @param a.field (character) field of A allele CN
#' @param b.field (character) field of B allele CN
#' @param ranges (logical) return actual ranges of (MR) cnloh, default TRUE. if FALSE returns just the start/end
#' @param return.type (character) data.table or GRanges
#' @param general (logical) CNLOH in the general sense - where both alleles change in CN and not necessarily have to go to zero
#' @param balanced.only (logical) only balanced rearrangemetns without change in total CN
#'
#' @return LOH starting points plus metadata column "balanced" indicated whether it is copy-neutral
grab_cnloh = function(gr,
                      a.field = "A",
                      b.field = "B",
                      ranges = TRUE,
                      return.type = "GRanges",
                      general = TRUE,
                      balanced.only = FALSE)
{

    if (!inherits(gr, 'GRanges')) {
        stop("gr must be GRanges")
    }

    if (!(a.field %in% names(values(gr))) | !(b.field %in% names(values(gr)))) {
        stop("a.field and b.field must be metadata columns of gr")
    }

    gr = copy(gr)

    gr = gr %Q% (width(gr) > 1)
    gr = gr %Q% ((!is.na(values(gr)[, a.field])) & (!is.na(values(gr)[, b.field])))

    ## sort GRanges
    seqlevels(gr) = sort(seqlevels(gr))
    gr = sort(gr)
    
    gr$A = values(gr)[[a.field]]
    gr$B = values(gr)[[b.field]]
    gr$total = gr$A + gr$B
    gr.dt = as.data.table(gr)

    ## get the previous seqname and previous CN
    gr.dt[, previous.seqnames := data.table::shift(seqnames, type = "lag")]
    gr.dt[, previous.A := data.table::shift(A, type = "lag")]
    gr.dt[, previous.B := data.table::shift(B, type = "lag")]
    gr.dt[, previous.total := data.table::shift(total, type = "lag")]

    ## check for LOH
    ## slightly more general definition
    if (general) {
        gr.dt[, start.loh := (seqnames == previous.seqnames) &
                    (((previous.A > A) & (previous.B < B)) | ((previous.A < A) & (previous.B > B)))]
        gr.dt[, end.loh := FALSE]
    } else {
        ## this is specific for cnloh
        ## e.g. the copy number of the major allele has to be 2
        gr.dt[, start.loh := (seqnames == previous.seqnames) &
                    (((previous.A > 0) & (A == 0)) | ((previous.B > 0) & (B == 0)))]
        gr.dt[, end.loh := (seqnames == previous.seqnames) &
                    (((previous.A == 0) & (A > 0)) | ((previous.B == 0) & (B > 0)))]
    }

    ## annotate each range with whether it's UPD or CNLOH
    if (general) {
        ## in general (general cnLOH) we don't care whether the total CN is 2 or not
        ## it only matters that there is now LOH
        gr.dt[, cnloh := (start.loh | end.loh) & (total == previous.total) & (A == 0 | B == 0)]
        ## in general LOH we also don't care if the total CN is 2
        gr.dt[, upd := (A == 0 | B == 0)]
    } else {
        gr.dt[, cnloh := (start.loh | end.loh) & (total == previous.total) & (total == 2) & (A == 0 | B == 0)]
        gr.dt[, upd := (A == 0 | B == 0) & (total == 2)]
    }
    ## if (ranges) {
        
    ## } else {
    ##     gr.dt[, cnloh := (start.loh | end.loh) & (total == previous.total)]
    ## }

    if (ranges) {
        loh.dt = gr.dt[, .(seqnames, start, end, strand = "*", balanced = cnloh,
                                  cnloh, upd,
                                  annotation = ifelse(cnloh, "cnloh", ifelse(upd, "upd", "other")))]
    } else {
        loh.dt = rbind(gr.dt[(start.loh), .(seqnames, start, end = start,
                                            strand = "-", balanced = cnloh)],
                       gr.dt[(end.loh), .(seqnames, start, end = start,
                                          strand = "+", balanced = cnloh)])
    }

    ## balanced.only will only give total cn-balanced instances
    ## e.g. "MR-CNLOH"
    if (balanced.only) {
        loh.dt = loh.dt[(balanced),]
    }

    if (return.type == "data.table") {
        return(loh.dt)
    }

    if (loh.dt[, .N]) {
        return(dt2gr(loh.dt, seqinfo = seqinfo(gr), seqlengths = seqlengths(gr)))
    }
    return(GRanges(seqinfo = seqinfo(gr), seqlengths = seqlengths(gr)))
}

#' @name standard_grl
#' @title standard_grl
#'
#' @param grl
#'
#' @return grl but filtered to only include standard chromosomes
standard_grl = function(grl = GRangesList()) {

    if (!length(grl)) {
        return(grl)
    }

    keep = sapply(seqnames(grl), function(x) {all(grepl('^(chr)*[0-9XY]+$', as.character(x)))})
    return(grl[keep])
}

#' @name gg2jab
#' @title gg2jab
#'
#' @description
#'
#' Hacky way to make a JaBbA-like output from a gGraph
#'
#' @param gg (gGraph)
#' @param purity (numeric) overwrite automated purity (default 1)
#' @param ploidy (numeric) overwrite automated purity ploidy calculation
#' @return list with names:
#' - segstats (GRanges, signed)
#' - adj (adjacency matrix)
#' - ab.edges (array)
#' - purity (numeric)
#' - ploidy (numeric)
#' - junctions (GRangesList)
gg2jab = function(gg, purity = NA, ploidy = NA) {

    if (!inherits(gg, 'gGraph')) {
        stop("Must supply gGraph")
    }

    ## create ab.edges object
    ab.edges = array(NA, dim = c(length(gg$junctions[type == "ALT"]), 3, 2),
                     dimnames = list(NULL, c('from', 'to', 'edge.ix'), c('+', '-')))
    ab.edges[, 1, 1] = gg$sedgesdt[sedge.id > 0][match(gg$junctions[type == "ALT"]$dt$edge.id, edge.id), from]
    ab.edges[, 2, 1] = gg$sedgesdt[sedge.id > 0][match(gg$junctions[type == "ALT"]$dt$edge.id, edge.id), to]
    ab.edges[, 3, 1] = gg$junctions[type == "ALT"]$dt$edge.id
    ab.edges[, 1, 2] = gg$sedgesdt[sedge.id < 0][match(gg$junctions[type == "ALT"]$dt$edge.id, edge.id), from]
    ab.edges[, 2, 2] = gg$sedgesdt[sedge.id < 0][match(gg$junctions[type == "ALT"]$dt$edge.id, edge.id), to]
    ab.edges[, 3, 2] = gg$junctions[type == "ALT"]$dt$edge.id

    ## prepare djacency matrix - a sparse matrix where the nonzero entries
    ## are the CNs if that field exists and 1 otherwise
    adj.dims = dim(gg$adj)
    if (!is.null(gg$sedgesdt$cn)) {
        adj = Matrix::sparseMatrix(i = gg$sedgesdt[, from],
                                   j = gg$sedgesdt[, to],
                                   x = gg$sedgesdt[, cn],
                                   dims = adj.dims)
    } else {
        adj = Matrix::sparseMatrix(i = gg$sedgesdt[, from],
                                   j = gg$sedgesdt[, to],
                                   x = 1,
                                   dims = adj.dims)
    }

    ## initalize output which is a list
    res = list(segstats = gg$gr,
               adj = adj,
               junctions = gg$junctions[type == "ALT"]$grl,
               ab.edges = ab.edges)

    ## calculate purity/ploidy if not supplied
    if (!is.na(ploidy)) {
        res$ploidy = ploidy
    } else {
        res$ploidy = weighted.mean(gg$nodes$dt[, cn], gg$nodes$dt[, width], na.rm = TRUE)
    }

    if (!is.na(purity)) {
        res$purity = purity
    } else {
        res$purity = 1
    }

    return(res)
}

#' @name egraph
#' @title egraph
#'
#' @description
#' make edge graph
#'
#' @param gg (gGraph)
#' @param altedges
#' @param thresh (numeric) distance threshold
#' @param verbose (logical) default FALSE
#' @param mc.cores (numeric) default 1
#'
#' @return adj (adjacency matrix of edge graph)
egraph = function(gg, altedges, thresh = 1e4, verbose = FALSE, chunksize = 1e30, mc.cores = 1) {

    if (verbose){
        message(sprintf('Computing junction graph across %s ALT edges with distance threshold %s',
                        length(altedges), thresh))
    }

    bp = grl.unlist(altedges$grl)[, c("grl.ix", "grl.iix")]
    bp.dt = gr2dt(bp)

    ix = split(1:length(bp), ceiling(runif(length(bp))*ceiling(length(bp)/chunksize)))
    ixu = unlist(ix)
    eps = 1e-9
    ij = do.call(rbind, split(1:length(bp), bp$grl.ix))
    xt.adj = old.adj = Matrix::sparseMatrix(1, 1, x = 0, dims = rep(length(bp), 2))

    if (!exists(".INF")) {
        .INF = pmax(sum(seqlengths(gg)), 1e9)
    }
    xt.adj[ixu, ] = do.call(rbind,
                            mclapply(ix,
                                     function(iix)
                                     {
                                         if (verbose>1)
                                             cat('.')
                                         tmpm =
                                             gr.dist(bp[iix],
                                                     gr.flipstrand(bp),
                                                     ignore.strand = FALSE)+eps
                                         return(as(tmpm, "Matrix"))
                                     },
                                     mc.cores = mc.cores))
    ## get back to marcin's version
    adj = xt.adj
    adj[which(is.na(as.matrix(adj)))] = 0
    adj[which(as.matrix(adj)>thresh)] = 0

    xt.adj[which(is.na(as.matrix(xt.adj)))] = .INF + 1
    ## two breakpoints of the same junction should be distance 1
    bp.pair = t(
        sapply(unique(bp$grl.ix),
               function(ix){
                   matrix(which(bp$grl.ix==ix), ncol=2, nrow=1)
               }))
    xt.adj[bp.pair] = 1

    if (verbose) {
        message("starting clustering")
    }

    ## do single linkage hierarchical clustering within `range`
    edist = as.dist(xt.adj)
    edist[which(is.na(edist))] = 0
    hcl = stats::hclust(edist, method = "single")
    hcl.lbl = cutree(hcl, h = thresh)
    bp.dt$hcl = hcl.lbl
    bp.hcl =
        bp.dt[,.(hcl.1 = .SD[grl.iix==1, hcl],
                 hcl.2 = .SD[grl.iix==2, hcl]),
              keyby=grl.ix]

    ## sometimes two breakpoints belong to diff hcl
    ## merge them!
    altedges$mark(hcl.1 = bp.hcl[.(seq_along(altedges)), hcl.1])
    altedges$mark(hcl.2 = bp.hcl[.(seq_along(altedges)), hcl.2])
    hcl.ig = igraph::graph_from_edgelist(
        bp.hcl[, unique(cbind(hcl.1, hcl.2))], directed = FALSE)
    hcl.comp = components(hcl.ig)
    altedges$mark(ehcl = as.integer(hcl.comp$membership)[bp.hcl[, hcl.1]])

    ## connect to MI's code
    adj[adj>thresh] = 0

    ## check bp pairs to see if they are actually reference connected (ignore.strand = TRUE)
    ## on the given graphs ...
    ## which if we have many graphs overlapping the same intervals
    ## may not actually be the case
    ## we only check connectivity using ref edges

    ## compute reference graph distance and
    ## remove any bp pairs that are farther away
    ## on the reference graph than on the
    ## linear reference

    ##FIX ME: can't handle when there are no reference edges
    refg = gg[, type == "REF"]##self[, type == 'REF']
    bpp = Matrix::which(adj!=0, arr.ind = TRUE)

    bpp1 = unique(bpp[,1])
    bpp2 = unique(bpp[,2])
    ## dref = gGnome:::pdist(bp[bpp[,1]], bp[bpp[,2]])
    dref.unique = gGnome:::pdist(bp[bpp1], bp[bpp2])
    ##drefg = diag(refg$dist(bp[bpp[,1]], bp[bpp[,2]]))
    drefg.unique = diag(refg$dist(bp[bpp1], bp[bpp2]))
    dref = dref.unique[match(bpp[, 1],bpp1)]
    drefg = drefg.unique[match(bpp[, 1],bpp1)]

    ix = drefg>dref
    if (any(ix))
        adj[bpp[ix,, drop = FALSE]] = FALSE
    if (verbose>1)
        cat('\n')

    adj = adj | Matrix::t(adj) ## symmetrize


    ## bidirected graph --> skew symmetric directed graph conversion
    ## split each junction (bp pair) into two nodes, one + and -
    ## arbitrarily call each bp1-->bp2 junction is "+" orientation
    ## then all odd proximities adjacent to bp1 will enter the "+"
    ## version of that junction and exit the "-" version

    ## new matrix will be same dimension as adj
    ## however the nodes will represents + and -
    ## orientation of junctions
    ## using the foollowing conversion

    ## i.e.
    ## bp2 --> bp1 + +
    ## bp2 --> bp2 + -
    ## bp1 --> bp1 - +
    ## bp1 --> bp2 - -

    ## we'll use the same indices just to keep things confusing
    junpos = bp1 = bp$grl.iix == 1
    junneg = bp2 = bp$grl.iix == 2

    adj2 = adj & FALSE ## clear out adj for new skew symmetric version
    adj2[junpos, junpos] = adj[bp2, bp1]
    adj2[junpos, junneg] = adj[bp2, bp2]
    adj2[junneg, junpos] = adj[bp1, bp1]
    adj2[junneg, junneg] = adj[bp1, bp2]

    return(adj2)
} 



#' @name ecycles
#' @title ecycles
#'
#' @description
#'
#' Return a data table of all ALT edges in a gGraph that is part of a simple cycle
#'
#' @param junctions (character) path to junctions file or GRangesList with field tier
#' @param jj (Junctions) if supplied, overrides junctions path
#' @param grl (GRangesList) if supplied, overrides junctions
#' @param tfield (character)
#' @param jabba_rds (character) path to JaBbA output
#' @param small.inv (logical) allow small inversions? default TRUE
#' @param thresh (numeric) distance threshold for reciprocality, default 1e3 (bp)
#' @param min.span (numeric) thres for removing small dups and dels (1e3 bp)
#' @param min.length (numeric) min number of unique edges in a cycle/path (default 2)
#' @param chunksize (numeric) default 1e30
#' @param tier2.only (logical) only keep cycles with at least 1 tier 2 junction (applies to short reads)
#' @param path (logical) simple paths instead of cycles? default FALSE
#' @param standard.only (logical) default TRUE
#' @param mc.cores (numeric) default 1
#' @param verbose (logical) default FALSE
#'
#' @return data.table with columns cycle.id, edge.id, class, junction
#' @export
ecycles = function(junctions = NULL,
                   jj = NULL,
                   grl = NULL,
                   tfield = 'tier',
                   jabba_rds = NULL,
                   small.inv = TRUE, 
                   thresh = 1e4,
                   min.span = 1e4,
                   min.length = 2, ## min number of junctions in path or cycle
                   chunksize = 1e30,
                   tier2.only = TRUE,
                   path = FALSE,
                   standard.only = TRUE,
                   mc.cores = 1,
                   verbose = FALSE) {

    ## empty output
    res = data.table(cycle.id = numeric(),
                     edge.id = numeric(),
                     class = character(),
                     junction = character())

    ## check if gGraph should be created from junctions
    if (check_class(jj, 'Junction')) {
        gg = gG(junctions = jj)
    } else if (check_class(grl, 'GRangesList')) {
        gg = gG(junctions = jJ(grl))
    } else {
        if (check_file(junctions)) {
            if (verbose) {
                message("reading junctions supplied in: ", junctions)
            }
            if (grepl(".rds$", junctions)) {
                jj = readRDS(junctions)
                if (inherits(jj, 'GRangesList')) {
                    jj = jJ(jj)
                }
            } else {
                jj = jJ(rafile = junctions)
            }
            if (verbose) {
                message("Creating gGraph from junctions")
            }
            gg = gG(junctions = jj)
        } else if (check_file(jabba_rds)) {
            if (verbose) {
                message("Reading gGraph from: ", jabba_rds)
            }
            gg = gG(jabba = jabba_rds)
        } else {
            stop("Must supply junctions or jabba_rds")
        }
    }


    ## remove small dups and dels
    espans = gg$junctions[type == "ALT"]$span
    eclass = gg$junctions[type == "ALT"]$dt$class
    if (small.inv) {
        if (verbose) {
            message("Allowing small inversions but removing small dups and dels under ",
                    min.span, " bp")
        }
        ggjuncs = gg$junctions[type == "ALT"][espans > min.span |
                                              eclass == "INV-like" |
                                              eclass == "TRA-like"]
    } else {
        if (verbose) {
            message("Not allowing small inversions, dups, or dels under ", min.span, " bp")
        }
        ggjuncs = gg$junctions[type == "ALT"][espans > min.span | eclass == "TRA-like"]
    }
    if (verbose){
        message("Number of junctions after removing low-span junctions: ", length(ggjuncs))
    }

    if (standard.only) {
        if (verbose) {
            message("Using only standard chromosomes")
        }
        ggjuncs.grl = standard_grl(ggjuncs$grl)
        ggjuncs = jJ(ggjuncs.grl)
    }
    gg = gG(junctions = ggjuncs)
    
    ## code copied from eclusters
    altedges = gg$edges[type == "ALT"]

    if (!length(altedges)) {
        if (verbose) {
            message("No ALT edges!")
        }
        return(res)
    }
    
    if (verbose){
        message('starting egraph')
    }

    adj2 = egraph(gg, altedges, thresh = thresh, mc.cores = mc.cores, verbose = verbose, chunksize = chunksize)

    if (verbose)
        message(sprintf('Created basic junction graph using distance threshold of %s', thresh))

    ## get nonzero entries adjacency matrix as a data table
    dt = as.data.table(Matrix::which(adj2, arr.ind = TRUE))
    ## iterate over unique in the data table
    unique.fr = unique(dt$row)
    ## create graph from the corresponding adjacency matrix outside of loop
    gr = graph.adjacency(adj2)
    ## list tracking nodes of each cycle
    all.cy = list()
    ## number of cycles in the graph
    num.cy = 0
    if (verbose) {
        message("Number of rows with outgoing edges: ", length(unique.fr))
    }
    for (rw in unique.fr) {
        fr = dt[row == rw, col]
        if (path) {
            ## get simple paths of length at least 2 FROM that node TO any node in the graph
            cy = igraph::shortest_paths(gr, from = rw, to = V(gr))$vpath
            sel = which(sapply(cy, length) >= min.length)
        } else {
            ## get simple cycles containing at least one other edge
            cy = igraph::shortest_paths(gr, from = fr, to = rw)$vpath
            ## only include loops with more than one node to exclude self loops
            sel = which(sapply(cy, function(x) {length(unique(x))}) >= min.length)
        }
        for (ix in sel) {
            num.cy = num.cy + 1
            if (path) {
                ## don't sort but recast node ID's as numeric
                all.cy[[num.cy]] = as.numeric(cy[[ix]])
            } else {
                ## if enumerating cycles it is important to sort to remove cyclic paths
                all.cy[[num.cy]] = sort(as.numeric(cy[[ix]]))
            }
        }
    }

    if (!length(all.cy)) {
        if (verbose) {
            message("No cycles/paths found!")
        }
        return(res)
    } 
    unique.cy = unique(all.cy)

    if (verbose) {
        message("Mapping cycle junctions to original junction ID's")
    }
    bp = grl.unlist(altedges$grl)[, c("grl.ix", "grl.iix")]
    ## don't sort to preserve edge order of paths
    ## jcl = lapply(unique.cy, function(x) unique(sort(bp$grl.ix[x]))) %>% unique
    jcl = lapply(unique.cy, function(x) bp$grl.ix[x]) %>% unique
    dcl = dunlist(unname(jcl)) %>% setnames(new = c('listid', 'edges'))
    dcl[, class := altedges$dt[dcl$edges, class]]
    dcl[, junction := grl.string(altedges$grl[dcl$edges])]

    ## transfer all junction metadata
    ## get column names nonintersecting with whatever is in dcl already
    selected.cols = setdiff(colnames(altedges$dt), colnames(dcl))
    dcl = as.data.table(cbind(dcl, altedges$dt[dcl$edges, ..selected.cols]))

    ## if tier is already provided then keep it, otherwise add tier 2 by default
    if (tfield %in% colnames(altedges$dt)) {
        dcl[, tier := altedges$dt[[tfield]][dcl$edges]]
    } else {
        dcl[, tier := 2]
    }
    if (tier2.only) {
        dcl[, t2.count := sum(tier == 2), by = 'listid']
        dcl = dcl[t2.count > 0,]
    }

    ## add length of path
    if (nrow(dcl)) {
        dcl[, pathlength := .N, by = listid]
        if (path) {
            ## add index of edge in path
            dcl[, pathid := 1:.N, by = listid]
            ## is this an internal nodei n the path?
            dcl[, internal := !(pathid == pathlength | pathid == 1)]
        } else {
            dcl[, pathid := 1:.N, by = listid]
            dcl[, internal := TRUE] ## all cycles are technically internal?
        }
    } else {
        dcl[, pathlength := numeric()]
        dcl[, pathid := numeric()]
        dcl[, internal := logical()]
    }
    
    return(dcl)
}

#' @name junction_coverage
#' @title junction_coverage
#'
#' @description
#'
#' get coverage around the breakpoints of each junction
#'
#' @param rafile (character)
#' @param cov_rds (character)
#' @param standard.only (logical) default TRUE
#' @param field (character) default foreground
#' @param pad (character) default 5e4
#' @param mask (character) path to coverage mask
#' @param verbose (logical)
#'
#' @return data.table with columns:
#' - bp.index (junction index)
#' - bp (either 1 or 2, indicating first or second breakpoint)
#' - breakpoint (gr.string of breakpoint)
#' - junction (grl.string of junction)
#' - fused (indicating if coverage is from fused or unfused side)
#' - cn (coverage)
junction_coverage = function(rafile = NULL,
                             cov_rds = NULL,
                             standard.only = TRUE,
                             field = "foreground",
                             pad = 5e4,
                             min.span = 1e4,
                             small.inv = TRUE, ## include small inversions?
                             mask = "~/projects/gGnome/files/lowmap.rds",
                             verbose = FALSE) {
    if (!check_file(rafile)) {
        stop("Must supply valid junctions file")
    }

    if (!check_file(cov_rds)) {
        stop("must supply valid coverage")
    }

    ## read junctions
    if (verbose) {
        message("Reading junctions from RAfile: ", rafile)
    }
    
    if (grepl(".rds$", rafile)) {
        jj = readRDS(rafile)
        if (inherits(jj, 'GRangesList')) {
            jj = jJ(jj)
        }
    } else {
        jj = jJ(rafile = rafile)
    }

    ## get just standard chroms
    grl = jj$grl
    if (standard.only) {
        grl = standard_grl(grl)
        jj = jJ(grl)
    }
    gg = gG(junctions = jj)

    ## get things with correct span
    jj = gg$junctions[gg$junctions$span > min.span | gg$junctions$dt$class == "INV-like"]
    gg = gG(junctions = jj)
    grl = gg$edges[type == "ALT"]$grl

    if (length(jj)) {

        if (verbose) {
            message("Number of junctions to check: ", length(grl))
        }
        
        piv = grl.pivot(gg$edges[type == "ALT"]$grl)
        altedges = gg$edges[type == "ALT"]
        ## get coverage for breakpoint 1
        ov1 = breakpoint_coverage(piv[[1]], cov_rds, field = field, pad = pad, mask = mask, verbose = verbose)
        ov1[, ":="(bp = 1, breakpoint = gr.string(piv[[1]])[bp.index], junction = grl.string(altedges$grl)[bp.index])]
        
        ## get coverage for breakpoint 2
        ov2 = breakpoint_coverage(piv[[2]], cov_rds, field = field, pad = pad, mask = mask, verbose = verbose)
        ov2[, ":="(bp = 2, breakpoint = gr.string(piv[[2]])[bp.index], junction = grl.string(altedges$grl)[bp.index])]

        return(rbind(ov1, ov2))
    }

    if (verbose) {
        message("No junctions!")
    }
    return(data.table())
}

#' @name breakpoint_coverage
#' @title breakpoint_coverage
#'
#' @description
#'
#' Get coverage bins within a certain distance of fused and unfused sides of breakpoints
#'
#' @param bps (GRanges) breakpoints
#' @param cov_rds (character) coverage file
#' @param field (character) default foreground
#' @param pad (numeric) default 5e4
#' @param mask (character) path to coverage mask
#' @param verbose (logical) default FALSE
#'
#' @return data.table with fields:
#' - cn (from coverage)
#' - bp.index (index from breakpoints)
breakpoint_coverage = function(bps = GRanges(),
                               cov_rds = NULL,
                               field = "foreground",
                               pad = 5e4,
                               mask = "~/projects/gGnome/files/lowmap.rds",
                               verbose = FALSE) {

    ## read coverage file
    if (!check_file(cov_rds)) {
        stop("Invalid coverage file")
    }
    cov = readRDS(cov_rds)

    if (!field %in% names(values(cov))) {
        stop("coverage metadata missing field: ", field)
    }

    ## read mask and check masked bins
    if (check_file(mask)) {
        mask.gr = readRDS(mask)
    } else {
        mask.gr = GRanges(seqlengths = seqlengths(cov))
    }

    ## remove bins with no coverage or masked
    cov = cov %Q% (!(cov %^% mask.gr))
    cov = cov %Q% (!is.na(values(cov)[[field]]))

    ## get fused and unfused sides of each breakpoint
    fused = resize(bps, width = pad, fix = 'start')
    unfused = resize(bps, width = pad, fix = 'end')
    fused$fused = 'fused'
    fused$bp.index = 1:length(fused)
    unfused$fused = 'unfused'
    unfused$bp.index = 1:length(unfused)
    bp.sides = c(fused, unfused)
    ov = gr.findoverlaps(bp.sides, cov, return.type = 'data.table')
    if (nrow(ov)) {
        ov[, fused := bp.sides$fused[query.id]]
        ov[, bp.index := bp.sides$bp.index[query.id]]
        ov[, cn := values(cov)[[field]][subject.id]]
        ov = ov[!is.na(cn) & cn > 0,]
    } else {
        res = ov
    }
    return(ov)
}

#' @name grab.hets.gtrack
#' @title grab.hets.gtrack
#'
#' @description
#'
#' create gTrack for hets
#'
#' @param hets (character) variants text file
#' @param jab (character)
#' @param ... additional params for gTrack
grab.hets.gtrack = function(hets = NULL,
                            jab = NULL,
                            purity = NULL,
                            ploidy = NULL,
                            max.ranges = 5e3,
                            ...) {

    if (!check_file(hets)) {
        stop("Invalid hets file")
    }

    args = list(...)
    
    hets.gr = grab.hets(agt.fname = hets)
    hets.gr$col = ifelse(hets.gr$allele == "major", "red", "blue")

    hets.gr$cn = hets.gr$count
    if (!is.null(purity) & !is.null(ploidy)) {

        ## use supplied purity and ploidy
        hets.gr$cn = skitools::rel2abs(hets.gr, field = "count", purity = purity, ploidy = ploidy, allele = TRUE)

    } else if (check_file(jab)) {
        jab.list = gG(jabba = jab)
        if (is.null(jab.list$meta$purity) | is.null(jab.list$meta$ploidy)) {
            stop("invalid jabba output missing purity/ploidy")
        }
        hets.gr$cn = skitools::rel2abs(hets.gr, field = "count", purity = jab.list$meta$purity, ploidy = jab.list$meta$ploidy, allele = TRUE)
    } 

    hets.gt = gTrack(hets.gr, y.field = "cn", circles = TRUE, lwd.border = 0.2, max.ranges = max.ranges, ...)
    hets.gt$legend.params = list(plot = FALSE)
    return(hets.gt)
}

#' @name grab.hets
#' @title grab.hets
#'
#' @description
#'
#' returns allele gtrack given sites.txt from het pileup
#'
#' @param agt.fname (character) path to sites.txt
#' @param min.frac (numeric) between 0 and 1, min frequency in normal to count as het site
#' 
#' @return allele gTrack
grab.hets = function(agt.fname = NULL,
                     min.frac = 0.2)
{
    if (is.null(agt.fname) || !file.exists(agt.fname)) {
        stop("agt.fname does not exist")
    }

    if (grepl(".txt$", agt.fname)) {
        agt.dt = fread(agt.fname)[alt.frac.n > min.frac & ref.frac.n > min.frac,]
        agt.dt[, ":="(alt.count.t = as.numeric(alt.count.t),
                      ref.count.t = as.numeric(ref.count.t),
                      alt.count.n = as.numeric(alt.count.n),
                      ref.count.n = as.numeric(ref.count.n))]
        
        ## add major and minor
        agt.dt[, which.major := ifelse(alt.count.t > ref.count.t, "alt", "ref")]
        agt.dt[, major.count := ifelse(which.major == "alt", alt.count.t, ref.count.t)]
        agt.dt[, minor.count := ifelse(which.major == "alt", ref.count.t, alt.count.t)]
    } else if (grepl(".rds$", agt.fname)) {
        agt.dt = grab.hets.from.maf(agt.fname)
    }

    ## melt the data frame
    agt.melted = rbind(agt.dt[, .(seqnames, start, end, count = major.count, allele = "major")],
                       agt.dt[, .(seqnames, start, end, count = minor.count, allele = "minor")]
                       )

    ## make GRanges
    agt.gr = dt2gr(agt.melted[, .(seqnames, start, end, count, allele)])
    agt.gr$count = as.numeric(agt.gr$count)

    return (agt.gr)
}

#' @name grab.cov.gtrack
#' @title grab.cov.gtrack
#'
#' @description
#'
#' create coverage gTrack
#'
#' @param cov (character) path to cov.rds
#' @param jab (character) path to jabba output
#' @param purity (numeric) estimated purity (overridden by JaBbA)
#' @param ploidy (numeric) estimated ploidy (overridden by JaBbA)
#' @param lwd.border (numeric)
#' @param max.ranges (numeric)
#' @param name (character)
#' @param mask (character) path to coverage mask
#' @param field (character) field in coverage file
#' @param ... additional params for gTrack
#'
#' @return gTrack of coverage
#' @export
grab.cov.gtrack = function(cov, jab = NULL,
                           purity = NULL,
                           ploidy = NULL,
                           field = "ratio",
                           lwd.border = 0.2,
                           max.ranges = 1e4,
                           name = "cov",
                           mask = "~/projects/gGnome/files/zc_stash/maskA_re.rds",###"~/projects/gGnome/files/zc_stash/wide.pon.mask.re.rds",##"~/projects/gGnome/files/zc_stash/pon.mask.re.rds", ##"~/projects/gGnome/files/zc_stash/maskA_re.rds",
                           chrsub = TRUE,
                           hidemask = FALSE,
                           ...) {
    if (grepl("txt$", cov) | grepl("csv$", cov)) {
        cov.gr = dt2gr(fread(cov))
    } else if (grepl("rds", cov)) {
        cov.gr = readRDS(cov)
    } else if (grepl("bw$", cov)) {
        cov.gr = rtracklayer::import.bw(cov)
        names(values(cov.gr)) = field
    } else {
        stop("invalid file name")
    }

    if (chrsub) {
        cov.gr = gr.nochr(cov.gr)
    }

    if (!field %in% names(values(cov.gr))) {
        stop("invalid field supplied")
    }

    if (!is.null(jab) && file.exists(jab)) {
        purity = readRDS(jab)$purity
        ploidy = readRDS(jab)$ploidy
    }

    ## NA anything that's masked
    if (check_file(mask)) {
        mask.gr = readRDS(mask)
        ## set colors if masked
        if (hidemask) {
            cov.gr = cov.gr %Q% (!gr.nochr(cov.gr) %^% gr.nochr(mask.gr))
        } else {
            values(cov.gr)[, "col"] = ifelse(gr.nochr(cov.gr) %^% gr.nochr(mask.gr), "red", "black")
        }
    } else {
        values(cov.gr)[, "col"] = "black"
    }


    if (!is.null(purity) && !is.null(ploidy)) {
        ##nomask = which(cov.gr$col == "black")
        cov.gr$cn = values(cov.gr)[[field]]
        values(cov.gr)[, "cn"] = skitools::rel2abs(gr.nochr(cov.gr),##[nomask],
                                                   field = field,
                                                   purity = purity,
                                                   ploidy = ploidy,
                                                   allele = FALSE)
    } else {
        cov.gr$cn = values(cov.gr)[[field]]
    }


    cov.gt = gTrack(cov.gr, y.field = "cn",
                    circles = TRUE,
                    lwd.border = lwd.border,
                    max.ranges = max.ranges,
                    name = name,
                    ...)
    cov.gt$legend.params = list(plot = FALSE)
    return(cov.gt)
}

#' @name get_consensus
#' @title get_consensus
#'
#' @description
#'
#' Get consensus breakends between sequenza and ascat
#' 
#' @param ascat (character) path to ascat segnemtation
#' @param sequenza (character) path to sequenza segmentation
#' @param breaks (character) 'ascat' or 'sequenza' (which breakpoints to use)
#' @param pad (numeric)
#' @param min.width (numeric) minimum width in bp for an SV default 5e4
#' @param verbose (logical)
#' 
#' @return GRanges of consensus loose ends
get_consensus = function(ascat = NA_character_,
                         sequenza = NA_character_,
                         breaks = 'ascat',
                         pad = 1e5,
                         min.width = 5e4,
                         max.gap = 5e4,
                         verbose = FALSE) {
    if (!file.exists(ascat) | !file.info(ascat)$size) {
        warning("ascat file not valid")
        ascat.loose = GRanges()
    } else {
        ascat.segs = grab_segs(bw = ascat,
                               field = 'score',
                               simplify = TRUE,
                               verbose = verbose)
        ascat.loose = grab_loose_from_seg(ascat.segs,
                                          min.width = min.width,
                                          max.gap = max.gap,
                                          return.type = 'GRanges')
    }
    
    if (!file.exists(sequenza) | !file.info(ascat)$size) {
        warning("sequenza file not valid")
        sequenza.loose = GRanges()
    } else {
        sequenza.segs = grab_segs(csv = sequenza,
                                  field = 'CNt',
                                  simplify = TRUE,
                                  verbose = verbose)
        sequenza.loose = grab_loose_from_seg(sequenza.segs,
                                             min.width = min.width,
                                             max.gap = max.gap,
                                             return.type = 'GRanges')
    }
    
    if (!(breaks %in% c('ascat', 'sequenza', 'all'))) {
        warning('breaks must be one of ascat, sequenza. using ascat by default')
        breaks = 'ascat'
    }

    ## filter for overlapping
    if (breaks == 'ascat') {
        keep.le = ascat.loose %^^% (sequenza.loose + pad)
        le = ascat.loose[keep.le]
    } else if (breaks == "sequenza") {
        keep.le = sequenza.loose %^^% (ascat.loose + pad)
        le = sequenza.loose[keep.le]
    } else {
        if (length(ascat.loose) == 0 & length(sequenza.loose) == 0) {
            le = GRanges()
        } else {
            le = grbind(ascat.loose, sequenza.loose)
        }
    }
    if (is.null(le)) {
        le = GRanges()
    }
    return(le)
}

#' @name cbs_correlation
#' @title cbs_correlation
#'
#' @param cov (character) path to coverage file
#' @param replication (character) path to replication timing GRanges
#' @param bin.size (numeric) median rebin
#'
#' @return list with three elements
#' - data: data.table of bins with foreground, replication timing
#' - pearson: pearson correlation
#' - spearman: spearman correlation
cbs_correlation = function(cov = NA_character_,
                           replication = "~/DB/fishHook/covariates/replication.timing.rds",
                           field = "foreground",
                           rep.field = "val",
                           log = TRUE,
                           bins.size = 1e6) {

    if (!file.exists(cov) && file.info(cov)$size) {
        stop("cov must be nonempty file")
    }
    
    cov.gr = readRDS(cov)

    if (is.null(values(cov.gr)[[field]])) {
        stop("supplied field must be in cov metadata: ", field)
    }

    cov.gr = cov.gr %Q% (!is.na(values(cov.gr)[[field]])) %Q% (values(cov.gr)[[field]] > 0)

    if (!length(cov.gr)) {
        stop("No nonzero/non-NA values for field: ", field)
    }

    if (log) {
        ## log transform the foreground
        values(cov.gr)[[field]] = log(values(cov.gr)[[field]])
        ## set infinite values to NA
        values(cov.gr)[[field]][is.infinite(values(cov.gr)[[field]])] = NA
    }

    rep.gr = readRDS(replication)

    if (is.null(values(rep.gr)[[rep.field]])) {
        stop("supplied rep.field must be in replication timing metadata: ", rep.field)
    }

    tiles = gr.tile(si2gr(seqinfo(cov.gr)), width = bins.size)
    out.gr = gr.val(tiles, cov.gr, val = field,
                    FUN = function(x, w, na.rm) {return(median(x, na.rm = TRUE))},
                    na.rm = TRUE,
                    default.val = NA)
    out.gr = gr.val(out.gr, rep.gr, val = rep.field, mean = TRUE, na.rm = TRUE)

    res = list(data = as.data.table(out.gr),
               pearson = cor(values(out.gr)[[field]], values(out.gr)[[rep.field]],
                             use = "na.or.complete", method = "pearson"),
               spearman = cor(values(out.gr)[[field]], values(out.gr)[[rep.field]],
                             use = "na.or.complete", method = "spearman"))

    return(res)
}

#' @name fused_unfused
#' @title fused_unfused
#' 
#' @param loose.dt (data.table)
#' @param jabba_rds (JaBbA)
#' @param verbose (logical)
#'
#' @return data.table with columns:
#' - loose.index
#' - fused.node.id
#' - unfused.node.id
#' - fused.side
fused_unfused = function(loose.dt,
                         jabba_rds = NA_character_,
                         verbose = TRUE) {

    if (!nrow(loose.dt)) {
        return(data.table(loose.index = numeric(), fused.node.id = numeric(), unfused.node.id = numeric()))
    }

    ## check that table has loose index. If it does not exist, add one.
    if (!("loose.index" %in% colnames(loose.dt))) {
        loose.dt[, loose.index := 1:.N]
    }

    if (!file.exists(jabba_rds)) {
        stop("supplied JaBbA file does not exist")
    }
    jab = gG(jabba = jabba_rds)
    sg = jab$gr##readRDS(jabba_rds)$segstats

    ## create GRanges of Loose Ends and overlap with JaBbA nodes
    loose.gr = dt2gr(loose.dt[, .(seqnames, start, end, strand)], seqlengths = seqlengths(sg))
    node.ends = gr.end(sg, ignore.strand = FALSE)[, c()]
    end.ov = gr.findoverlaps(loose.gr, node.ends, first = TRUE, ignore.strand = TRUE)
    flank.gr = flank(loose.gr, width = 1, start = TRUE) ## include DOWNSTREAM bases
    node.starts = gr.start(sg, ignore.strand = FALSE)[, c()]
    start.ov = gr.findoverlaps(flank.gr, node.starts, first = TRUE, ignore.strand = TRUE)
    

    ## create output datatable
    res = loose.dt[, .(loose.index)]
    res[, loose.subject.id := end.ov$subject.id[match(loose.index, end.ov$query.id)]]
    res[, flank.subject.id := start.ov$subject.id[match(loose.index, start.ov$query.id)]]
    if ("snode.id" %in% names(values(sg))) {
        res[, fused.node.id := abs(sg$snode.id[loose.subject.id])]
        res[, unfused.node.id := abs(sg$snode.id[flank.subject.id])]
    } else {
        res[, fused.node.id := sg$node.id[loose.subject.id]]
        res[, unfused.node.id := sg$node.id[flank.subject.id]]
    }

    return(res)
}

#' @name gg_loose_end
#' @title gg_loose_end
#'
#' @param jabba_rds (character) path to jabba output
#' @param return.type (character) one of GRanges, data.table
#' @param std.chrs (logical) return only loose ends on standard chromosomes
#' 
#' @return loose ends of a JaBbA graph as GRanges
#' if gGraph is not balanced, will produce a warning and return empty GRanges
gg_loose_end = function(jabba_rds = NA_character_,
                        return.type = "GRanges",
                        std.chrs = TRUE) {
    
    if (!check_file(jabba_rds)) {
        stop("file does not exist: ", jabba_rds)
    }
    
    gg = gG(jabba = jabba_rds)

    if (is.null(gg$nodes$dt$loose.cn.left) | is.null(gg$nodes$dt$loose.cn.right)) {
        stop("gg missing field loose.cn.left and loose.cn.right. please check if graph if balanced.")
    }

    ll = gr2dt(gr.start(gg$nodes[!is.na(cn) & loose.cn.left>0]$gr))[, ":="(lcn = loose.cn.left,
                                                                           strand = "+")]
    lr = gr2dt(gr.end(gg$nodes[!is.na(cn) & loose.cn.right>0]$gr))[, ":="(lcn = loose.cn.right,
                                                                           strand = "-")]
    if (nrow(ll)) {
        ll[, ":="(lcn = loose.cn.left, strand = "+")]
    }

    if (nrow(lr)) {
        lr[, ":="(lcn = loose.cn.right, strand = "-")]
    }

    all.loose.dt = rbind(ll, lr, fill = TRUE)

    if (all.loose.dt[, .N]) {
        stdchrs = c(as.character(1:22), "X", "Y")
        all.loose.dt = all.loose.dt[(as.character(seqnames) %in% stdchrs) |
                                    (as.character(seqnames) %in% paste0("chr", stdchrs)),]
    }


    if (return.type == "GRanges") {
        all.loose.gr = dt2gr(all.loose.dt, seqlengths = seqlengths(gg))
        return(all.loose.gr)
    }
    return(all.loose.dt)
}

#' @name grab_loose
#' @title grab_loose
#'
#' @param jabba_rds (character) path to JaBbA
#' @param normal_cov (character) path to coverage file with gc-adjusted counts
#' @param tumor_cov (character) path to dryclean with background
#' @param tumor.field (character) background coverage field (background)
#' @param norm.field (character) normal count field (reads.corrected)
#' @param mask (character) path to coverage mask
#' @param bp.mask (character) path to coverage mask
#' @param mask_pad (numeric) pad on distance to mask (in practice just for annotation purposes)
#' @param norm_thresh (numeric) 0.6 (from filter.loose)
#' @param bg_thresh (numeric) 0.9 threshold for germline CN change
#' @param min_size (numeric) minimum size in BP
#' @param subset (logical) subset loose ends to include only valid ones before returning
#' @param purity (numeric) tumor purity
#' @param ploidy (numeric) tumor ploidy
#' @param chrsub (logical) remove chr prefix of coverage? default TRUE but should be FALSE if running hg38
#' @param verbose (logical)
grab_loose = function(jabba_rds,
                      normal_cov = "/dev/null",
                      tumor_cov = "/dev/null",
                      tumor.field = "foreground",
                      norm.field = 'foreground',
                      ## cov.pad = 5e4,
                      mask = "~/projects/gGnome/files/zc_stash/tile.mask.rds",
                      bp.mask = "/dev/null", ##"~/projects/gGnome/files/zc_stash/common.bps.rds",
                      mask_pad = 0,
                      norm_thresh = 0.6,
                      bg_thresh = 0.3,
                      min_size = 0,
                      purity = NULL,
                      ploidy = NULL,
                      chrsub = TRUE,
                      verbose = TRUE) {

    if (!check_file(jabba_rds)) {
        stop("JaBbA file path not valid")
    }

    if (!check_file(normal_cov)) {
        if (verbose) {
            message("Coverage file not supplied!")
        }
        use.coverage = FALSE
    } else {
        if (verbose) {
            message("Reading normal coverage")
        }
        cov.gr = readRDS(normal_cov)
        if (chrsub) { cov.gr = gr.nochr(cov.gr) }
        use.coverage = TRUE
    }

    if (!check_file(tumor_cov)) {
        if (verbose) {
            message("Background file not supplied!")
        }
        use.bg = FALSE
    } else {
        if (verbose) {
            message("Reading background coverage")
        }
        bg.gr = readRDS(tumor_cov)
        if (chrsub) { bg.gr = gr.nochr(bg.gr) }
        use.bg = TRUE
    }

    if (use.coverage) {
        if (!(norm.field %in% names(values(cov.gr)))) {
            stop("supplied field not in coverage metadata: ", norm.field)
        }
        sel = which(!is.infinite(values(cov.gr)[, norm.field]) & !is.na(values(cov.gr)[, norm.field]))
        cov.gr = cov.gr[sel]
    }

    if (use.bg) {
        if (!(tumor.field %in% names(values(bg.gr)))) {
            stop("supplied field not in coverage metadata: ", tumor.field)
        }
        sel = which(!is.infinite(values(bg.gr)[, tumor.field]) & !is.na(values(bg.gr)[, tumor.field]))
        bg.gr = bg.gr[sel]
    }

    ## pull loose ends from jabba_rds
    if (verbose) {
        message("Grabbing loose ends from JaBbA input")
    }
    gg = gG(jabba = jabba_rds)
    le.gr = gg$loose %Q% (terminal == FALSE) 
    le.dt = as.data.table(le.gr)

    if (!nrow(le.dt)) {
        if (verbose) {
            message("No loose ends!")
        }
        return(le.dt)
    }
    
    le.gr$loose.index = 1:length(le.gr)

    ## check if loose end overlaps mask
    if (check_file(mask)) {

        if (verbose) {
            message("Checking overlaps with mask")
        }
        mask.gr = readRDS(mask)
        overlap.mask = le.gr %^% (gr.stripstrand(mask.gr) + mask_pad)
        le.dt[, mask := overlap.mask]

        ## also remove masked coverage bins!
        if (use.coverage) {
            cov.gr = cov.gr %Q% (!(cov.gr %^% gr.stripstrand(mask.gr)))
        }

        ## remove masked background bins
        if (use.bg) {
            bg.gr = bg.gr %Q% (!(bg.gr %^% gr.stripstrand(mask.gr)))
        }
        
    } else {
        le.dt[, mask := FALSE]
    }

    if (check_file(bp.mask)) {

        if (verbose) {
            message("Checking overlaps with common breakpoints!")
        }
        
        bp.mask.gr = readRDS(bp.mask)
        le.dt[, bp.mask := le.gr %^% gr.stripstrand(bp.mask.gr)]
    } else {
        le.dt[, bp.mask := FALSE]
    }

    ## get node ids of fused and unfused loose ends
    if (verbose) {
        message("Identifying fused and unfused sides of each loose end")
    }

    ## grab gGraph
    gg = readRDS(jabba_rds)
    if (!inherits(gg, 'gGraph')) {
        gg = gGnome::gG(jabba = jabba_rds)
    }

    ## re2labs transformation
    if (use.coverage) {
        if (verbose) { message("rel2abs transforming normal coverage") }
        cov.gr$abscn = skitools::rel2abs(cov.gr, field = norm.field, purity = 1, ploidy = 2)
    }

    if (use.bg) {
        if (verbose) { message("rel2abs transforming tumor coverage") }
        ##mean.ploidy = (gg$meta$purity * gg$meta$ploidy) + (1 - gg$meta$purity) * 2
        ## if (is.null(purity)) {
        ##     purity = gg$meta$purity
        ## }
        ## if (is.null(ploidy)) {
        ##     ploidy = gg$meta$ploidy
        ## }
        
        purity = gg$meta$purity
        ploidy = gg$meta$ploidy
        message("Purity: ", purity)
        message("Ploidy: ", ploidy)
        
        bg.gr$abscn = skitools::rel2abs(bg.gr,
                                        field = tumor.field,
                                        purity = purity,##gg$meta$purity,
                                        ploidy = ploidy)##gg$meta$ploidy)
    }
    
    fused.unfused.dt = fused_unfused(le.dt, jabba_rds)

    ## get flanking normal coverage (borrowed from filter.loose)
    fused.melted.dt = melt.data.table(fused.unfused.dt, id.vars = "loose.index",
                                   measure.vars = c("fused.node.id", "unfused.node.id"),
                                   variable.name = "fused",
                                   value.name = "node.id")[!is.na(node.id),]

    ## fix this!! use the max of the node width or 10 kbp
    ## grab the strand from the loose end
    fused.melted.dt[, loose.strand := le.dt[, strand][match(loose.index, le.dt[, loose.index])]]

    segs.dt = gg$nodes$dt[fused.melted.dt$node.id, .(seqnames, start, end)]
    segs.dt[, loose.strand := fused.melted.dt[, loose.strand]]
    segs.dt[, fused := fused.melted.dt[, fused]]
    segs.dt[loose.strand == "+" & fused == "fused.node.id", ":="(end = pmax(end, start + 1e4))]
    segs.dt[loose.strand == "+" & fused == "unfused.node.id", ":="(start = pmax(pmin(start, end - 1e4), 1))]
    segs.dt[loose.strand == "-" & fused == "fused.node.id", ":="(start = pmax(pmin(start, end - 1e4), 1))]
    segs.dt[loose.strand == "-" & fused == "unfused.node.id", ":="(end = pmax(end, start + 1e4))]
    
    ## segs.gr = gg$nodes$gr[fused.melted.dt$node.id, c()]
    segs.gr = trim(dt2gr(segs.dt))
    names(segs.gr) = NULL
    segs.gr$loose.index = fused.melted.dt$loose.index
    segs.gr$fused = fused.melted.dt$fused == "fused.node.id"

    if (use.coverage) {

        ## use breakpoint coverage function?
        ## only keep coverage within 1e5 of loose end?
        ## no, don't do this, use the entire node, if possible
        ## ie use the max of the entire node, vs. 10 kbp if the node for some reason is tiny
        sub.cov.gr = cov.gr ##%Q% (cov.gr %^% (le.gr))## + 5e4))
        ## sides.cov.dt = gr.findoverlaps(cov.gr, segs.gr, return.type = 'data.table')
        sides.cov.dt = gr.findoverlaps(cov.gr, segs.gr, return.type = 'data.table')
        if (nrow(sides.cov.dt)) {

            ## use rel2abs transformation
            ## sides.cov.dt[, normal.cov := values(sub.cov.gr)[query.id, "abscn"]]
            ## use foreground (this should be always positive and can be logged)
            ## sides.cov.dt[, normal.fg := values(sub.cov.gr)[query.id, norm.field]]

            sides.cov.dt[, normal.cov := values(cov.gr)[query.id, "abscn"]]
            ## use foreground (this should be always positive and can be logged)
            sides.cov.dt[, normal.fg := values(cov.gr)[query.id, norm.field]]
            sides.cov.dt[, fused := values(segs.gr)$fused[subject.id]]
            sides.cov.dt[, loose.index := values(segs.gr)$loose.index[subject.id]]
        }

    }

    if (use.bg) {

        ## add background segmentation here, which is hopefully not too noisy
        sub.bg.gr = bg.gr ##%Q% (bg.gr %^% (le.gr + 5e4))
        ## sides.cov.dt = gr.findoverlaps(cov.gr, segs.gr, return.type = 'data.table')
        bg.cov.dt = gr.findoverlaps(bg.gr, segs.gr, return.type = 'data.table')

        if (nrow(bg.cov.dt)) {
            ## bg.cov.dt[, bg.cov := sub.bg.gr$abscn[query.id]]
            ## bg.cov.dt[, fused := segs.gr$fused[subject.id]]
            bg.cov.dt[, bg.cov := bg.gr$abscn[query.id]]
            bg.cov.dt[, fused := segs.gr$fused[subject.id]]
            bg.cov.dt[, loose.index := values(segs.gr)$loose.index[subject.id]]
        }
    }

    if (use.coverage) {
        ## compute mean coverage for fused and unfused sides

        if (sides.cov.dt[, .N]) {
            mean.cov.dt = sides.cov.dt[,
                                       .(mean.normal.cov = mean(.SD$normal.cov, na.rm = T)),
                                       by = .(loose.index, fused)]
            mean.cov.dt[, fused := ifelse(fused, 'fused', 'unfused')]
            mean.cov.dt = dcast.data.table(mean.cov.dt, loose.index ~ fused, value.var = "mean.normal.cov")

            ## also do KS-test
            sides.cov.dt[, both.sides := sum(.SD$fused == FALSE) > 3 &
                               sum(.SD$fused == TRUE) > 3, by = loose.index]
            if (sides.cov.dt[(both.sides), .N]) {
                ks.cov.dt = sides.cov.dt[(both.sides),
                                         .(w = wilcox.test(.SD$normal.fg[.SD$fused == TRUE],
                                                           .SD$normal.fg[.SD$fused == FALSE])$p.value),
                                         by = .(loose.index)]
            } else {
                ks.cov.dt = data.table(loose.index = numeric(), ks = numeric(), t = numeric())
            }
            
            if ("fused" %in% colnames(mean.cov.dt) & "unfused" %in% colnames(mean.cov.dt)) {
                mean.cov.dt[, norm.change := norm_thresh < abs(fused - unfused)]
            } else {
                mean.cov.dt[, norm.change := NA]
            }

            le.dt$norm.fused = mean.cov.dt$fused[match(le.dt$loose.index, mean.cov.dt$loose.index)]
            le.dt$norm.unfused = mean.cov.dt$unfused[match(le.dt$loose.index, mean.cov.dt$loose.index)]
            le.dt$norm.change = mean.cov.dt$norm.change[match(le.dt$loose.index, mean.cov.dt$loose.index)]

            med.x = as.data.table(cov.gr)[grepl("(chr)*X$", seqnames),
                                          median(abscn, na.rm = TRUE)]
            if (verbose) { message("Median X abs CN: ", med.x) }
            male = med.x < 1.5
            if (male) {
                if (verbose) { message("Adjusting X for male sample...") }
                le.dt[grepl("(chr)*X$", seqnames),
                      norm.change := 0.5 * norm_thresh < abs(norm.fused - norm.unfused)]
            }

            le.dt$norm.w = ks.cov.dt$w[match(le.dt$loose.index, ks.cov.dt$loose.index)]
            le.dt$norm.t = ks.cov.dt$t[match(le.dt$loose.index, ks.cov.dt$loose.index)]

            ## idk, NA handling?
            le.dt[is.na(norm.fused), norm.change := TRUE]
            le.dt[is.na(norm.unfused), norm.change := TRUE]

            ## track logp value from wilcoxon and t test
            ## track log fold change and difference between fused and unfused
            le.dt[, logw := -log10(norm.w)]
            le.dt[, norm.diff := norm.fused - norm.unfused]
        } else {
            ## if loose ends overlap zero coverage... this is problematic
            le.dt$norm.change = TRUE
            le.dt$norm.fused = NA
            le.dt$norm.unfused = NA
        }
    } else {
        le.dt$norm.change = FALSE
        le.dt$norm.fused = NA
        le.dt$norm.unfused = NA
    }

    if (use.bg) {


        ## compute mean coverage for fused and unfused sides
        if (bg.cov.dt[, .N]) {
            mean.cov.dt = bg.cov.dt[, .(mean.normal.cov = mean(bg.cov, na.rm = T)),
                                       by = .(loose.index, fused)]
            mean.cov.dt[, fused := ifelse(fused, 'fused', 'unfused')]
            mean.cov.dt = dcast.data.table(mean.cov.dt, loose.index ~ fused, value.var = "mean.normal.cov")

            ## run tests
            bg.cov.dt[, both.sides := sum(.SD$fused == FALSE) > 3 &
                               sum(.SD$fused == TRUE) > 3,
                         by = loose.index]
            
            if (bg.cov.dt[(both.sides), .N]) {
                ks.cov.dt = bg.cov.dt[(both.sides),
                                      .(w = wilcox.test(.SD$bg.cov[.SD$fused == TRUE],
                                                        .SD$bg.cov[.SD$fused == FALSE])$p.value),
                                      by = .(loose.index)]
            } else {
                ks.cov.dt = data.table(loose.index = numeric(), ks = numeric(), t = numeric())
            }
            
            ## determine whether it's higher than normal beta
            if ("fused" %in% colnames(mean.cov.dt) & "unfused" %in% colnames(mean.cov.dt)) {
                mean.cov.dt[, bg.change := bg_thresh < fused - unfused] ## signed differences
            } else {
                mean.cov.dt[, bg.change := NA]
            }

            ## mark normal change and copy normal fused and unfused coverage
            le.dt$tumor.change = le.dt$loose.index %in% mean.cov.dt[(bg.change), loose.index]

            ## add p values
            le.dt$tumor.w = ks.cov.dt$w[match(le.dt$loose.index, ks.cov.dt$loose.index)]
            le.dt[, tumor.logw := -log10(tumor.w)]

            ## if there is no normal coverage, we cannot confidently tell that there is no change in normal CN
            le.dt$tumor.fused = mean.cov.dt$fused[match(le.dt$loose.index, mean.cov.dt$loose.index)]
            le.dt$tumor.unfused = mean.cov.dt$unfused[match(le.dt$loose.index, mean.cov.dt$loose.index)]
            le.dt[is.na(tumor.fused), tumor.change := TRUE]
            le.dt[is.na(tumor.unfused), tumor.change := TRUE]
            le.dt[, tumor.diff := tumor.fused - tumor.unfused]
        } else {
            le.dt$tumor.change = TRUE
            le.dt$tumor.fused = NA
            le.dt$tumor.unfused = NA
        }
    } else {
        le.dt$tumor.change = NA
        le.dt$tumor.fused = NA
        le.dt$tumor.unfused = NA
    }

    ## ## get CN of each side of loose end
    fused.unfused.dt[, fused.cn := gg$nodes$dt$cn[fused.node.id]]
    fused.unfused.dt[, unfused.cn := gg$nodes$dt$cn[unfused.node.id]]
    fused.unfused.dt[, fused.lower := (fused.cn <= unfused.cn)]

    
    ## merge these with loose ends
    le.dt = merge.data.table(le.dt, fused.unfused.dt, by = "loose.index", all.x = TRUE)

    ## BP mask is still important - the CNP might be absent in this particular matched normal but present
    le.dt[, keep := (!norm.change)]

    ## check whether each loose end is a putative false positive...
    le.dt[, possible.fp := (!fused.lower) & (tumor.fused < fused.cn) &
                           (tumor.unfused > unfused.cn) & (tumor.fused - tumor.unfused < 0.2)]


    ## add stringified loose end
    loose.end.str = gr.string(dt2gr(le.dt[, .(seqnames, start, end, strand)]))
    le.dt[, loose.end := paste0(seqnames, ":", start, strand)]

    return(le.dt)
}
           

#' @name pp_tile
#' @title pp_tile
#'
#' @description tiles an input JaBbA graph or karyograph into bins of size tile.width
#'
#' Creates a data.table representing tiled genome with column $cn
#' This can be used as input to pp_plot to QC purity/ploidy estimates
#'
#' @param seg (character) path to input segmentation
#' @param cov (character) path to coverage file
#' @param kag (character) path to karyograph
#' @param jab (character) path to JaBbA output
#' @param purity (numeric) overrides pp in jab
#' @param ploidy (numeric) overrides pp in jab
#' @param field (character) default "foreground" (signal name in cov)
#' @param mask (character) default ~/projects/gGnome/files/new.mask.rds
#' @param tile.width (numeric) default 1e5
#'
#' @return data.table of tiled genome with mean
pp_tile = function(seg = NA_character_,
                   cov = NA_character_,
                   kag = NA_character_,
                   jab = NA_character_,
                   purity = NA,
                   ploidy = NA,
                   min.quantile = 1e-3,
                   min.val = 0,
                   max.quantile = 0.999,
                   max.val = 100,
                   field = "foreground",
                   mask = "~/projects/gGnome/files/new.mask.rds",
                   return.type = "data.table",
                   tile.width = 1e5) {

    

    if (file.exists(kag)) {
        message("Reading karyograph...")
        seg.gr = gG(jabba = kag)$nodes$gr[, "cnmle"]
        seg.gr$cn = seg.gr$cnmle
    } else {
        if (!file.exists(seg)) {
            stop("seg not valid")
        }
        if (!file.exists(cov)) {
            stop("cov not valid")
        }
        if (is.null(purity) | is.na(purity) | is.null(ploidy) | is.na(ploidy)) {
            if (is.null(jab)) {
                stop("JaBbA must be supplied if purity/ploidy are not provided")
            }
            if (!file.exists(jab)) {
                stop("JaBbA must be supplied if purity/ploidy are not provided")
            }
            message("Reading purity and ploidy from JaBbA output")
            jab.list = readRDS(jab)
            purity = jab.list$purity
            ploidy = jab.list$ploidy
        }
        message("Using purity: ", purity)
        message("Using ploidy: ", ploidy)

        seg.gr = readRDS(seg)
        cov.gr = readRDS(cov)

        if (!field %in% names(values(cov.gr))) {
            stop("cov missing field in metadata: ", field)
        }

        cov.gr$rel.cn = values(cov.gr)[[field]]

        ## mask shit out first
        if (file.exists(mask)) {
            mask.gr = readRDS(mask)
        } else {
            mask.gr = GRanges(seqlengths = seqlengths(cov.gr))
        }

        na.ix = which(cov.gr %^% mask.gr)
        cov.gr$rel.cn[na.ix] = NA

        min.rel.cn = quantile(cov.gr$rel.cn, min.quantile, na.rm = TRUE)
        max.rel.cn = quantile(cov.gr$rel.cn, max.quantile, na.rm = TRUE)

        cov.gr$rel.cn[which(cov.gr$rel.cn <= min.rel.cn)] = NA
        cov.gr$rel.cn[which(cov.gr$rel.cn >= max.rel.cn)] = NA

        cov.gr$rel.cn[which(cov.gr$rel.cn <= min.val)] = NA
        cov.gr$rel.cn[which(cov.gr$rel.cn >= max.val)] = NA

        seg.gr = gr.val(seg.gr[, c()],
                        cov.gr[, "rel.cn"] %Q% (!is.na(rel.cn)),
                        mean = TRUE, na.rm = FALSE,
                        val = "rel.cn")
        
        seg.gr$nbins = seg.gr %N% (cov.gr %Q% (!is.na(rel.cn)))
        seg.gr$rel.cn[which(seg.gr$nbins < 3)] = NA
        seg.gr$cn = rel2abs(seg.gr, field = "rel.cn", purity = purity, ploidy = ploidy)

    }
    tile.gr = gr.val(gr.tile(seg.gr, width = tile.width),
                     seg.gr,
                     mean = TRUE,
                     na.rm = FALSE,
                     val = "cn")

    tile.gr$cnmle = pmax(0, round(tile.gr$cn))

    if (return.type == "data.table") {
        return(as.data.table(tile.gr))
    }

    return(tile.gr)

}


#' @name plot_pp
#' @title plot_pp
#'
#' @description create purity ploidy plot from tile output
plot_pp = function(tile.dt,
                   min.cn = -0.5,
                   max.cn = 7.5,
                   bins = 400,
                   ytrans = TRUE,
                   vlines = TRUE,
                   xlab = "Estimated CN",
                   ylab = "Count") {

    tmp = tile.dt[cn > min.cn & cn < max.cn,]
    
    pt = ggplot(tmp, aes(x = cn)) +
    geom_histogram(fill = "cornflowerblue", bins = bins, alpha = 0.8) +
    scale_x_continuous(breaks = 0:floor(max.cn),
                       labels = 0:floor(max.cn) %>% as.character)

    if (ytrans) {
        pt = pt + scale_y_continuous(trans = "log1p", breaks = c(100, 1000, 10000, 100000))
    }

    if (vlines) {
        pt = pt + geom_vline(xintercept = 0:floor(max.cn), color = "red", linetype = "longdash")
    }

    pt = pt + labs(x = xlab, y = ylab) + theme_bw()

    return(pt)
}


#' @name grab_jabba_epgaps
#' @title grab_jabba_epgaps
#'
#' @description
#'
#' create named vector (or keyed table) of jabba epgaps
#' 
#' @param pairs (keyed data.table)
#' @param mc.cores (numeric) default 8
#' @param return.table (logical) return data.table? default FALSE
#'
#' @return either named list of data.table
grab_jabba_epgaps = function(pairs = NULL, mc.cores = 8, return.table = FALSE) {

    if (!inherits(pairs, "data.table") || is.null(key(pairs))) {
        stop("pairs is not a kayed data.table")
    }

    if (key(pairs) != "pair") {
        stop("name of key must be $pair")
    }

    if (is.null(pairs$jabba_rds)) {
        stop("jabba_rds missing from pairs")
    }

    eg = mclapply(setNames(nm = pairs[file.exists(jabba_rds), pair]),
                  function (p, tbl) {
                      return(gG(jabba = pairs[p, jabba_rds])$nodes$dt$epgap[1])
                  }, pairs, mc.cores = 32)

    if (return.table) {
        return(data.table(pair = names(eg),
                          jabba_epgap = unlist(eg),
                          key = "pair"))
    }
    return(unlist(eg))
}

#' @name junctions_analysis
#' @title junctions_analysis
#'
#' @param gg (character) path to JaBbA gGraph
#' @param gs (character) path to SV GRanges
#' @param method (character) default 'jabba', others not implemented yet
#' @param max.dist (numeric) maximum distance to consider junction overlap, default 1kb
#' @param id (character)
#' @param mask (character) path to coverage mask (junctions in these regions are not considered.)
#' @param verbose (logical) default FALSE
#'
#' @param res (list with elements)
#' - junctions.jj: Junctions object
#' - junctions.res: data.table with columns pair, precision, recall
junctions_analysis = function(gg = NULL,
                              gs = "~/projects/JaBbA_pipeline/db/hcc1954.sv.gs.rds",
                              max.dist = 1e3,
                              method = "jabba",
                              id = "sample",
                              mask = "~/projects/gGnome/files/new.mask.rds",
                              verbose = FALSE) {

    if (verbose) {
        message("reading gold standard segments")
    }

    gs.cncp = readRDS(gs) ## should be GRangesList

    if (grepl(pattern = "ascat", x = method, ignore.case = TRUE)) {
        stop("method not implemented yet")
    } else if (grepl(pattern = "jabba", x = method, ignore.case = TRUE)) {
        if (is.null(gg)) {
            stop("gg must be provided")
        }
        if (!file.exists(gg)) {
            stop("gg must be provided")
        }
        junctions.grl = gG(jabba = gg)$junctions[type == "ALT" & cn > 0]$grl
    } else if (grepl(pattern = "sequenza", x = method, ignore.case = TRUE)) {
        stop("method not implemented yet")
    } else {
        stop("method not implemented yet")
    }

    if (file.exists(mask) && file.info(mask)$size) {
        mask.gr = readRDS(mask)
    } else {
        mask.gr = GRanges(seqlengths = seqlengths(gs.cncp))
    }

    ## pivot breakpoints and remove any overlapping with a masked region
    bps = grl.pivot(gs.cncp)
    bp1.ov = bps[[1]] %^% mask.gr
    bp2.ov = bps[[2]] %^% mask.gr
    gs.cncp = gs.cncp[(!bp1.ov) & (!bp2.ov)]

    junctions.bps = grl.pivot(junctions.grl)
    junctions.bp1.ov = junctions.bps[[1]] %^% mask.gr
    junctions.bp2.ov = junctions.bps[[2]] %^% mask.gr
    junctions.grl = junctions.grl[(!junctions.bp1.ov) & (!junctions.bp2.ov)]

    ## create junctions objects
    gs.jj = jJ(gs.cncp)
    junctions.jj = jJ(junctions.grl)

    ## get the overlap
    merge.jj = merge.Junction(x = gs.jj, y = junctions.jj, pad = max.dist, all = TRUE, cartesian = FALSE)

    ## calculate some stats
    n.gold = merge.jj$dt[!is.na(x), .N]
    n.call = merge.jj$dt[!is.na(y), .N]
    n.tp = merge.jj$dt[!is.na(x) & !is.na(y), .N]
    precision = n.tp / n.call
    recall = n.tp / n.gold
    f1 = 2 * (precision * recall) / (precision + recall)
    
    ## result
    junctions.res = data.table(pair = id,
                               n.gold = n.gold,
                               n.call = n.call,
                               n.tp = n.tp,
                               precision = precision,
                               recall = recall,
                               f1 = f1)

    res = list(junctions.jj = merge.jj, junctions.res = junctions.res)
    return(res)
}

#' @name cnloh_analysis
#' @title cnloh_analysis
#'
#' @description
#'
#' Compare the CNLOH points in two graphs with allelic CN annotations
#' "CNLOH" here is defined broadly: it is any site where the CN of parental alleles change in opposite directions
#' However, the total CN stays the same
#'
#' @param gg (character) path to gGraph under consideration
#' @param gs (character) path to gold standard segmentation under consideration
#' @param gg.fields (character) character vector with length two giving fields for the two alleles in gg
#' @param gs.fields (character) character vector with length two giving fields for the two alleles in gs
#' @param max.dist (character) max allowable distance to true positive call (default 1e5, generous)
#' @param new (logical) use "new" way where we benchmark based on locus overlap (default TRUE)
#' @param id (character) sample ID
#' @param method (character) input method (default gg, one of c("gr", "csv", "df", "gg"))
#' @param exclude (character) character vector of chromosomes to exclude
#' @param tile.width (numeric) tile width for considering cnloh (default 1e4)
#' @param mask (character) path to coverage mask
#' @param confusion (logical) if TRUE returns a data.table in cnloh.res that can be used to make a confusion matrix with geom_tile
#' @param verbose (logical) default FALSE
#'
#' @return list with two elements, cnloh.res with summary statistics and cnloh.gr with ranges of segments
cnloh_analysis = function(gg, gs,
                          gg.fields = c("cn.high", "cn.low"),
                          gs.fields = c("acn", "bcn"),
                          max.dist = 1e5,
                          new = TRUE,
                          id = "sample",
                          method = "gg",
                          exclude = c("X", "Y"),
                          tile.width = 1e4,
                          mask = "~/projects/gGnome/files/zc_stash/maskA_re.rds",
                          confusion = FALSE,
                          verbose = FALSE)
{
    
    gs.jab = gGnome::gG(jabba = gs)

    if (method == "gg") {
        if (new) {
            ## ASSUMES MELTED GG
            jab = gG(jabba = gg)
            method.dt = dcast.data.table(jab$nodes$dt[, .(seqnames, start, end, strand, allele, cn)],
                                         formula = seqnames + start + end + strand ~ allele,
                                         value.var = "cn")
            method.dt[, cn.high := ifelse(major >= minor, major, minor)]
            method.dt[, cn.low := ifelse(major >= minor, minor, major)]
            method.gr = dt2gr(method.dt[, .(seqnames, start, end, strand, cn.high, cn.low)])

            if (confusion) {
                method.cnloh = grab_cnloh(method.gr,
                                          a.field = "cn.high", b.field = "cn.low",
                                          ranges = TRUE,
                                          general = FALSE,
                                          balanced.only = FALSE) %>% gr.stripstrand
            } else {
                method.cnloh = grab_cnloh(method.gr,
                                          a.field = "cn.high", b.field = "cn.low",
                                          ranges = TRUE,
                                          general = FALSE, ## hmmm
                                          balanced.only = TRUE) %>% gr.stripstrand
            }

            ## method.dt = merge.data.table(adt, bdt, by = 
            ## bdt = gg$nodes$dt[(allele == "minor"), .(seqnames, start, end, strand, allele, bcn = cn)]
        } else {
            gg.jab = gGnome::gG(jabba = gg)
            method.cnloh.jj = gg.jab$junctions[(cnloh) & cn > 0]

            if (length(method.cnloh.jj)) {
                method.cnloh = gr.stripstrand(grl.pivot(method.cnloh.jj$grl)[[1]])
            } else {
                method.cnloh = GRanges()
            }
        }
    } else {

        if (grepl(pattern = "jabba", x = method, ignore.case = TRUE)) {
                method.gr = grab_allelic_segs(gg = gg,
                                              field = c("cn.high", "cn.low"),
                                              simplify = FALSE,
                                              cnloh = FALSE,
                                              loh = FALSE,
                                              melt = FALSE,
                                              verbose = verbose)
            } else if (method == "gr") {
                method.gr = grab_allelic_segs(gr = gg,
                                              field = gg.fields,
                                              cnloh = FALSE,
                                              loh = FALSE,
                                              simplify = FALSE,
                                              melt = FALSE,
                                              verbose = verbose)
            } else if (method == "csv") {
                method.gr = grab_allelic_segs(csv = gg,
                                              field = gg.fields,
                                              cnloh = FALSE,
                                              loh = FALSE,
                                              simplify = FALSE,
                                              melt = FALSE,
                                              verbose = verbose)
            } else if (method == "df") {
                method.gr = grab_allelic_segs(df = gg,
                                              field = gg.fields,
                                              cnloh = FALSE,
                                              loh = FALSE,
                                              simplify = FALSE,
                                              melt = FALSE,
                                              verbose = verbose)
            } else {
                stop("method not implemented yet")
            }

        if (new) {
            if (confusion) {
                ## setting general = FALSE forces cnloh to be a subset of loh and forces total CN to be 2
                ## setting balanced.only = FALSE returns UPD ranges in addition to CNLOH
                method.cnloh = grab_cnloh(gr = method.gr,
                                          a.field = "cn.high", b.field = "cn.low",
                                          ranges = TRUE,
                                          general = FALSE,
                                          balanced.only = FALSE) %>% gr.stripstrand
            } else {
                method.cnloh = grab_cnloh(gr = method.gr,
                                          a.field = "cn.high", b.field = "cn.low",
                                          ranges = TRUE,
                                          general = FALSE,
                                          balanced.only = TRUE) %>% gr.stripstrand
            }
        } else {
            method.cnloh = grab_cnloh(gr = method.gr,
                                      a.field = "cn.high", b.field = "cn.low",
                                      general = TRUE,
                                      balanced.only = TRUE) %>% gr.stripstrand
        }
    }

    ## if we're comparing ranges, grab segments the same way
    if (new) {
        gs.gr = grab_allelic_segs(gg = gs,
                                  field = c("acn", "bcn"),
                                  simplify = FALSE,
                                  cnloh = FALSE,
                                  loh = FALSE,
                                  melt = FALSE,
                                  verbose = verbose)
        if (confusion) {
            gs.cnloh = grab_cnloh(gr = gs.gr, a.field = "cn.high", b.field = "cn.low",
                                  ranges = TRUE,
                                  general = FALSE,
                                  balanced.only = FALSE) %>% gr.stripstrand
        } else {
            gs.cnloh = grab_cnloh(gr = gs.gr, a.field = "cn.high", b.field = "cn.low",
                                  ranges = TRUE,
                                  general = FALSE,
                                  balanced.only = TRUE) %>% gr.stripstrand
        }
    } else {

        ## otherwise we just want the location of the junction
        ## some graphs may not have CNLOH annotation - in this case just ignore
        if ("cnloh" %in% names(gs.jab$junctions$dt)) {
            gs.cnloh.jj = gs.jab$junctions[(cnloh) & cn > 0]
        } else {
            gs.cnloh.jj = jJ()
        }

        if (length(gs.cnloh.jj)) {
            gs.cnloh = gr.stripstrand(grl.pivot(gs.cnloh.jj$grl)[[1]])
        } else {
            gs.cnloh = GRanges()
        }

    }

    ## exclude segments in no-go chromosomes
    if (length(method.cnloh)) {
        method.cnloh = method.cnloh %Q% ((!(as.character(seqnames(method.cnloh)) %in% exclude)) &
                                         GenomicRanges::start(method.cnloh) > 0)
    }
    
    if (length(gs.cnloh)) {
        gs.cnloh = gs.cnloh %Q% ((!(as.character(seqnames(gs.cnloh)) %in% exclude)) & 
                                 GenomicRanges::start(gs.cnloh) > 1)
    }

    if (new) {

        if (check_file(mask)) {
            mask.gr = readRDS(mask)
        } else {
            mask.gr = GRanges()
        }
        
        tiles = gr.tile(si2gr(si = hg_seqlengths()), width = tile.width)

        if (confusion) {

            ## browser()
            ## if we're trying to return a confusion matrix there should be a nice 'annotation' column
            ## overlap with gold standard to transfer annotation
            ov.gs.dt = gr.findoverlaps(tiles, gs.cnloh, return.type = "data.table")
            ov.gs.dt[, annotation := values(gs.cnloh)[subject.id, "annotation"]]
            gs.label.dt = ov.gs.dt[, .(annotation = ifelse(any(annotation == "cnloh", na.rm = TRUE),
                                                           "cnloh",
                                                           ifelse(any(annotation == "upd", na.rm = TRUE), "upd", "other"))),
                                   by = query.id]

            ## overlap with query to transfer annotation
            ov.gr.dt = gr.findoverlaps(tiles, method.cnloh, return.type = "data.table")
            ov.gr.dt[, annotation := values(method.cnloh)[subject.id, "annotation"]]
            gr.label.dt = ov.gr.dt[, .(annotation = ifelse(any(annotation == "cnloh", na.rm = TRUE),
                                                           "cnloh",
                                                    ifelse(any(annotation == "upd", na.rm = TRUE), "upd", "other"))),
                                   by = query.id]

            ## transfer back to tiles
            values(tiles)[, "gs"] = gs.label.dt[, annotation][match(1:length(tiles), gs.label.dt[, query.id])]
            values(tiles)[, "query"] = gr.label.dt[, annotation][match(1:length(tiles), gr.label.dt[, query.id])]

            tiles = tiles %Q% ((tiles %O% mask.gr) < 0.8)
            tiles = tiles %Q% (as.character(seqnames(tiles)) != "M" & as.character(seqnames(tiles)) != "Y")
            tiles = tiles %Q% ((!is.na(values(tiles)[, "gs"])) & (!is.na(values(tiles)[, "query"])))

            ## compute grid counts
            cn.cor = as.data.table(tiles)[, .(n.tiles = .N), by = .(gs, query)]
            ## convert this to mega base pairs using tile.width
            cn.cor[, mbp := n.tiles * (tile.width / 1e6)]

            cnloh.dt = as.data.table(tiles)
            
        } else {
            values(tiles)[, "gs"] = tiles %^% gs.cnloh
            values(tiles)[, method] = tiles %^% method.cnloh
            tiles = tiles %Q% ((tiles %O% mask.gr) < 0.8)
            tiles = tiles %Q% (as.character(seqnames(tiles)) != "M" & as.character(seqnames(tiles)) != "Y")
            tiles = tiles %Q% ((!is.na(values(tiles)[, "gs"])) & (!is.na(values(tiles)[, method])))

            cn.cor = data.table(pair = id,
                                n.cnloh = sum(values(tiles)[, method]),
                                n.cnloh.gs = sum(values(tiles)[, "gs"]),
                                tp = sum(values(tiles)[, method] & values(tiles)[, "gs"]),
                                fn = sum(values(tiles)[, "gs"] & (!values(tiles)[, method])),
                                fp = sum(values(tiles)[, method] & (!values(tiles)[, "gs"])),
                                tn = sum((!values(tiles)[, method]) & (!values(tiles)[, "gs"])))

            
            cn.cor[, precision.cnloh := ifelse(tp > 0, tp / (tp + fp), 0)]
            cn.cor[, recall.cnloh := ifelse(tp > 0, tp / (tp + fn), 0)]

            ## set to NA if there are simply not any CNLOH
            cn.cor[, recall.cnloh := ifelse(n.cnloh.gs == 0, 1, recall.cnloh)]
            cn.cor[, precision.cnloh := ifelse(n.cnloh.gs == 0 & fp == 0, 1, precision.cnloh)]

            cn.cor[, f1.cnloh := ifelse(precision.cnloh > 0 & recall.cnloh > 0,
                                        2 * precision.cnloh * recall.cnloh / (precision.cnloh + recall.cnloh),
                                        0)]

            cnloh.dt = as.data.table(tiles)[(gs) | (get(method))]
        }
        if (cnloh.dt[, .N]) {
            cnloh.dt[, pair := id]
        }
        cn.cor[, pair := id]
        out = list(cnloh.res = cn.cor, cnloh.gr = cnloh.dt)
    } else {
        out = compare_cncp(gr = method.cnloh,
                           gs = gs.cnloh,
                           max.dist = max.dist,
                           verbose = verbose)

        out.dt = out$cncp.res
        out.dt$pair = id

        out.gr = out$cncp.gr
        if (length(out.gr)) {
            out.gr$pair = id
        }

        out = list(cnloh.res = out.dt, cnloh.gr = out.gr)
    }

    return(out)
}
                          

#' @name cncp_analysis
#' @title cncp_analysis
#'
#' @description
#'
#' This is a wrapper function for CN change point analysis.
#' It compares a CN segmentation of a genome to a set of gold-standard break points.
#'
#' @param seg (character) path to seg file (can be BigWig, GRanges, or gGraph)
#' @param gs (character) path to gold standard SVs (as GRangesList)
#' @param method (character) analysis method (one of jabba, ascat, sequenza, default jabba)
#' @param max.dist (numeric) bp threshold under which a called CNCP is considered to overlap a gold-standard CNCP
#' @param id (character) sample ID
#' @param simplify (logical) simplify the input segments when calling CNCP?
#' @param loose.only (logical) only consider the loose ends
#' @param mask (character) masked ranges (ignore breakpoints in this region)
#' @param verbose (logical) default FALSE
#'
#' @return list with entries:
#' - $cncp.gr: GRanges of CN change points in the segmentation + distance to nearest gold-standard point
#' - $cncp.res: containing F1 score, precision, and recall of CN change points
cncp_analysis = function(seg = NULL,
                         gs = "~/projects/JaBbA_pipeline/db/hcc1954.sv.gs.rds",
                         method = "jabba",
                         id = "sample",
                         max.dist = 1e4,
                         mask = "~/projects/gGnome/files/zc_stash/maskA_re.rds",
                         field = "cn",
                         simplify = FALSE,
                         loose.only = FALSE,
                         ignore.strand = FALSE,
                         verbose = FALSE) {

    if (verbose) {
        message("reading gold standard segments")
    }

    ## let's allow input as either GRanges or GRangesList
    gs.tmp = readRDS(gs)
    if (inherits(gs.tmp, "GRanges")) {
        gs.cncp = gs.tmp
    } else {
        gs.cncp = stack(readRDS(gs))
    }

    if (grepl(pattern = "ascat", x = method, ignore.case = TRUE)) {
        method.gr = grab_segs(bw = seg, field = "score",
                              simplify = simplify,
                              verbose = verbose)
    } else if (grepl(pattern = "jabba", x = method, ignore.case = TRUE)) {
        if (!loose.only) {
            ## method.gr = grab_segs(gg = seg,
            ##                       simplify = simplify,
            ##                       verbose = verbose)
            gg = gG(jabba = seg)
            loose.gr = gg$loose %Q% ((terminal == FALSE) & (cn > 0))
            jj.gr = unlist(gg$junctions[type == "ALT" & cn > 0]$grl)
            method.gr = grbind(loose.gr, jj.gr)
        } else {
            ## if only considering loose ends
            ## then just use $loose function of gGnome
            ## and grab the non-terminal loose ends
            gg = gG(jabba = seg)
            method.gr = gg$loose %Q% (terminal == FALSE & cn > 0)
        }
    } else if (grepl(pattern = "sequenza", x = method, ignore.case = TRUE)) {
        method.gr = grab_segs(csv = seg, field = "CNt",
                              simplify = simplify,
                              verbose = verbose)
    } else if (method == "gr") {
        method.gr = grab_segs(gr = seg, field = field, simplify = TRUE, verbose = verbose)
    } else if (method == "df") {
        method.gr = grab_segs(df = seg, field = field, simplify = TRUE, verbose = verbose)
    } else if (method == "csv") {
        method.gr = grab_segs(csv = seg, field = field, simplify = TRUE, verbose = verbose)
    } else if (method == "grl") {
        grl = readRDS(seg)
        if (length(grl)) {
            ## make these breakends unique
            method.gr = stack(grl)
            method.gr = reduce(gr.stripstrand(method.gr) + max.dist)
        } else {
            method.gr = GRanges()
        }
    } else {
        stop("method not implemented yet")
    }

    if (file.exists(mask) && file.info(mask)$size) {
        mask.gr = readRDS(mask)
    } else {
        mask.gr = GRanges(seqlengths = seqlengths(method.gr))
    }

    if (loose.only | grepl(pattern = "jabba", x = method, ignore.case = TRUE)) {
        method.cncp = method.gr
    } else if (method == "grl") {
          method.cncp = method.gr
    } else {
        method.cncp = grab_cncp(gr = method.gr, verbose = verbose, include.na = FALSE, ignore.strand = FALSE)
    }

    ## remove masked breakpoints
    method.cncp = method.cncp %Q% ((method.cncp %O% mask.gr) < 0.9)
    gs.cncp = gs.cncp %Q% (!(gs.cncp %^% mask.gr))

    if (ignore.strand) {
        method.cncp = gr.stripstrand(method.cncp)
        gs.cncp = gr.stripstrand(gs.cncp)
    }

    out = compare_cncp(gr = method.cncp,
                       gs = gs.cncp,
                       max.dist = max.dist,
                       verbose = verbose,
                       ignore.strand = ignore.strand)

    ## cncp.res = bp.analysis(method.cncp, gs.cncp, pad = max.dist)
    out.gr = out$cncp.gr
    if (length(out.gr)) {
        out.gr$pair = id
    }

    out.dt = out$cncp.res
    out.dt$pair = id

    ## method.cncp$pair = id
    ## cncp.res[, pair := id]

    return(list(cncp.gr = out.gr, cncp.res = out.dt))
}

#' @name compare_jabba_cn
#' @title compare_jabba_cn
#'
#' @description
#'
#' This is a function that compares segment CNs in two jabba outputs
#'
#' @param gg1 (character) path to JaBbA 1
#' @param gg2 (character) path to JaBbA 2
#' @param field1 (character) default cn
#' @param field2 (character) default cn
#' @param tile.width (numeric) resolution at which to tile the genome for correlation. Default 1e4.
#' @param min.width (numeric) min tile width to consider for RMSE, etc.
#' @param id (character) sample ID
#' @param verbose (logical) default FALSE
#'
#' @return list with entries $cncp.gr $cncp.res
compare_jabba_cn = function(gg1 = NULL, gg2 = NULL,
                            field1 = "cn", field2 = "cn",
                            id = "sample",
                            tile.width = 1e4,
                            min.width = 1e4,
                            verbose = FALSE) {

    if (verbose) {
        message("reading JaBbAs")
    }

    gg1.gr = grab_segs(gg = gg1, field = field1, simplify = TRUE, verbose = verbose)
    gg2.gr = grab_segs(gg = gg2, field = field2, simplify = TRUE, verbose = verbose)

    ## remove centromeric regions
    gg1.gr = remove_centromeres(gg1.gr)
    gg2.gr = remove_centromeres(gg2.gr)

    tiles = gr.tile(gg1.gr, width = tile.width)
    cn.comp = compare_cn(gr = gg1.gr, tiles = tiles, nm = "gg1", verbose = verbose)
    cn.comp = compare_cn(gr = gg2.gr, tiles = cn.comp, nm = "gg2", verbose = verbose)
    cn.cor = compute_cn_correlation(gr = cn.comp,
                                    method = "gg1",
                                    gs = "gg2",
                                    min.width = min.width,
                                    verbose = verbose)

    cn.comp$pair = id
    cn.cor[, pair := id]
    return(list(cn.res = cn.cor, cn.seg = cn.comp))
}

#' @name cn_analysis
#' @title cn_analysis
#'
#' @description
#'
#' This is a function that compares segment CNs to a gold standard.
#'
#' @param seg (character) path to seg file
#' @param field (character) field(s) containing CN (or allelic CN) values
#' @param gr (GRanges) allow GRanges input directly
#' @param gs (character) path to gold standard segs
#' @param gs.field (character) CN/allelic CN fields in the gold standard
#' @param gs.method (character) data type of gold standard (gg or gr)
#' @param method (character) analysis method (one of jabba, ascat)
#' @param id (character) sample ID
#' @param tile.width (numeric) resolution at which to tile the genome for correlation. Default 1e4.
#' @param min.width (numeric) min tile width to consider for RMSE, etc.
#' @param mask (character) path to coverage mask
#' @param allelic (logical) run allelic (as opposed to total CN) analysis ? default FALSE
#' @param cnloh (logical) CNLOH analysis? (e.g. get precision, recall, F1 for cnloh detection). default FALSE
#' @param loh (logical) LOH analysis? (e.g. get precision, recall, F1 for cnloh detection). default FALSE
#' @param cnv (logical) focus analysis exclusively on CNV regions (as opposed to genome-wide) (default FALSE)
#' @param verbose (logical) default FALSE
#'
#' @return list with entries $cn.gr $cn.res
cn_analysis = function(seg = NULL,
                       field = "cn",
                       gr = NULL,
                       gs = "~/projects/JaBbA_pipeline/db/hcc1954.cn.array.new.rds",
                       gs.field = "cn",
                       gs.method = "gr",
                       method = "jabba",
                       id = "sample",
                       tile.width = 1e4,
                       min.width = 1e3,
                       mask = "~/projects/gGnome/files/zc_stash/maskA_re.rds",
                       allelic = FALSE,
                       cnloh = FALSE,
                       loh = FALSE,
                       cnv = FALSE,
                       verbose = FALSE) {

    if (verbose) {
        message("reading gold standard segments")
    }

    if (cnloh | loh) {
        if (gs.method == "gr") {
            gs.gr = grab_allelic_segs(gr = gs,
                                      field = gs.field,
                                      simplify = TRUE,
                                      cnloh = cnloh,
                                      loh = loh,
                                      verbose = verbose)
        } else if (gs.method == "gg") {
            gs.gr = grab_allelic_segs(gg = gs,
                                      field = gs.field,
                                      simplify = TRUE,
                                      cnloh = cnloh,
                                      loh = loh,
                                      verbose = verbose)
        } else {
            stop("invalid option supplied for gs.method")
        }

        if (is.null(gr)) {

            if (grepl(pattern = "jabba", x = method, ignore.case = TRUE)) {
                method.gr = grab_allelic_segs(gg = seg,
                                              field = c("cn.high", "cn.low"),
                                              simplify = TRUE,
                                              cnloh = cnloh,
                                              loh = loh,
                                              verbose = verbose)
            } else if (method == "gr") {
                method.gr = grab_allelic_segs(gr = seg, field = field,
                                              cnloh = cnloh,
                                              loh = loh,
                                              simplify = TRUE,
                                              verbose = verbose)
            } else if (method == "csv") {
                method.gr = grab_allelic_segs(csv = seg, field = field,
                                              cnloh = cnloh,
                                              loh = loh,
                                              simplify = TRUE,
                                              verbose = verbose)
            } else if (method == "df") {
                method.gr = grab_allelic_segs(df = seg, field = field,
                                              cnloh = cnloh,
                                              loh = loh,
                                              simplify = TRUE,
                                              verbose = verbose)
            } else {
                stop("method not implemented yet")
            }
        } else if (inherits(gr, 'GRanges')) {
            stop("Not implemented yet!")
        } else {
            stop("invalid input for gr")
        }
    } else {

        if (!allelic) {
            if (gs.method == "gr") {
                gs.gr = grab_segs(gr = gs, field = gs.field, simplify = TRUE, verbose = verbose)
            } else if (gs.method == "gg") {
                gs.gr = grab_segs(gg = gs, field = gs.field, simplify = TRUE, verbose = verbose)
            } else {
                stop("invalid option supplied for gs.method")
            }

            if (is.null(gr)) {

                if (grepl(pattern = "ascat", x = method, ignore.case = TRUE)) {
                    method.gr = grab_segs(bw = seg, field = "score",
                                          simplify = TRUE,
                                          verbose = verbose)
                } else if (grepl(pattern = "jabba", x = method, ignore.case = TRUE)) {
                    method.gr = grab_segs(gg = seg, field = "cn",
                                          simplify = TRUE,
                                          verbose = verbose)
                } else if (grepl(pattern = "sequenza", x = method, ignore.case = TRUE)) {
                    method.gr = grab_segs(csv = seg, field = "CNt",
                                          simplify = TRUE,
                                          verbose = verbose)
                } else if (method == "gr") {
                    method.gr = grab_segs(gr = seg, field = field, simplify = TRUE, verbose = verbose)
                } else if (method == "csv") {
                    method.gr = grab_segs(csv = seg, field = field, simplify = TRUE, verbose = verbose)
                } else if (method == "df") {
                    method.gr = grab_segs(df = seg, field = field, simplify = TRUE, verbose = verbose)
                } else {
                    stop("method not implemented yet")
                }
            } else if (inherits(gr, 'GRanges')) {
                if ('cn' %in% names(values(gr))) {
                    method.gr = gr[, 'cn']
                } else {
                    stop('supplied granges missing field cn')
                }
            } else {
                stop("invalid input for gr")
            }
        } else {
            if (gs.method == "gr") {
                gs.gr = grab_allelic_segs(gr = gs, field = gs.field, simplify = TRUE, verbose = verbose)
            } else if (gs.method == "gg") {
                gs.gr = grab_allelic_segs(gg = gs, field = gs.field, simplify = TRUE, verbose = verbose)
            } else {
                stop("invalid option supplied for gs.method")
            }

            if (is.null(gr)) {

                if (grepl(pattern = "jabba", x = method, ignore.case = TRUE)) {
                    method.gr = grab_allelic_segs(gg = seg,
                                                  field = c("cn.high", "cn.low"),
                                                  simplify = TRUE,
                                                  verbose = verbose)
                } else if (method == "gr") {
                    method.gr = grab_allelic_segs(gr = seg, field = field, simplify = TRUE, verbose = verbose)
                } else if (method == "csv") {
                    method.gr = grab_allelic_segs(csv = seg, field = field, simplify = TRUE, verbose = verbose)
                } else if (method == "df") {
                    method.gr = grab_allelic_segs(df = seg, field = field, simplify = TRUE, verbose = verbose)
                } else {
                    stop("method not implemented yet")
                }
            } else if (inherits(gr, 'GRanges')) {
                stop("Not implemented yet!")
                ## if ('cn' %in% names(values(gr))) {
                ##     method.gr = gr[, 'cn']
                ## } else {
                ##     stop('supplied granges missing field cn')
                ## }
            } else {
                stop("invalid input for gr")
            }
        }
    }

    ## remove centromeric regions
    ## method.gr = remove_centromeres(method.gr)

    ## remove tiles overlapping with coverage mask if provided
    if (file.exists(mask) && file.info(mask)$size) {
        mask.gr = readRDS(mask)
    } else {
        mask.gr = GRanges(seqlengths = seqlengths(method.gr))
    }

    ## tile the gold standard genome
    tiles = gr.tile(si2gr(si = hg_seqlengths()), width = tile.width)

    if (allelic) {
        tiles.high = tiles[, c()]
        tiles.low = tiles[, c()]
        values(tiles.high)[, "allele"] = "major"
        values(tiles.low)[, "allele"] = "major"
        tiles = c(tiles.high, tiles.low)
    }

    ## get the mask fraction of each tile
    values(tiles)[, "mask"] = tiles %O% mask.gr

    ## only pick segments greater than minimum width
    names(method.gr) = NULL
    method.gr = method.gr %Q% (width(method.gr) > min.width)
    gs.gr = gs.gr %Q% (width(gs.gr) > min.width)

    ## if cnv-specific analysis is desired, these need to be identified
    if (cnv) {
        if (verbose) { message("Identifying regions with CNVs") }
        if (allelic) {
            gs.dt = as.data.table(gs.gr)[, mean.cn := round(weighted.mean(x = cn,
                                                                          w = (end - start),
                                                                          na.rm = TRUE)),
                                         by = .(seqnames, allele)]
        } else {
            gs.dt = as.data.table(gs.gr)[, mean.cn := round(weighted.mean(x = cn,
                                                                          w = (end - start),
                                                                          na.rm = TRUE)),
                                         by = .(seqnames)]
        }
        cnv.gr = (gs.gr %Q% (abs(cn - gs.dt[, mean.cn]) >= 1))
        tiles = tiles %Q% (tiles %^% cnv.gr)
    } 


    ## overlap the gold standard and query ranges against tiles
    if (!(cnloh | loh)) {
        if (allelic) {
            gs.ov.dt = gr.findoverlaps(query = tiles, subject = gs.gr,
                                       return.type = "data.table", by = "allele")
            method.ov.dt = gr.findoverlaps(query = tiles, subject = method.gr,
                                           return.type = "data.table", by = "allele")
        } else {
            gs.ov.dt = gr.findoverlaps(query = tiles, subject = gs.gr, return.type = "data.table")
            method.ov.dt = gr.findoverlaps(query = tiles, subject = method.gr, return.type = "data.table")
        }

        gs.ov.dt[, gs := values(gs.gr)[, "cn"][subject.id]]
        method.ov.dt[, method := values(method.gr)[, "cn"][subject.id]]

        gs.dt = gs.ov.dt[, .(gs = weighted.mean(x = gs, w = end - start, na.rm = TRUE)), by = query.id]
        method.dt = method.ov.dt[, .(method = weighted.mean(x = method, w = end - start, na.rm = TRUE)), by = query.id]

        values(tiles)[, "gs"] = gs.dt[, gs][match(1:length(tiles), gs.dt[, query.id])]
        values(tiles)[, method] = method.dt[, method][match(1:length(tiles), method.dt[, query.id])]
    } else {
        values(tiles)[, "gs"] = tiles %^% gs.gr
        values(tiles)[, method] = tiles %^% method.gr
    }
    
    tiles = tiles %Q% (mask < 0.8)##(tiles %O% mask.gr) < 0.8)
    tiles = tiles %Q% (as.character(seqnames(tiles)) != "M" & as.character(seqnames(tiles)) != "Y")
    tiles = tiles %Q% ((!is.na(gs)) & (!is.na(values(tiles)[, method])))

    if (!(cnloh | loh)) {
        cn.cor = compute_cn_correlation(gr = tiles,
                                        method = method,
                                        gs = "gs",
                                        min.width = 0, ## keep all tiles
                                        verbose = verbose)
    } else {
        tiles = tiles %Q% (as.character(seqnames(tiles)) != "X")
        if (cnloh) {
            cn.cor = data.table(pair = id,
                                n.cnloh = sum(values(tiles)[, method]),
                                n.cnloh.gs = sum(values(tiles)[, "gs"]),
                                tp = sum(values(tiles)[, method] & values(tiles)[, "gs"]),
                                fn = sum(values(tiles)[, "gs"] & (!values(tiles)[, method])),
                                fp = sum(values(tiles)[, method] & (!values(tiles)[, "gs"])),
                                tn = sum((!values(tiles)[, method]) & (!values(tiles)[, "gs"])))

            
            cn.cor[, precision.cnloh := ifelse(tp > 0, tp / (tp + fp), 0)]
            cn.cor[, recall.cnloh := ifelse(tp > 0, tp / (tp + fn), 0)]

            ## set to NA if there are simply not any CNLOH
            cn.cor[, recall.cnloh := ifelse(n.cnloh.gs == 0, 1, recall.cnloh)]
            cn.cor[, precision.cnloh := ifelse(n.cnloh.gs == 0 & fp == 0, 1, precision.cnloh)]

            cn.cor[, f1.cnloh := ifelse(precision.cnloh > 0 & recall.cnloh > 0,
                                        2 * precision.cnloh * recall.cnloh / (precision.cnloh + recall.cnloh),
                                        0)]
        } else {
            cn.cor = data.table(pair = id,
                                n.loh = sum(values(tiles)[, method]),
                                n.loh.gs = sum(values(tiles)[, "gs"]),
                                tp = sum(values(tiles)[, method] & values(tiles)[, "gs"]),
                                fn = sum(values(tiles)[, "gs"] & (!values(tiles)[, method])),
                                fp = sum(values(tiles)[, method] & (!values(tiles)[, "gs"])),
                                tn = sum((!values(tiles)[, method]) & (!values(tiles)[, "gs"])))

            cn.cor[, precision.loh := ifelse(tp > 0, tp / (tp + fp), 0)]
            cn.cor[, recall.loh := ifelse(tp > 0, tp / (tp + fn), 0)]

            ## set to NA if there are simply not any CNLOH
            cn.cor[, recall.loh := ifelse(n.loh.gs == 0, 1, recall.loh)]
            cn.cor[, precision.loh := ifelse(n.loh.gs == 0 & fp == 0, 1, precision.loh)]

            cn.cor[, f1.loh := ifelse(precision.loh > 0 & recall.loh > 0,
                                      2 * precision.loh * recall.loh / (precision.loh + recall.loh),
                                      0)]
        }
    }

    ## add total non-NA width?
    values(tiles)[, "pair"] = id
    cn.cor[, pair := id]
    cn.cor[, query.width := sum(width(tiles[!is.na(values(tiles)[, method])]))]
    return(list(cn.res = cn.cor, cn.seg = tiles))
}

#' @name check_class
#' @title check_class
#'
#' @description
#'
#' Test if object is an instance of some class
#' But also if object is null/NA, which breaks inherits()
#'
#' @param obj
#' @param class (character)
#'
#' @return logical, TRUE if obj is an instance of class
check_class = function(obj = NULL, class = character()) {

    if (is.null(obj)) {
        return(FALSE)
    }

    return(inherits(x = obj, what = class))
}

#' @name check_file
#' @title check_file
#'
#' @description
#'
#' Check if file supplied is nonempty
#'
#' @param fn
#'
#' @return logical, TRUE if file exists and is nonempty
check_file = function(fn = NULL) {
    if (is.null(fn)) {
        return(FALSE)
    }

    if (!file.exists(fn)) {
        return(FALSE)
    }

    if (!file.info(fn)$size) {
        return(FALSE)
    }

    return(TRUE)
}

#' @name grab_centromeres
#' @title grab_centromeres
#'
#' @param pad centromeres?
#' 
#' @return GRanges of centromeres
grab_centromeres = function(pad = 0) {
    bands.td = gTrack::karyogram()
    bands.td$height=5
    bands = bands.td@data
    bands = grl.unlist(do.call(`GRangesList`, bands))
    centromeres.gr = reduce((bands %Q% (stain == "acen")) + pad)
    return(centromeres.gr)
}

#' @name remove_centromeres
#' @title remove_centromeres
#'
#' @description
#'
#' Remove segments in an input GRanges overlapping centromeres, preserving metadata.
#'
#' @param segs (GRanges)
#'
#' @return input GRanges with ranges overlapping centromeres removed.
remove_centromeres = function(segs, pad = 0) {

    ## read centromeric regions from gTrack
    bands.td = gTrack::karyogram()
    bands.td$height=5
    bands = bands.td@data
    bands = grl.unlist(do.call(`GRangesList`, bands))
    eligible = bands %Q% (stain != "acen") ## excluding CENTROMERE
    eligible = gr.reduce(eligible + 1e4) - 1e4
    starts = start(eligible)
    ends = end(eligible)
    sl = seqlengths(eligible)
    starts[starts > 1] = starts[starts > 1] + pad
    ends[ends < sl[as.character(seqnames(eligible))]] = ends[ends < sl[as.character(seqnames(eligible))]] - pad
    new.eligible = GRanges(seqnames = as.character(seqnames(eligible)),
                       ranges = IRanges(start = starts, end = ends))
    segs = segs %*% new.eligible[, c()]
    return(segs)
}

#' @name overlap_breakpoints_with_jabba
#' @title overlap_breakpoints_with_jabba
#'
#' @description
#'
#' Overlap a GRanges of loose ends with
#' the breakpoints in a JaBbA graph
#'
#' Determine whether each loose end overlaps a junction, a loose end,
#' or overlaps an (optional) coverage mask.
#'
#' @param gr (GRanges) of loose ends
#' @param jabba_rds (character) path to JaBbA output
#' @param mask (character) path to coverage mask
#' @param min.footprint (numeric) min SV size (default 1e4)
#' @param bp.pad (numeric) distance in bp from junction/loose end/mask (default 1e4)
#' @param mask.pad (numeric) distance in bp from junction/loose end/mask (default 1e4)
#' @param verbose (logical)
#'
#' @return GRanges with same length as input with three metadata columns:
#' $junction (logical)
#' $loose (logical)
#' $mask (logical)
overlap_breakpoints_with_jabba = function(gr = GRanges(), jabba_rds = NA_character_,
                                          mask = "~/projects/gGnome/files/new.mask.rds",
                                          min.footprint = 0,
                                          bp.pad = 1e4,
                                          mask.pad = 1e4,
                                          verbose = FALSE) {
    if (!inherits(gr, 'GRanges')) {
        stop("gr must be GRanges")
    }
    if (!file.exists(jabba_rds) || !file.info(jabba_rds)$size) {
        stop("JaBbA output not valid file")
    }
    if (!length(gr)) {
        if (verbose) {
            message("Returning empty GRanges")
        }
        return(gr)
    }

    if (verbose) {
        message("Copying input to avoid mutation")
    }

    gr = copy(gr)

    if (verbose) {
        message("Loading JaBbA input")
    }

    gg = gG(jabba = jabba_rds)

    if (verbose) {
        message("Grabbing loose ends from JaBbA")
    }
    
    gg.loose.gr = gg_loose_end(jabba_rds = jabba_rds, return.type = 'GRanges')

    if (verbose) {
        message("Grabbing junction breakpoints from JaBbA")
    }
    bp.gr = GRanges()
    if (length(gg$junctions)) {
        all.alt.junctions = gg$junctions[type == "ALT" & cn > 0]
    } else {
        all.alt.junctions = jJ()
    }
    if (length(all.alt.junctions)) {
        large.sv.ix = which(all.alt.junctions$span > min.footprint)
        if (length(large.sv.ix)) {
            bp.gr = all.alt.junctions[large.sv.ix]$breakpoints
        }
    }

    mask.gr = GRanges()
    if (file.exists(mask) && file.info(mask)$size) {
        if (verbose) {
            message("Reading coverage mask")
        }
        mask.gr = readRDS(mask)
    }

    if (length(gg.loose.gr)) {
        if (verbose) {
            message("Overlapping with JaBbA loose ends")
        }
        gr$loose = gr %^^% (gr.stripstrand(gg.loose.gr) + bp.pad)
    } else {
        if (verbose) {
            message("No unmasked loose ends in supplied JaBbA")
        }
        gr$loose = FALSE
    }

    if (length(bp.gr)) {
        if (verbose) {
            message("Overlapping with JaBbA ALT junctions")
        }
        gr$junction = gr %^^% (gr.stripstrand(bp.gr) + bp.pad)
    } else {
        if (verbose) {
            message("No CN > 0 ALT junctions with required span")
        }
        gr$junction = FALSE
    }

    if (length(mask.gr)) {
        if (verbose) {
            message("Overlapping with coverage mask")
        }
        gr$mask = gr %^% (gr.stripstrand(mask.gr) + mask.pad)
    } else {
        gr$mask = FALSE
    }
    return(gr)
}

#' @name grab_loose_from_seg
#' @title grab_loose_from_seg
#' 
#' @description
#'
#' get the location and direction of copy number increases
#' e.g. if the neighboring node on the left is lower, then strand is +
#' if the neighboring node on the right is lower, then strand is -
#' (opposite direction from paper, but consistent with loose ends from JaBbA)
#'
#' @param gr (GRanges) with metadata field $cn
#' @param min.width (numeric) minimum SV width (default 0?)
#' @param max.gap (numeric) maximum gap between non-NA segments
#' @param verbose (logical)
#'
#' @return GRanges each with width 1 giving CN change points
grab_loose_from_seg = function(gr = NULL,
                               min.width = 0,
                               max.gap = 1e5,
                               return.type = "data.table",
                               verbose = FALSE) {

    if (is.null(gr$cn)) {
        stop("gr missing $cn metadata")
    }

    if (verbose) {
        message("Sorting gr...")
    }

    gr = gr %Q% (!is.na(cn))
    ## gr = sortSeqlevels(gr)
    ## gr = sort(gr)

    dt = as.data.table(gr)
    dt = dt[order(start, decreasing = FALSE), .SD, by = seqnames]

    dt[, previous.cn := data.table::shift(cn, type = "lag")]
    dt[, next.cn := data.table::shift(cn, type = "lead")]
    dt[, previous.seqnames := data.table::shift(seqnames, type = "lag")]
    dt[, next.seqnames := data.table::shift(seqnames, type = "lead")]

    dt[seqnames == previous.seqnames & cn > previous.cn, bp.start := TRUE]
    dt[seqnames == next.seqnames & cn > next.cn, bp.end := TRUE]

    ## annotate with width
    dt[, previous.width := data.table::shift(width, type = 'lag')]
    dt[, next.width := data.table::shift(width, type = 'lead')]
    dt[, bp.width := pmin(width, previous.width)]
    dt[, bp.width := pmin(width, next.width)]

    ## make sure that not preceded/followed by a big NA segment
    dt[, previous.end := data.table::shift(end, type = "lag")]
    dt[, previous.gap := start - previous.end]
    dt[, next.start := data.table::shift(start, type = "lead")]
    dt[, next.gap := next.start - end]

    bps = rbind(dt[(bp.start) & bp.width > min.width & previous.gap < max.gap,
                   .(seqnames, start, end = start, width = 1, strand = "+")],
                 dt[(bp.end) & bp.width > min.width & next.gap < max.gap,
                    .(seqnames, start = end, end, width = 1, strand = "-")])

    if (return.type == "data.table") {
        return(bps)
    }
    return(dt2gr(bps, seqlengths = seqlengths(gr)))
}

#' @name grab_allelic_segs
#' @name grab_allelic_segs
#'
#' @description
#'
#' get allelic copy number as a melted data table
#'
#' @param gg (character) path to jabba_rds
#' @param gr (character) path to GRanges (.rds file)
#' @param csv (character) path to delimited file (.txt)
#' @param df (character) path to data.frame object (.rds)
#' @param field (character) fields of major and minor alleles
#' @param cnloh (logical) return just CNLOH loci? default FALSE
#' @param loh (logical) return just LOH loci? default FALSE
#' @param simplify (logical) simplify segs by cn?
#' @param melt (logical) melt ranges? if not returns with fields cn.high and cn.low
#' @param verbose (logical) print stuff or nah
#'
#' @return (melted) GRanges with values haplotype and cn if CNLOH is FALSE and just cnloh segs if TRUE
grab_allelic_segs = function(gg = "/dev/null",
                             gr = "/dev/null",
                             csv = "/dev/null",
                             df = "/dev/null",
                             field = c("cn.high", "cn.low"),
                             cnloh = FALSE,
                             loh = FALSE,
                             simplify = FALSE,
                             melt = TRUE,
                             verbose = FALSE) {

    if (length(field) != 2) {
        stop("Field must be a character vector of length 2")
    }
    
    if (check_file(gg)) {

        if (verbose) {
            message("Reading gGraph")
        }
        this.gg = gG(jabba = gg)
        segs = gr.stripstrand(this.gg$nodes$gr[, field])
    } else if (check_file(gr)) {

        if (verbose) {
            message("Reading GRanges")
        }
        
        this.gr = readRDS(gr)
        segs = gr.stripstrand(this.gr[, field])
    } else if (check_file(df)) {

        if (verbose) {
            message("Grabbing dataframe")
        }

        this.df = readRDS(df)
        segs = GRanges(seqnames = this.df[, 1],
                       ranges = IRanges(start = this.df[, 2], end = this.df[, 3]))

        values(segs) = this.df[, field]
        
    } else if (check_file(csv)) {

        if (verbose) {
            message("Reading csv")
        }

        dt = fread(csv)
        sn = grep("^([Cc]hr*)|(seqnames)", colnames(dt), value = T)[1]
        st = grep("^[Ss]tart", colnames(dt), value = T)[1]
        ed = grep("^[Ee]nd", colnames(dt), value = T)[1]
        setnames(dt, c(sn, st, ed), c("seqnames", "start", "end"))
        dt[, seqnames := as.character(seqnames)]
        dt[seqnames == "23", seqnames := "X"]
        dt[seqnames == "24", seqnames := "Y"]
        this.csv = dt2gr(dt[end > start,])
        segs = gr.stripstrand(this.csv[, field])
        
    } else {
        stop("no valid file supplied")
    }

    ## manually set cn.low and cn.high so that they represent major and minor
    values(segs)[, "cn.high"] = pmax(values(segs)[, field[1]], values(segs)[, field[2]])
    values(segs)[, "cn.low"] = pmin(values(segs)[, field[1]], values(segs)[, field[2]])

    ## remove any NA ranges (this can happen as sometimes there are places
    ## with LogR but no SNVs
    segs = segs %Q% (!is.na(cn.high)) %Q% (!is.na(cn.low))

    ## identify CNLOH segments and return if that is desired
    if (cnloh) {
        if (verbose) { message("Identifying CNLOH segments") }
        segs = segs %Q% (cn.high == 2 & cn.low == 0)
        return(segs)
    } else if (loh) {
        if (verbose) { message("Identifying LOH segments") }
        segs = segs %Q% (cn.high > 0 & cn.low == 0)
        return(segs)
    }

    if (melt) {

        ## separate into major and minor, then melt ranges
        major.segs = segs[, "cn.high"]
        names(values(major.segs)) = "cn"
        values(major.segs)[, "allele"] = "major"

        minor.segs = segs[, "cn.low"]
        names(values(minor.segs)) = "cn"
        values(minor.segs)[, "allele"] = "minor"

        melted.segs = c(major.segs, minor.segs)
        names(melted.segs) = NULL
        return(melted.segs)
    }

    return(segs)
    
}
    
    



#' @name grab_segs
#' @title grab_segs
#'
#' @description
#'
#' read CN segmentation from supplied gGraph, GRanges, bigWig, or delimited file coercible to GRanges
#' 
#' @param gg (character) path to jabba_rds
#' @param gr (character) path to GRanges (.rds file)
#' @param bw (character) path to bigWig (.bw)
#' @param csv (character) path to delimited file coercible to GRanges
#' @param df (character) path to data.frame object
#' @param simplify (logical) default TRUE. Merge adjacent ranges with the same copy number?
#' @param verbose (logical)
#' 
#' @return GRanges with metadata field $cn
grab_segs = function(gg = "/dev/null",
                     gr = "/dev/null",
                     bw = "/dev/null",
                     csv = "/dev/null",
                     df = "/dev/null",
                     field = "cn",
                     simplify = TRUE,
                     verbose = FALSE) {

    if (check_file(gg)) {

        if (verbose) {
            message("Reading gGraph")
        }
        this.gg = gG(jabba = gg)
        segs = gr.stripstrand(this.gg$nodes$gr[, field])
        
    } else if (check_file(gr)) {

        if (verbose) {
            message("Reading GRanges")
        }
        
        this.gr = readRDS(gr)
        segs = gr.stripstrand(this.gr[, field])
        values(segs)[, "cn"] = copy(values(segs)[, field])
        
    } else if (check_file(df)) {

        if (verbose) {
            message("Grabbing dataframe")
        }

        this.df = readRDS(df)
        segs = GRanges(seqnames = this.df[, 1],
                        ranges = IRanges(start = this.df[, 2], end = this.df[, 3]),
                        cn = this.df[, field])
        
    } else if (check_file(bw)) {

        if (verbose) {
            message("Reading bigWig")
        }
        
        this.bw = import.bw(bw)
        segs = gr.stripstrand(this.bw[, field])
    } else if (check_file(csv)) {

        if (verbose) {
            message("Reading csv")
        }

        dt = fread(csv)
        sn = grep("^([Cc]hr*)|(seqnames)", colnames(dt), value = T)[1]
        st = grep("^[Ss]tart", colnames(dt), value = T)[1]
        ed = grep("^[Ee]nd", colnames(dt), value = T)[1]
        setnames(dt, c(sn, st, ed), c("seqnames", "start", "end"))
        dt[, seqnames := as.character(seqnames)]
        dt[seqnames == "23", seqnames := "X"]
        dt[seqnames == "24", seqnames := "Y"]
        this.csv = dt2gr(dt[end > start,])
        segs = gr.stripstrand(this.csv[, field])
        
    } else {
        stop("no valid file supplied")
    }

    if (is.null(segs$cn)) {
        segs$cn = values(segs)[[field]]
        values(segs)[[field]] = NULL
    }

    if (verbose) {
        message("excluding centromeric segments")
    }

    ## fill in gaps for segments that are not perfectly adjacent
    segs.gaps = gaps(segs)
    segs.gaps = segs.gaps %Q% (strand(segs.gaps)=="*")
    if (length(segs.gaps)) {
        segs.gaps$cn = NA
        segs = c(segs, segs.gaps)
    }

    if (!simplify) {
        return(segs)
    }

    if (verbose) {
        message("Simplifying segments")
    }

    new.segs = simplify_segs(gr = segs[, "cn"], sort.gr = TRUE, verbose = verbose)

    return(new.segs)
}

#' @name simplify_segs
#' @title simplify_segs
#' 
#' @description
#'
#' simplify CN segmentation of GRanges object. Must have field $cn.
#'
#' @param gr (GRanges)
#' @param sort.gr (logical) default TRUE
#' @param verbose (logical) default FALSE
#'
#' @return simplifed GRanges with NA segments removed
simplify_segs = function(gr = NULL, sort.gr = TRUE, verbose = FALSE) {

    if (is.null(gr)) {
        stop("gr must be supplied")
    }

    if (!inherits(gr, 'GRanges')) {
        stop("gr must be GRanges")
    }

    if (is.null(gr$cn)) {
        stop("gr must have field $cn")
    }

    new.segs = unlist(reduce(split(gr, ~ cn)))
    new.segs$cn = as.numeric(names(new.segs))
    names(new.segs) = NULL

    ## reorder
    if (!sort.gr) {
        return(new.segs)
    }

    if (verbose) {
        message("sorting simplified ranges...")
    }
    
    new.segs = GenomeInfoDb::sortSeqlevels(new.segs)
    new.segs = sort(new.segs)

    return(new.segs)
}

#' @name grab_cncp
#' @title grab_cncp
#' 
#' @description
#'
#' copy number changepoints from simplifed CN segmentation
#'
#' @param gr (GRanges) with metadata field $cn
#' @param include.na (logical) include segments with NA CN as CN changepoint
#' @param verbose (logical)
#'
#' @return GRanges each with width 1 giving CN change points
grab_cncp = function(gr = NULL,
                     ignore.strand = FALSE,
                     include.na = TRUE,
                     verbose = FALSE) {

    if (is.null(gr$cn)) {
        stop("gr missing $cn metadata")
    }

    if (verbose) {
        message("Sorting gr...")
    }

    ## gr = GenomeInfoDb::sortSeqlevels(gr)
    ## gr = sort(gr)
    dt = as.data.table(gr)[order(start)][order(seqnames)]

    if (ignore.strand) {
        dt[, shifted.cn := data.table::shift(cn)]
        dt[, shifted.seqnames := data.table::shift(seqnames)]

        if (include.na) {
            cncp.dt = dt[(shifted.cn != cn |
                          is.na(shifted.cn) |
                          is.na(cn)) &
                         shifted.seqnames == seqnames,]
        } else {
            cncp.dt = dt[shifted.cn != cn &
                         !is.na(cn) &
                         !is.na(shifted.cn) &
                         shifted.seqnames == seqnames,]
        }
        cncp.gr = GRanges(seqnames = cncp.dt[, seqnames],
                          ranges = IRanges(start = cncp.dt[, start],
                                           width = 1),
                          cn = cncp.dt[, cn])
    } else {

        if (!include.na) {
            dt = dt[!is.na(cn)]
        }
        
        dt[, ":="(prev.cn = data.table::shift(cn, type = "lag"),
                  next.cn = data.table::shift(cn, type = "lead"),
                  prev.seqnames = data.table::shift(seqnames, type = "lag"),
                  next.seqnames = data.table::shift(seqnames, type = "lead"))]

        if (include.na) {
            cncp.dt = rbind(dt[(cn > prev.cn | is.na(cn)) & (seqnames == prev.seqnames),
                               .(seqnames, start, end = start, strand = "+", cn = cn - prev.cn)],
                            dt[(cn > next.cn | is.na(cn)) & (seqnames == next.seqnames),
                               .(seqnames, start = end, end, strand = "-", cn = cn - next.cn)])
        } else {
            cncp.dt = rbind(dt[(cn > prev.cn | is.na(cn)) & (seqnames == prev.seqnames),
                               .(seqnames, start, end = start, strand = "+", cn = cn - prev.cn)],
                            dt[(cn > next.cn | is.na(cn)) & (seqnames == next.seqnames),
                               .(seqnames, start = end, end, strand = "-", cn = cn - next.cn)])            
        }
        cncp.gr = GRanges(seqnames = cncp.dt[, seqnames],
                          ranges = IRanges(start = cncp.dt[, start],
                                           width = 1),
                          strand = cncp.dt[, strand],
                          cn = cncp.dt[, cn])
    }

    return(cncp.gr)
}

#' @name grab_bad_segs
#' @title grab_bad_segs
#'
#' @description
#'
#' Return GRanges with bad segments (high difference between fields)
#'
#' @param jabba_rds (character) path to jabba output
#' @param gs (character) path to gold standard GRanges
#' @param kag (character) path to karyograph - compare to cn.old
#' @param field1 (character) default cn
#' @param field2 (character) default cn
#' @param bad.thresh (character) numeric, default 1
#' @param mask (character) path to mask
#' @param return.type (character) either GRanges or data.table
#' @param verbose (logical)
#' @param GRanges with columns cn, cn.gs, and diff (or data.table)
grab_bad_segs = function(jabba_rds,
                         gs = NULL,
                         kag = NULL,
                         field1 = "cn",
                         field2 = "cn",
                         bad.thresh = 1,
                         mask = "~/projects/gGnome/files/new.mask.rds",
                         return.type = "data.table",
                         verbose = FALSE) {

    if (!check_file(jabba_rds)) {
        stop("Bad jabba_rds")
    }

    gg = gG(jabba = jabba_rds)
    
    if (check_file(gs)) {
        gs.gr = readRDS(gs)
        if (!inherits(gs.gr, 'GRanges')) {
            stop("must supply GRanges in gs")
        }
    } else if (check_file(kag)) {
        gs.gr = gG(jabba = kag)$nodes$gr[, field2]
    } else {
        stop("either kag or gs must be supplied")
    }

    if (field2 %in% names(values(gs.gr))) {
        gs.gr$cn.gs = values(gs.gr)[[field2]]
    } else {
        stop("not found in gs metadata: ", field2)
    }

    if (check_file(mask)) {
        mask.gr = readRDS(mask)
    } else {
        warning("Invalid mask supplied: ", mask)
        mask.gr = GRanges(seqlengths = seqlengths(gg))
    }

    gs.gr = gr.tile(gs.gr, 1e4) %$% gs.gr[, "cn.gs"]
    gs.gr = gs.gr %Q% (!(gs.gr %^% mask.gr))
    ov = gr.findoverlaps(gg$nodes$gr, gs.gr, return.type = "data.table")
    ov[, cn := gg$nodes$dt[, get(field1)][query.id]]
    ov[, cn.gs := gs.gr$cn.gs[subject.id]]
    ov[, diff := abs(cn - cn.gs)]


    if (return.type == 'GRanges') {
        return(dt2gr())
    }

    return(ov[diff > bad.thresh,])
}
                         

#' @name compare_cn
#' @title compare_cn
#' 
#' @description
#'
#' compare CN of two GRanges objects, which each must have field $cn
#'
#' @param gr (GRanges) test CN
#' @param tiles (GRanges) tiles to evaluate CN fit over
#' @param nm (character) metadata field to store output, default cn
#' @param verbose (logical)
#'
#' @return data.table of overlapping ranges
compare_cn = function(gr = NULL, tiles = NULL, nm = "cn", verbose = FALSE) {

    if (!inherits(tiles, 'GRanges')) {
        stop("tiles must be GRanges")
    }
    
    tiles = copy(tiles)

    ## normalize seqinfo
    if (verbose) {
        message("Normalizing seqinfo")
    }

    sl = hg_seqlengths(chr = FALSE)
    tiles = tiles %Q% (seqnames %in% names(sl))
    seqlevels(tiles) = names(sl)
    tiles = gr.fix(tiles, sl, drop = TRUE)
    
    ## make these consistent
    gr = gr[which(as.character(seqnames(gr)) %in% unique(as.character(seqnames(tiles))))]
    ## seqinfo(gr) = seqinfo(tiles)
    ## seqlengths(gr) = seqlengths(tiles)

    gr.ov = gr.findoverlaps(query = gr, subject = tiles, qcol = "cn",
                            return.type = "data.table")


    gr.ov[, width := end - start]
    gr.ov = gr.ov[!is.na(width) & width > 0, .(cn = weighted.mean(cn, width, na.rm = T)), by = subject.id]
    values(tiles)[[nm]] = gr.ov$cn[match(1:length(tiles), gr.ov[, subject.id])]

    return(tiles)
}
    
#' @name compare_cncp
#' @title compare_cncp
#' 
#' @description
#'
#' compare locations of CN change points given by two GRanges objects
#'
#' @param gr (GRanges) test CN change points
#' @param gs (GRanges) gold standard CN change points
#' @param max.dist (numeric) max allowable distance in BP for a match
#' @param verbose (logical)
#'
#' @return list with cncp.res and cncp.gr
compare_cncp = function(gr = NULL, gs = NULL, max.dist = 1e3, verbose = FALSE, ignore.strand = FALSE) {

    if (verbose) {
        message("Computing pairwise distances between CN change points")
    }

    if (length(gr) & length(gs)) {
        ov = gr.dist(gr1 = gr, gr2 = gs, ignore.strand = ignore.strand)
        ov.dt = data.table(which(ov<max.dist, arr.ind=T))
        setnames(ov.dt, c("row", "col"), c("gr.id", "gs.id"))
        ov.dt[, dist := ov[cbind(gr.id, ov.dt$gs.id)]]
    } else {
        ov.dt = data.table(gr.id = numeric(), gs.id = numeric(), dist = numeric())
    }
    
    if (verbose) {
        message("computing match stats")
    }

    res = data.table(n.cncp = length(gr),
                     n.cncp.gs = length(gs),
                     tp = length(unique(ov.dt[, gr.id])))

    res[, fp := n.cncp - tp]
    res[, fn := n.cncp.gs - length(unique(ov.dt[, gs.id]))]
    res[, precision.cncp := ifelse(tp > 0, tp / (tp + fp), 0)]
    res[, recall.cncp := ifelse(tp > 0, tp / (tp + fn), 0)]
    res[, f1.cncp := ifelse(precision.cncp > 0 & recall.cncp > 0, 2 * precision.cncp * recall.cncp / (precision.cncp + recall.cncp), 0)]

    ## prepare query GRanges and distance to closest gold standard GRanges
    cncp.gr = gr[, c()]

    if (nrow(ov.dt)) {
        closest.dt = ov.dt[, .(gs.id = .SD$gs.id[which.min(.SD$dist)],
                                   dist = min(.SD$dist, na.rm = TRUE)),
                           by = gr.id]
        cncp.gr$closest.gs = closest.dt$gs.id[match(1:length(cncp.gr), closest.dt$gr.id)]
        cncp.gr$dist = closest.dt$dist[match(1:length(cncp.gr), closest.dt$gr.id)]
    } 

    out = list(cncp.res = res,
               cncp.gr = cncp.gr)

    return(out)
}

#' @name compute_cn_correlation
#' @title compute_cn_correlation
#'
#' @description
#'
#' compute Pearson and Spearman correlation and RMSE
#'
#' @param gr (GRanges)
#' @param method (character) metadata field of method results
#' @param gs (character) metadata field of gold standard results
#' @param min.width (numeric) ignore segs shorter than this
#' @param verbose (logical)
#'
#' @return data.table with columns pearson.cn and spearman.cn
compute_cn_correlation = function(gr = NULL,
                                  method = "cn",
                                  gs = "gs",
                                  min.width = 1e3,
                                  verbose = FALSE) {

    if (is.null(values(gr)[[method]])) {
        stop("$method missing")
    }

    if (is.null(values(gr)[[gs]])) {
        stop("$gs missing")
    }

    if (verbose) {
        message("Computing RMSE")
    }

    sel = which(width(gr) >= min.width)
    gr = gr[sel]

    ## remove anything with NA values
    sel = which((!is.na(values(gr)[, method])) & (!is.na(values(gr)[, gs])))
    gr = gr[sel]

    cn.err = values(gr)[[method]] - values(gr)[[gs]]
    sq.err = cn.err^2
    mse = mean(sq.err, na.rm = TRUE)
    rmse = sqrt(mse)

    if (verbose) {
        message("Computing CN correlation")
    }
    
    res = data.table(pearson.cn = cor(values(gr)[[method]],
                                      values(gr)[[gs]],
                                      use = "na.or.complete",
                                      method = "pearson"),
                     spearman.cn = cor(values(gr)[[method]],
                                      values(gr)[[gs]],
                                      use = "na.or.complete",
                                      method = "spearman"),
                     rmse = rmse)

    return(res)
}

#' @name grab.hets.from.maf
#' @title grab.hets.from.maf
#'
#' @description
#'
#' get hets into format needed by phased.binstats from HMF maf_approx
#'
#' @param agt.fname (character)
#' @param min.frac (numeric)
#'
#' @param GRanges
grab.hets.from.maf = function(agt.fname, min.frac = 0.2) {

    if (!file.exists(agt.fname)) {
        stop("invalid file")
    }
    
    gr = readRDS(agt.fname)

    if (!inherits(gr, "GRanges")) {
        stop("rds file must contain GRanges object.")
    }

    dt = as.data.table(gr)
    dt[, major.count := ifelse(alt.count.t > ref.count.t, alt.count.t, ref.count.t)]
    dt[, minor.count := ifelse(alt.count.t > ref.count.t, ref.count.t, alt.count.t)]

    ## out.dt = rbind(dt[, .(seqnames, start, end, count = major.count, allele = "major")],
    ##                dt[, .(seqnames, start, end, count = minor.count, allele = "minor")])

    ## return(dt2gr(out.dt))
    return(dt)
}

#' @name reads2clouds
#' @title reads2clouds
#'
#' @param reads (data.table) must have columns seqnames, start, BX
#' @param thresh (numeric) distance between barcodes for merging to clouds
#' @param read.size (numeric) base pairs of average read?
#'
#' @return data.table representing read clouds
#' modified version of reads to clouds
reads2clouds = function(reads, thresh = 5e3, read.size = 101)
{
    reads = gr2dt(reads)
    setkeyv(reads, c("seqnames", "start"))
    ## make sure end is there, if not add it
    if (is.null(reads$end)) {
        reads[, end := start + read.size]
    }
    reads[, bx.diff := c((start-data.table::shift(end))[-1], NA), by = .(seqnames, BX)]
    reads[, rl := label.runs(bx.diff<thresh | is.na(bx.diff)), by = .(seqnames, BX)]
    reads[, rl.last := data.table::shift(rl), by = .(seqnames, BX)]
    reads[is.na(rl), rl := ifelse(is.na(rl.last), -(1:.N), rl.last)] ## label remaining loners
    reads[, rll := paste(seqnames, BX, rl, sep = '_')]
    reads = reads[, .(start = start[1], end = end[.N], nreads = .N), by = .(seqnames, BX, rll)]
    return(reads)
}

#' @name valid_walks_legacy
#' @title valid_walks_legacy
#'
#' @param junctions (gGnome::Junction)
#' @param pad (numeric) how much to pad the sides of walks (default 2e4)
#' @param max_dist (numeric) default 1 Mbp
#'
#' @return gWalks
#' @export
valid_walks_legacy = function(junctions,
                       pad = 2e4,
                       max_dist = 1e6) {

    bps = stack(junctions$grl)[, c()]

    ## all acceptable permutations of breakends
    breakend.permutations = list(c(1,2,3,4),
                                 c(2,1,3,4),
                                 c(1,2,4,3),
                                 c(2,1,4,3))

    ## track valid permutations
    valid.permutations = c(FALSE, FALSE, FALSE, FALSE)
    for (ix in 1:4) {
        ## consider the i'th permutation
        perm = breakend.permutations[[ix]]
        permuted.bps = bps[, c()][perm]
        ## the contig of the middle two breakends must be the same
        if (as.character(seqnames(permuted.bps)[2]) == as.character(seqnames(permuted.bps)[3])) {
            ## the strand of the middle two breakends must be opposite
            if (as.character(strand(permuted.bps)[2]) != as.character(strand(permuted.bps)[3])) {
                ## the positive stranded breakend must have lower genomic coordinate
                if (as.character(strand(permuted.bps[2])) == "+" &&
                    start(permuted.bps)[2] < start(permuted.bps)[3] &&
                    start(permuted.bps)[3] - start(permuted.bps)[2] < max_dist) {
                    valid.permutations[ix] = TRUE
                } else if (as.character(strand(permuted.bps[2])) == "-" &&
                           start(permuted.bps)[2] > start(permuted.bps)[3] &&
                           start(permuted.bps)[2] - start(permuted.bps)[3] < max_dist) {
                    valid.permutations[ix] = TRUE
                }
            }
        }
    }

    make_seg = function(bp, start = TRUE) {
        strand.bp = as.character(strand(bp))
        if (start) {
            if (strand.bp == "-") {
                gr = GRanges(seqnames = seqnames(bp),
                             ranges = IRanges(start = pmax(1, start(bp) - pad + 1),
                                              end = start(bp)),
                             strand = "+",
                             seqlengths = seqlengths(bp))
            } else if (strand.bp == "+") {
                gr = GRanges(seqnames = seqnames(bp),
                             ranges = IRanges(start = end(bp),
                                              end = pmin(end(bp) + pad - 1,
                                                         seqlengths(bp)[as.character(seqnames(bp))])),
                             strand = "-",
                             seqlengths = seqlengths(bp))
            }
        } else {
            if (strand.bp == "-") {
                gr = GRanges(seqnames = seqnames(bp),
                             ranges = IRanges(start = pmax(1, start(bp) - pad + 1),
                                              end = start(bp)),
                             strand = "-",
                             seqlengths = seqlengths(bp))
            } else if (strand.bp == "+") {
                gr = GRanges(seqnames = seqnames(bp),
                             ranges = IRanges(start = end(bp),
                                              end = pmin(end(bp) + pad - 1,
                                                         seqlengths(bp)[as.character(seqnames(bp))])),
                             strand = "+",
                             seqlengths = seqlengths(bp))
            }
        }
        return(gr)
    }
    
    ## make one-junction walks
    gr1 = make_seg(bps[1, c()], start = TRUE)
    gr2 = make_seg(bps[2, c()], start = FALSE)
    gr3 = make_seg(bps[3, c()], start = TRUE)
    gr4 = make_seg(bps[4, c()], start = FALSE)

    out.grl = GRangesList(c(gr1, gr2),c(gr3, gr4))
    
    ## make two-junction walks if any are valid
    if (any(valid.permutations)) {
        which.valid = which(valid.permutations)
        for (ix in which.valid) {
            permutation = breakend.permutations[[ix]]
            permuted.bps = bps[permutation, c()]
            ## the first segment is the first breakend, extended to the end of the chromosome
            gr1 = make_seg(permuted.bps[1, c()], start = TRUE)
            if (as.character(strand(permuted.bps[2])) == "+") {
                ## the second segment is the middle two breakends
                ## specifically the piece between them
                gr2 = GRanges(seqnames = seqnames(permuted.bps)[2],
                              ranges = IRanges(start = start(permuted.bps[2]),
                                               end = end(permuted.bps[3])),
                              strand = "+",
                              seqlengths = seqlengths(bps))
            } else {
                gr2 = GRanges(seqnames = seqnames(permuted.bps)[2],
                              ranges = IRanges(start = start(permuted.bps[3]),
                                               end = end(permuted.bps[2])),
                              strand = "-",
                              seqlengths = seqlengths(bps))
            }
            ## the third segment is the last breakend, extended to the end of the chromosome
            gr3 = make_seg(permuted.bps[4, c()], start = FALSE)
            gr = c(gr1, gr2, gr3)
            out.grl[[length(out.grl) + 1]] <- gr
        }
    }

    ## disjoin at breakpoints
    out = gW(grl = out.grl, disjoin = TRUE)##gr.stripstrand(bps))

    ## set cis or trans annotation
    annot = ifelse(1:length(out) > 2, "cis", "trans")
    out$set(orientation = annot)

    return(out)
}


#' @name valid_walks
#' @title valid_walks
#'
#' @param junctions (gGnome::Junction)
#' @param pad (numeric) how much to pad the sides of walks (default 2e4)
#' @param max_dist (numeric) default 1 Mbp
#'
#' @return gWalks
#' @export
valid_walks = function(junctions,
                       pad = 2e4,
                       max_dist = 1e6) {

    bps = stack(junctions$grl)[, c()]

    ## all acceptable permutations of breakends
    breakend.permutations = list(c(1,2,3,4),
                                 c(2,1,3,4),
                                 c(1,2,4,3),
                                 c(2,1,4,3))

    intrachromosomal = length(unique(as.character(seqnames(bps)))) == 1

    ## track valid permutations
    valid.permutations = c(FALSE, FALSE, FALSE, FALSE)
    for (ix in 1:4) {
        ## consider the i'th permutation
        perm = breakend.permutations[[ix]]
        permuted.bps = bps[, c()][perm]
        ## the contig of the middle two breakends must be the same
        if (as.character(seqnames(permuted.bps)[2]) == as.character(seqnames(permuted.bps)[3])) {
            ## the strand of the middle two breakends must be opposite
            if (as.character(strand(permuted.bps)[2]) != as.character(strand(permuted.bps)[3])) {
                ## the positive stranded breakend must have lower genomic coordinate
                if (as.character(strand(permuted.bps[2])) == "+" &&
                    start(permuted.bps)[2] < start(permuted.bps)[3] &&
                    start(permuted.bps)[3] - start(permuted.bps)[2] < max_dist) {
                    valid.permutations[ix] = TRUE
                } else if (as.character(strand(permuted.bps[2])) == "-" &&
                           start(permuted.bps)[2] > start(permuted.bps)[3] &&
                           start(permuted.bps)[2] - start(permuted.bps)[3] < max_dist) {
                    valid.permutations[ix] = TRUE
                }
                ## the spanning segment cannot pass both other breakends
                if (intrachromosomal) {
                    if ((pmax(start(permuted.bps[2]), start(permuted.bps[3])) > pmax(start(permuted.bps[4]), start(permuted.bps[1]))) &
                        (pmin(start(permuted.bps[2]), start(permuted.bps[3])) < pmin(start(permuted.bps[4]), start(permuted.bps[1])))) {
                        valid.permutations[ix] = FALSE
                    }
                }
            }
        }
    }

    ## create a graph with just these junctions
    wgg = gG(junctions = junctions)
    snodes.gr = wgg$gr

    make_seg = function(bp, start = TRUE, noextra = FALSE) {
        strand.bp = as.character(strand(bp))
        if (start) {
            ## browser()
            shifted.bp = GenomicRanges::shift(gr.flipstrand(bp),
                               ifelse(as.character(strand(bp)) == "-", -1, 1))
            ## get the gap node
            gap.node = snodes.gr[(snodes.gr %^^% shifted.bp)] %>% reduce
            ## get the node after the gap node
            if (noextra) {
                if (width(gap.node) > pad) {
                    gr = GenomicRanges::resize(gap.node, width = pad, fix = "end")
                } else {
                    gr = gap.node
                }
            } else {
                shifted.gap.node = GenomicRanges::shift(gr.start(gap.node[, c()], ignore.strand = FALSE),
                                                        ifelse(as.character(strand(bp)) == "-", -1, 1))
                next.gap.node = snodes.gr[(snodes.gr %^^% shifted.gap.node)]
            ## resize if its big
                if (length(next.gap.node) && width(next.gap.node) > pad) {
                    next.gap.node = GenomicRanges::resize(next.gap.node, width = pad, fix = "end")
                }
                gr = gUtils::grbind(gap.node, next.gap.node) %>% reduce
            }
        } else {
            ## browser()
            shifted.bp = GenomicRanges::shift(bp, ifelse(as.character(strand(bp)) == "-", -1, 1))
            ## get the gap node
            gap.node = snodes.gr[(snodes.gr %^^% shifted.bp)] %>% reduce
            if (noextra) {
                if (width(gap.node) > pad) {
                    gr = GenomicRanges::resize(gap.node, width = pad, fix = "start")
                } else {
                    gr = gap.node
                }
            } else {
                shifted.gap.node = GenomicRanges::shift(gr.end(gap.node[, c()], ignore.strand = FALSE),
                                                        ifelse(as.character(strand(bp)) == "-", -1, 1))
                next.gap.node = snodes.gr[(snodes.gr %^^% shifted.gap.node)]
                ## resize if its big
                if (length(next.gap.node) && width(next.gap.node) > pad) {
                    next.gap.node = GenomicRanges::resize(next.gap.node, width = pad, fix = "start")
                }
                gr = gUtils::grbind(gap.node, next.gap.node) %>% reduce
            }
            ## ## shift the gap node
            ## shifted.gap.node = GenomicRanges::shift(gap.node[, c()],
            ##                                         ifelse(as.character(strand(bp)) == "-", -1, 1))
            ## gr = snodes.gr[(snodes.gr %^^% shifted.bp) | (snodes.gr %^^% shifted.gap.node)] %>% reduce
            ## shifted.bp = GenomicRanges::resize(shifted.bp, fix = "start", width = pad)
            ## gr = snodes.gr[(snodes.gr %^^% shifted.bp)] %>% reduce
            ## if (strand.bp == "-") {
            ##     gr = GRanges(seqnames = seqnames(bp),
            ##                  ranges = IRanges(start = pmax(1, start(bp) - pad + 1),
            ##                                   end = start(bp)),
            ##                  strand = "-",
            ##                  seqlengths = seqlengths(bp))
            ## } else if (strand.bp == "+") {
            ##     gr = GRanges(seqnames = seqnames(bp),
            ##                  ranges = IRanges(start = end(bp),
            ##                                   end = pmin(end(bp) + pad - 1,
            ##                                              seqlengths(bp)[as.character(seqnames(bp))])),
            ##                  strand = "+",
            ##                  seqlengths = seqlengths(bp))
            ## }
        }
        return(gr)
    }
    
    ## make one-junction walks
    gr1 = make_seg(bps[1, c()], start = TRUE)
    gr2 = make_seg(bps[2, c()], start = FALSE)
    gr3 = make_seg(bps[3, c()], start = TRUE)
    gr4 = make_seg(bps[4, c()], start = FALSE)

    ## get locus annotation
    values(gr1)[, "locus"] = 1
    values(gr2)[, "locus"] = 2
    values(gr3)[, "locus"] = 1
    values(gr4)[, "locus"] = 2

    out.grl = GRangesList(c(gr1, gr2),c(gr3, gr4))
    
    ## make two-junction walks if any are valid
    ## we know what the orientation is based on the number of valid permutations
    if (any(valid.permutations)) {
        which.valid = which(valid.permutations)
        for (ix in which.valid) {
            permutation = breakend.permutations[[ix]]
            permuted.bps = bps[permutation, c()]
            ## the first segment is the first breakend, extended to the end of the chromosome
            if (sum(valid.permutations) > 1) {
                gr1 = make_seg(permuted.bps[1, c()], start = TRUE)
            } else {
                gr1 = make_seg(permuted.bps[1, c()],
                               start = TRUE,
                               noextra = TRUE)
            }
            values(gr1)[, "locus"] = 1
            if (as.character(strand(permuted.bps[2])) == "+") {
                ## the second segment is the middle two breakends
                ## specifically the piece between them
                gr2 = GRanges(seqnames = seqnames(permuted.bps)[2],
                              ranges = IRanges(start = start(permuted.bps[2]),
                                               end = end(permuted.bps[3])),
                              strand = "+",
                              seqlengths = seqlengths(bps))
            } else {
                gr2 = GRanges(seqnames = seqnames(permuted.bps)[2],
                              ranges = IRanges(start = start(permuted.bps[3]),
                                               end = end(permuted.bps[2])),
                              strand = "-",
                              seqlengths = seqlengths(bps))
            }
            values(gr2)[, "locus"] = 2
            ## the third segment is the last breakend, extended to the end of the chromosome
            if (sum(valid.permutations) > 1) {
                gr3 = make_seg(permuted.bps[4, c()], start = FALSE)
            } else {
                gr3 = make_seg(permuted.bps[4, c()],
                               start = FALSE,
                               noextra = TRUE)
            }
            values(gr3)[, "locus"] = 1
            gr = c(gr1, gr2, gr3)
            out.grl[[length(out.grl) + 1]] <- gr
        }
    }

    ## disjoin at breakpoints
    out = gW(grl = out.grl, graph = wgg, disjoin = FALSE)##gr.stripstrand(bps))

    ## set cis or trans annotation
    orientation = unlist(lapply(1:length(out),
                                function(walk.ix) {
                                    edges.dt = out[walk.ix]$edges$dt[, .(alt = class != "REF",
                                                                    edge.ix = 1:.N)]
                                    if (edges.dt[(alt), .N] == 1) {
                                        return("trans")
                                    }
                                    if (edges.dt[(alt), abs(diff(edge.ix))] == 1) {
                                        return("cis")
                                    }
                                    return("trans")
                                })
                         )
    ## annot = ifelse(1:length(out) > 2, "cis", "trans")
    out$set(orientation = orientation)

    return(out)
}

#' @name haplotype_support
#' @title haplotype_support
#'
#' @param wks (gWalk) haplotypes
#' @param reads (GRanges) putative supporting reads
#'
#' @return gWalk with n_supp_bx annotation giving the number of supporting barcodes per walk
haplotype_support = function(wks, reads, min.consec = 3, return.bx = TRUE)
{
    if (is.null(wks$dt$n.alt)) {stop("wks must have n.alt annotation")}
    ## initialized null output
    out = wks$dt[, .(walk.id, n_supp_bx = 0)]
    single.junction.walks = wks[n.alt == 1]
    ## get read support for single junctions
    tumor.res = lapply(1:length(single.junction.walks),
                       function(ix)
                       {
                           res = readsupport:::score.walks(wks = single.junction.walks[ix]$grl,
                                                           reads = reads,
                                                           pad = 0,
                                                           use.discordant = TRUE,
                                                           raw = TRUE)
                           tumor.res = as.data.table(reshape2::melt(as.matrix(res$sc)))[(value > 0), .(BX = Var1, walk.id = ix)]
                           return(tumor.res)
                       })
    tumor.res = rbindlist(tumor.res,fill = TRUE)
    ## get read support for all walks, then subset for "unique barcodes" supporting each walk
    res = readsupport:::score.walks(wks = wks$grl,
                                    reads = reads,
                                    pad = 0,
                                    use.discordant = TRUE,
                                    raw = TRUE)
    all.tumor.res = as.data.table(reshape2::melt(as.matrix(res$sc)))[(value > 0), .(BX = Var1, walk.id = Var2)]
    if (tumor.res[, .N])
    {
        tumor.res[, n.alt := single.junction.walks$dt$n.alt[match(walk.id, single.junction.walks$dt$walk.id)]]
        tumor.res[, edge.id := as.numeric(single.junction.walks$dt$alt.edges[match(walk.id, single.junction.walks$dt$walk.id)])]
        bx.by.junction.res = tumor.res[, .(alt.edges = paste(sort(edge.id), collapse = " ")), by = BX]
        final.walks.edges.dt = lapply(1:length(wks),
                                      function(ix)
                                      {
                                          wks[ix]$edges$dt[, ":="(walk.id = ix)]
                                      }) %>% rbindlist(fill = TRUE)
        final.walks.edges.dt[, n1.positive.gap := wks$graph$nodes$dt[n1, positive.gap]]
        final.walks.edges.dt[, n2.positive.gap := wks$graph$nodes$dt[n2, positive.gap]]
        final.walks.edges.dt[, alt.consecutive := class != "REF" & data.table::shift(class) != "REF"]
        ## get clusters of alt junctions
        final.walks.edges.dt[, cum.alt := cumsum(is.na(data.table::shift(class)) |
                                                 (data.table::shift(class) != class)),
                             by = walk.id]
        cliques = final.walks.edges.dt[(cum.alt > 0), .(alt.edges = .SD[class != "REF", paste(sort(edge.id), collapse = " ")]), by  = .(cum.alt, walk.id)]
        nbx = lapply(cliques$alt.edges, function(eds) {bx.by.junction.res[alt.edges == eds, .N]}) %>% unlist
        cliques[, n_supp_bx := nbx]
        out = cliques[, .(n_supp_bx = sum(n_supp_bx)), by = walk.id]
    }
    out[, n_supp_ubx := 0]
    if (all.tumor.res[, .N])
    {
        all.tumor.res[, nwalks := .SD[, length(unique(walk.id))], by = BX]
        if (all.tumor.res[nwalks == 1,.N])
        {
            nwalks.dt = all.tumor.res[(nwalks == 1), .(n_supp_ubx = .SD[, length(unique(BX))]), by = walk.id]
            out[, n_supp_ubx := nwalks.dt[, n_supp_ubx][match(walk.id, nwalks.dt$walk.id)]]
            out[is.na(n_supp_ubx), n_supp_ubx := 0]
        }
    }
    return(out)
}


#' @name cluster2walks
#' @title cluster2walks
#'
#' @description
#' get all walks associated with an edge cluster
#' 
#' @param junctions (Junction) involved rearrangements
#' @param gaps (GRanges) signed gaps
#' @param verbose (logical) default FALSE
cluster2walks = function(junctions, gaps, verbose = FALSE)
{
    ## helper function for embedding circular walks into linear ones
    ## returns a GRangesList corresponding to the embedded walk
    ## params:
    ## circular.grl: length one GRangesList for circular walk
    ## linear.grl: length one GRangesList for linear walk
    ## circular.ix: which segment of the circular walk intersects the linear one?
    ## linear.ix: which segment of the linear walk intersects the circular one?
    ## flip.circular: flip strand of the circular walk
    circular2linear = function(circular.grl, linear.grl, circular.ix, linear.ix, flip.circular)
    {
        circular.gr = circular.grl[[1]]
        linear.gr = linear.grl[[1]]
        if (flip.circular) {
            circular.gr = rev(gr.flipstrand(circular.gr))
            circular.ix = length(circular.gr) - circular.ix + 1
        }
        gr = GRanges(seqlengths = seqlengths(linear.gr))
        ## add linear walk up to and not including the intersectin node
        if (linear.ix > 1) {
            gr = c(gr, linear.gr[1:linear.ix-1])
        }
        ## add circular walk from intersecting node to the end
        if (circular.ix <= length(circular.gr)) {
            gr = c(gr, circular.gr[circular.ix:length(circular.gr)])
        }
        ## add circular walk from start up to and not including the intersecting node
        if (circular.ix > 1) {
            gr = c(gr, circular.gr[1:circular.ix - 1])
        }
        ## add linear walk from intersecting node to the end
        if (linear.ix <= length(linear.gr)) {
            gr = c(gr, linear.gr[linear.ix:length(linear.gr)])
        }
        out = GRangesList(gr, compress = FALSE)
        return(out)
    }
    ## create a new gGraph
    new.gg = gG(junctions = junctions)
    ## mark nodes to reflect whether they are part of some (positive) gap
    ## filter out width one ranges
    gaps = gaps %Q% (width(gaps) > 1)
    new.gg$nodes$mark(gap = (new.gg$nodes$gr %O% gaps) > 0.99)
    new.gg$nodes$mark(positive.gap = (new.gg$nodes$gr %OO% gaps) > 0.99)
    ## get preliminary walks
    prelim.walks = new.gg$walks()
    ## get segments involved in each walk
    all.walks.grl.dt = as.data.table(stack(prelim.walks$grl))
    all.walks.grl.dt[, grl.iix := 1:.N, by = walk.id]
    ## get edges involved in each walk
    all.walks.edges.dt = lapply(1:length(prelim.walks),
                                function(ix)
                                {
                                    prelim.walks[ix]$edges$dt[, ":="(walk.id = ix,
                                                                     circular = prelim.walks$dt[ix, circular])]
                                }) %>% rbindlist(fill = TRUE)
    ## check to see if there are circular walks
    if (any(prelim.walks$dt[, circular]))
    {
        ## cartesian join of circular and linear walks by UNSIGNED node id
        merged.walks.dt = merge.data.table(all.walks.grl.dt[(circular), .(walk.id, snode.id, node.id, grl.iix)],
                                           all.walks.grl.dt[(!circular), .(walk.id, snode.id, node.id, grl.iix)],
                                           by = "node.id",
                                           suffixes = c(".circular", ".linear"),
                                           allow.cartesian = TRUE)
        ## we need to remove pairs of walks which have overlapping junctions
        ## this is so that per walk each junction appears at most once
        merged.walks.edges.dt = merge.data.table(all.walks.edges.dt[(circular) & (class != "REF"),
                                                                    .(walk.id, edge.id)],
                                                 all.walks.edges.dt[(!circular) & (class != "REF"),
                                                                    .(walk.id, edge.id)],
                                                 by = "edge.id",
                                                 suffixes = c(".circular", ".linear"),
                                                 allow.cartesian = TRUE)
        merged.walks.dt[, pair.ix := paste(walk.id.circular, walk.id.linear)]
        merged.walks.edges.dt[, pair.ix := paste(walk.id.circular, walk.id.linear)]
        embed.dt = merged.walks.dt[!(pair.ix %in% merged.walks.edges.dt$pair.ix),]
        ## extract circular and linear walks
        circular.walks = prelim.walks[embed.dt[, walk.id.circular]]
        linear.walks = prelim.walks[embed.dt[, walk.id.linear]]
        ## get GRangeList of embedded walks
        grl.lst = lapply(1:embed.dt[, .N],
                         function(ix)
                         {
                             return(circular2linear(circular.grl = circular.walks[ix]$grl,
                                                    linear.grl = linear.walks[ix]$grl,
                                                    circular.ix = embed.dt[ix, grl.iix.circular],
                                                    linear.ix = embed.dt[ix, grl.iix.linear],
                                                    flip.circular = embed.dt[ix, snode.id.circular == -snode.id.linear]))
                         })
        embed.grl = do.call("c", grl.lst)
    } else {
        embed.grl = GRangesList()
    }
    ## get the longest linear walk
    all.walks.edges.dt[(!circular), nedges := .N, by = walk.id]
    all.walks.edges.dt[(!circular), n.alt := .SD[class != "REF", .N], by = walk.id] 
    longest.linear.walks = all.walks.edges.dt[(!circular),
                                              .(longest.walk.id = .SD[which.max(nedges)[1], walk.id],
                                                class = .SD[1, class]),
                                              by = .(edge.id, n.alt)][class != "REF"]
    ## select those walks plus
    linear.walk.ids = longest.linear.walks[class != "REF", unique(longest.walk.id)]
    linear.walks.grl = prelim.walks[linear.walk.ids]$grl
    ## prelim.walks$graph$edges$dt[class != "REF", .(bp1, bp2, class)]
    if (length(embed.grl))
    {
        final.walks = gW(grl = c(linear.walks.grl, embed.grl), graph = prelim.walks$graph)
    }
    else
    {
        final.walks = gW(grl = linear.walks.grl, graph = prelim.walks$graph)
    }
    ## remove non-unique walks
    uwalks = which(!duplicated(grl.string(final.walks$grl)))
    final.walks = final.walks[uwalks]
    ## add gap status, cis/trans, edge.id, number of ALT  annotations
    final.walks.edges.dt = lapply(1:length(final.walks),
                                  function(ix)
                                  {
                                      final.walks[ix]$edges$dt[, ":="(walk.id = ix)]
                                  }) %>% rbindlist(fill = TRUE)
    final.walks.edges.dt[, n1.positive.gap := final.walks$graph$nodes$dt[n1, positive.gap]]
    final.walks.edges.dt[, n2.positive.gap := final.walks$graph$nodes$dt[n2, positive.gap]]
    final.walks.edges.dt[, alt.consecutive := class != "REF" & data.table::shift(class) != "REF"]
    walk.orientation.dt = final.walks.edges.dt[,.(cis = any(alt.consecutive & n1.positive.gap, na.rm = TRUE),
                                                  n.cis = sum(alt.consecutive & n1.positive.gap, na.rm = TRUE)),
                                               by = walk.id]
    walk.alt.edges.dt = final.walks.edges.dt[, .(n.alt = .SD[class != "REF", .N],
                                                 alt.edges = .SD[class != "REF", paste(sort(edge.id),
                                                                                       collapse = ' ')],
                                                 keep = .SD[class != "REF", length(unique(edge.id)) == length(edge.id)]),
                                             by = walk.id]
    final.walks$set(cis = walk.orientation.dt[match(final.walks$dt$walk.id, walk.id), cis])
    final.walks$set(n.cis = walk.orientation.dt[match(final.walks$dt$walk.id, walk.id), n.cis])
    final.walks$set(n.alt = walk.alt.edges.dt[match(final.walks$dt$walk.id, walk.id), n.alt])
    final.walks$set(alt.edges = walk.alt.edges.dt[match(final.walks$dt$walk.id, walk.id), alt.edges])
    final.walks$set(keep = walk.alt.edges.dt[match(final.walks$dt$walk.id, walk.id), keep])
    ## formatting
    final.walks$set(name = final.walks$dt$walk.id)
    final.walks$set(col = ifelse(final.walks$dt$cis, "red", "blue"))
    return(final.walks[keep == TRUE])
}

