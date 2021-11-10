#' @import data.table
#' @import GenomicRanges
#' @import ggplot2
#' @import gTrack
#' @import gGnome
#' @import gUtils
#' @importFrom igraph graph.adjacency shortest_paths
#' @importFrom rtracklayer import import.bw
#' @importFrom skitools rel2abs ppng

#' @name ecycles
#' @title ecycles
#'
#' @description
#'
#' Return a data table of all ALT edges in a gGraph that is part of a simple cycle
#'
#' @param junctions (character) path to junctions file or GRangesList
#' @param jabba_rds (character) path to JaBbA output
#' @param thresh (numeric) distance threshold for reciprocality, default 1e3 (bp)
#' @param min.span (numeric) thres for removing small dups and dels (1e3 bp)
#' @param chunksize (numeric) default 1e30
#' @param mc.cores (numeric) default 1
#' @param verbose (logical) default FALSE
#'
#' @return data.table with columns cycle.id, edge.id, class, junction
#' @export
ecycles = function(junctions = NULL, jabba_rds = NULL, thresh = 1e3,
                   min.span = 1e3,
                   chunksize = 1e30, mc.cores = 1, verbose = FALSE) {

    ## empty output
    res = data.table(cycle.id = numeric(),
                     edge.id = numeric(),
                     class = character(),
                     junction = character())

    ## check if gGraph should be created from junctions
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


    ## remove small dups and dels
    espans = gg$junctions[type == "ALT"]$span
    eclass = gg$junctions[type == "ALT"]$dt$class
    ggjuncs = gg$junctions[type == "ALT"][espans > min.span | eclass == "INV-like" | eclass == "TRA-like"]
    if (verbose){
        message("Number of junctions after removing small dups and dels: ", length(ggjuncs))
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
    
    bp = grl.unlist(altedges$grl)[, c("grl.ix", "grl.iix")]
    bp.dt = gr2dt(bp)

    ix = split(1:length(bp), ceiling(runif(length(bp))*ceiling(length(bp)/chunksize)))
    ixu = unlist(ix)
    eps = 1e-9
    ij = do.call(rbind, split(1:length(bp), bp$grl.ix))
    xt.adj = old.adj = Matrix::sparseMatrix(1, 1, x = 0, dims = rep(length(bp), 2))

    if (verbose){
        message(sprintf('Computing junction graph across %s ALT edges with distance threshold %s', length(altedges), thresh))
    }

    if (!exists(".INF")){
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

    ## do single linkage hierarchical clustering within `range`
    hcl = stats::hclust(as.dist(xt.adj), method = "single")
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
    for (rw in unique.fr) {
        fr = dt[row == rw, col]
        cy = igraph::shortest_paths(gr, from = fr, to = rw)$vpath
        sel = which(sapply(cy, length) > 1)
        for (ix in sel) {
            num.cy = num.cy + 1
            all.cy[[num.cy]] = sort(as.numeric(cy[sel][[ix]]))
        }
    }

    if (!length(all.cy)) {
        if (verbose) {
            message("No cycles found!")
        }
        return(res)
    } 
    unique.cy = unique(all.cy)
    jcl = lapply(unique.cy, function(x) unique(sort(bp$grl.ix[x]))) %>% unique
    dcl = dunlist(unname(jcl)) %>% setnames(new = c('listid', 'edges'))
    dcl[, class := altedges$dt[dcl$edges, class]]
    dcl[, junction := grl.string(altedges$grl[dcl$edges])]
    return(dcl)
}


#' @name grab_cov_gtrack
#' @title grab_cov_gtrack
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
                           ...) {
    if (grepl("txt", cov)) {
        cov.gr = dt2gr(fread(cov))
    } else if (grepl("rds", cov)) {
        cov.gr = readRDS(cov)
    } else {
        stop("invalid file name")
    }

    if (!field %in% names(values(cov.gr))) {
        stop("invalid field supplied")
    }

    if (!is.null(jab) && file.exists(jab)) {
        purity = readRDS(jab)$purity
        ploidy = readRDS(jab)$ploidy
    }

    if (!is.null(purity) && !is.null(ploidy)) {
        cov.gr$cn = skitools::rel2abs(cov.gr, field = field, purity = purity, ploidy = ploidy, allele = FALSE)
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
        stop("ascat file not valid")
    }
    if (!file.exists(sequenza) | !file.info(ascat)$size) {
        stop("sequenza file not valid")
    }
    if (!(breaks %in% c('ascat', 'sequenza'))) {
        warning('breaks must be one of ascat, sequenza. using ascat by default')
        breaks = 'ascat'
    }
    ascat.segs = grab_segs(bw = ascat,
                           field = 'score',
                           simplify = TRUE,
                           verbose = verbose)
    sequenza.segs = grab_segs(csv = sequenza,
                              field = 'CNt',
                              simplify = TRUE,
                              verbose = verbose)

    ## get loose ends from each segmentation alg
    ascat.loose = grab_loose_from_seg(ascat.segs,
                                      min.width = min.width,
                                      max.gap = max.gap,
                                      return.type = 'GRanges')
    sequenza.loose = grab_loose_from_seg(sequenza.segs,
                                         min.width = min.width,
                                         max.gap = max.gap,
                                         return.type = 'GRanges')

    ## filter for overlapping
    if (breaks == 'ascat') {
        keep.le = ascat.loose %^^% (sequenza.loose + pad)
        le = ascat.loose[keep.le]
    } else {
        keep.le = sequenza.loose %^^% (ascat.loose + pad)
        le = sequenza.loose[keep.le]
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
#' 
#' @return loose ends of a JaBbA graph as GRanges
#' if gGraph is not balanced, will produce a warning and return empty GRanges
gg_loose_end = function(jabba_rds = NA_character_, return.type = "GRanges") {
    
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
#' @param coverage (character) path to coverage file with gc-adjusted counts
#' @param tumor.field (character) coverage field (ratio)
#' @param norm.field (character) normal count field (norm.counts)
#' @param mask (character) path to coverage mask
#' @param mask_pad (numeric) pad on distance to mask
#' @param bad_thresh (numeric) default NA (skip this)
#' @param norm_thresh (numeric) 0.6 (from filter.loose)
#' @param corr_thresh (numeric) 0.5 autocorrelation threshold
#' @param min_size (numeric) minimum size in BP
#' @param subset (logical) subset loose ends to include only valid ones before returning
#' @param verbose (logical)
grab_loose = function(jabba_rds,
                      coverage = "/dev/null",
                      tumor_coverage = coverage,
                      norm.field = 'norm.counts',
                      tumor.field = 'foreground',
                      cov.pad = 1e5,
                      mask = "~/projects/gGnome/files/new.mask.rds",
                      mask_pad = 1e4,
                      bad_thresh = NA,##0.4,
                      norm_thresh = 0.6,
                      corr_thresh = 1, ## autocorrelation threshold
                      min_size = 5e4,
                      subset = FALSE,
                      verbose = TRUE) {

    if (!check_file(jabba_rds)) {
        stop("JaBbA file path not valid")
    }

    jabba_rds = file.path(dirname(jabba_rds), "jabba.raw.rds")
    jabba_simple_rds = file.path(dirname(jabba_rds), "jabba.simple.rds")

    if (!check_file(coverage)) {
        if (verbose) {
            message("Coverage file not supplied!")
        }
        use.coverage = FALSE
    } else {
        if (verbose) {
            message("Reading normal coverage")
        }
        cov.gr = readRDS(coverage)
        use.coverage = TRUE
    }

    if (!check_file(tumor_coverage)) {
        if (verbose) {
            message("Tumor coverage file not valid")
        }
        use.tumor.coverage = FALSE
    } else {
        if (verbose) {
            message("Reading tumor coverage.")
        }
        use.tumor.coverage = TRUE
        tum.cov.gr = readRDS(tumor_coverage)
    }

    if (use.coverage) {
        if (!(norm.field %in% names(values(cov.gr)))) {
            stop("supplied field not in coverage metadata: ", norm.field)
        }
        cov.gr = cov.gr %Q% (!is.infinite(values(cov.gr)[[norm.field]]))
    }

    if (use.tumor.coverage) {
        if (!(tumor.field %in% names(values(tum.cov.gr)))) {
            stop("supplied field not in tumor_coverage metadata: ", tumor.field)
        }
        tum.cov.gr = tum.cov.gr %Q% (!is.infinite(values(tum.cov.gr)[[tumor.field]]))
    }

    ## pull loose ends from jabba_rds
    if (verbose) {
        message("Grabbing loose ends from JaBbA input")
    }
    le.gr = gg_loose_end(jabba_rds, return.type = 'GRanges')
    le.dt = as.data.table(le.gr)

    if (!nrow(le.dt)) {
        if (verbose) {
            message("No loose ends!")
        }
        return(le.dt)
    }
    
    le.gr$loose.index = 1:length(le.gr)

    ## check if loose end overlaps mask
    if (file.exists(mask)) {
        mask.gr = readRDS(mask)
        overlap.mask = le.gr %^% (gr.stripstrand(mask.gr) + mask_pad)
        le.dt[, mask := overlap.mask]
    }

    ## get node ids of fused and unfused loose ends
    if (verbose) {
        message("Identifying fused and unfused sides of each loose end")
    }
    
    fused.unfused.dt = fused_unfused(le.dt, jabba_rds)
    simple.fused.unfused.dt = fused_unfused(le.dt, jabba_simple_rds)

    ## grab gGraph
    gg = gG(jabba = jabba_rds)
    gg.simple = gG(jabba = jabba_simple_rds)

    ## get flanking normal coverage (borrowed from filter.loose)
    fused.melted.dt = melt.data.table(fused.unfused.dt, id.vars = "loose.index",
                                   measure.vars = c("fused.node.id", "unfused.node.id"),
                                   variable.name = "fused",
                                   value.name = "node.id")[!is.na(node.id),]
    
    segs.gr = gg$nodes$gr[fused.melted.dt$node.id, c()]
    segs.gr$loose.index = fused.melted.dt$loose.index
    segs.gr$fused = fused.melted.dt$fused == "fused.node.id"

    
    ## extend fused and unfused size 100 mbp in each direction, then overlap with coverage
    sides.gr = gr.findoverlaps(segs.gr, le.gr + cov.pad, by = 'loose.index')
    values(sides.gr) = cbind(values(sides.gr), values(le.gr[sides.gr$subject.id]))
    sides.gr$fused = segs.gr$fused[sides.gr$query.id]

    if (use.coverage) {

        ## overlap with normal coverage
        sides.cov.dt = gr.findoverlaps(cov.gr, sides.gr, return.type = 'data.table')
        sides.cov.dt[, normal.cov := values(cov.gr)[[norm.field]][query.id]]
        sides.cov.dt[, fused := values(sides.gr)$fused[subject.id]]
        sides.cov.dt[, loose.index := values(sides.gr)$loose.index[subject.id]]

    }

    if (use.tumor.coverage) {

        ## get tumor coverage to compute autocorrelation
        tumor.cov.dt = gr.findoverlaps(tum.cov.gr, sides.gr, return.type = 'data.table')
        tumor.cov.dt[, tumor.cov := values(tum.cov.gr)[[tumor.field]][query.id]]
        tumor.cov.dt[, fused := values(sides.gr)$fused[subject.id]]
        tumor.cov.dt[, loose.index := values(sides.gr)$loose.index[subject.id]]

        ## compute waviness at loose ends
        if(verbose) {
            message("Calculating waviness around loose end")
        }
        ## autocorrelation.dt = tumor.cov.dt[, .(waviness = max(.waviness(start[fused], tumor.cov[fused]),
        ##                                                      .waviness(start[!fused], tumor.cov[!fused]),
        ##                                                      na.rm=T)), by=loose.index]

        tumor.cov.dt[, count := .N, by = .(loose.index, fused)]
        autocorrelation.dt = tumor.cov.dt[count > 1,
                                          .(autocorr = max(as.numeric(acf(tumor.cov)$acf)[-1],
                                                           na.rm = TRUE)),
                                          by = .(fused, loose.index)][, .(autcorr = max(autocorr, na.rm = T)),
                                                                      by = loose.index]

        ## add autocorrelation info
        le.dt[, autocorrelation := autocorrelation.dt$autcorr[match(loose.index, autocorrelation.dt$loose.index)]]
        le.dt[, highcorr := autocorrelation >= corr_thresh]

    } else {

        if (verbose) {
            message("Skipping autocorrelation calculation as no coverage supplied")
        }
        le.dt[, autocorrelation := NA]
        le.dt[, highcorr := FALSE]
    }

    if (use.coverage) {
        ## compute mean coverage for fused and unfused sides
        mean.cov.dt = sides.cov.dt[, .(mean.normal.cov = mean(normal.cov, na.rm = T)),
                                   by = .(loose.index, fused)]
        mean.cov.dt[, fused := ifelse(fused, 'fused', 'unfused')]
        mean.cov.dt = dcast.data.table(mean.cov.dt, loose.index ~ fused, value.var = "mean.normal.cov")

        ## determine whether it's higher than normal beta but only consider only diploid chr
        autosomal.chrs = grep("(^(chr)*[0-9]+$)", names(seqlengths(cov.gr)), value = TRUE)
        norm.beta = mean(as.data.table(cov.gr)[seqnames %in% autosomal.chrs,][[norm.field]], na.rm = TRUE) * 0.5 ## factor of 0.5 is for allelic CN
        mean.cov.dt[, norm.change := norm_thresh * norm.beta < abs(fused - unfused)]

        ## mark normal change and copy normal fused and unfused coverage
        le.dt$norm.change = le.dt$loose.index %in% mean.cov.dt[(norm.change), loose.index]
        le.dt$norm.fused = mean.cov.dt$fused[match(le.dt$loose.index, mean.cov.dt$loose.index)]
        le.dt$norm.unfused = mean.cov.dt$unfused[match(le.dt$loose.index, mean.cov.dt$loose.index)]
    } else {
        le.dt$norm.change = FALSE
        le.dt$norm.fused = NA
        le.dt$norm.unfused = NA
    }
    
    ## get size of loose end
    simple.fused.unfused.dt[, fused.size := gg.simple$nodes$dt$width[fused.node.id]]
    simple.fused.unfused.dt[, unfused.size := gg.simple$nodes$dt$width[unfused.node.id]]
    simple.fused.unfused.dt[, size := pmin(fused.size, unfused.size, na.rm = TRUE)]
    simple.fused.unfused.dt[, small := size <= min_size]

    ## ## get CN of each side of loose end
    fused.unfused.dt[, fused.cn := gg$nodes$dt$cn[fused.node.id]]
    fused.unfused.dt[, unfused.cn := gg$nodes$dt$cn[unfused.node.id]]
    fused.unfused.dt[, fused.lower := (fused.cn <= unfused.cn)]

    ## check whether this sample overlaps a bad part of the coverage
    ## fused.unfused.dt[, fused.cn.old := gg$nodes$dt$cn.old[fused.node.id]]
    ## fused.unfused.dt[, unfused.cn.old := gg$nodes$dt$cn.old[unfused.node.id]]

    ## merge these with loose ends
     le.dt = merge.data.table(le.dt, fused.unfused.dt, by = "loose.index", all.x = TRUE)
    le.dt = merge.data.table(le.dt,
                             simple.fused.unfused.dt[, .(loose.index, fused.size, unfused.size, size, small)],
                             by = "loose.index", all.x = TRUE)

    ## ## indicate whether to "keep"
    ## le.dt[, coverage.bad := FALSE]
    ## if (!is.na(bad_thresh)) {
    ##     le.dt[!is.na(fused.cn.old) & fused.size > unfused.size,
    ##           coverage.bad := fused.cn - fused.cn.old > bad_thresh]
    ##     le.dt[!is.na(unfused.cn.old) & fused.size < unfused.size,
    ##           coverage.bad := unfused.cn.old - unfused.cn > bad_thresh]
    ## } else {
    ##     le.dt[, coverage.bad := FALSE]
    ## }
    
    le.dt[, keep := (!(mask) & !(norm.change) & !(fused.lower) & !(small) & !(highcorr))]

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
                              gs = NULL,
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
#' @param mask (character) masked ranges (ignore breakpoints in this region)
#' @param verbose (logical) default FALSE
#'
#' @return list with entries:
#' - $cncp.gr: GRanges of CN change points in the segmentation + distance to nearest gold-standard point
#' - $cncp.res: containing F1 score, precision, and recall of CN change points
cncp_analysis = function(seg = NULL,
                         gs = NULL,
                         method = "jabba",
                         id = "sample",
                         max.dist = 1e4,
                         mask = "~/projects/gGnome/files/new.mask.rds",
                         simplify = FALSE,
                         verbose = FALSE) {

    if (verbose) {
        message("reading gold standard segments")
    }

    gs.cncp = stack(readRDS(gs))

    if (grepl(pattern = "ascat", x = method, ignore.case = TRUE)) {
        method.gr = grab_segs(bw = seg, field = "score",
                              simplify = simplify,
                              verbose = verbose)
    } else if (grepl(pattern = "jabba", x = method, ignore.case = TRUE)) {
        method.gr = grab_segs(gg = seg,
                              simplify = simplify,
                              verbose = verbose)
    } else if (grepl(pattern = "sequenza", x = method, ignore.case = TRUE)) {
        method.gr = grab_segs(csv = seg, field = "CNt",
                              simplify = simplify,
                              verbose = verbose)
    } else {
        stop("method not implemented yet")
    }

    if (file.exists(mask) && file.info(mask)$size) {
        mask.gr = readRDS(mask)
    } else {
        mask.gr = GRanges(seqlengths = seqlengths(method.gr))
    }

    method.cncp = grab_cncp(gr = method.gr, verbose = verbose, include.na = TRUE)

    ## remove masked breakpoints
    method.cncp = method.cncp %Q% (!(method.cncp %^% mask.gr))
    gs.cncp = gs.cncp %Q% (!(gs.cncp %^% mask.gr))
    
    out = compare_cncp(gr = method.cncp,
                       gs = gs.cncp,
                       max.dist = max.dist,
                       verbose = verbose)

    ## cncp.res = bp.analysis(method.cncp, gs.cncp, pad = max.dist)
    out.gr = out$cncp.gr
    out.gr$pair = id

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
                            min.width = 1e3,
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
#' @param gr (GRanges) allow GRanges input directly
#' @param gs (character) path to gold standard segs
#' @param method (character) analysis method (one of jabba, ascat)
#' @param tile.width (numeric) resolution at which to tile the genome for correlation. Default 1e4.
#' @param min.width (numeric) min tile width to consider for RMSE, etc.
#' @param id (character) sample ID
#' @param mask (character) path to coverage mask
#' @param verbose (logical) default FALSE
#'
#' @return list with entries $cncp.gr $cncp.res
cn_analysis = function(seg = NULL,
                       gr = NULL,
                       gs = NULL,
                       method = "jabba",
                       id = "sample",
                       tile.width = 1e4,
                       min.width = 1e3,
                       mask = "~/projects/gGnome/files/new.mask.rds",
                       verbose = FALSE) {

    if (verbose) {
        message("reading gold standard segments")
    }

    gs.gr = grab_segs(gr = gs, field = "cn", simplify = TRUE, verbose = verbose)

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

    ## remove centromeric regions
    method.gr = remove_centromeres(method.gr)

    ## remove tiles overlapping with coverage mask if provided
    if (file.exists(mask) && file.info(mask)$size) {
        mask.gr = readRDS(mask)
    } else {
        mask.gr = GRanges(seqlengths = seqlengths(method.gr))
    }
    
    tiles = gr.tile(gs.gr, width = tile.width)
    
    cn.comp = compare_cn(gr = method.gr, tiles = tiles, nm = method, verbose = verbose)
    cn.comp = compare_cn(gr = gs.gr, tiles = cn.comp, nm = "gs", verbose = verbose)

    cn.comp = cn.comp %Q% (!(cn.comp %^% mask.gr))
    
    cn.cor = compute_cn_correlation(gr = cn.comp,
                                    method = method,
                                    gs = "gs",
                                    min.width = min.width,
                                    verbose = verbose)

    cn.comp$pair = id
    cn.cor[, pair := id]
    return(list(cn.res = cn.cor, cn.seg = cn.comp))
}



#' @name check_file
#' @title check_file
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
#' @param simplify (logical) default TRUE. Merge adjacent ranges with the same copy number?
#' @param verbose (logical)
#' 
#' @return GRanges with metadata field $cn
grab_segs = function(gg = "/dev/null",
                     gr = "/dev/null",
                     bw = "/dev/null",
                     csv = "/dev/null",
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
        setnames(dt, colnames(dt)[1:3], c("seqnames", "start", "end"))
        this.csv = dt2gr(dt)
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

    ## fill in gaps for segments taht are not perfectly adjacent
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

    new.segs = simplify_segs(gr = segs, sort.gr = TRUE, verbose = verbose)

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
    
    new.segs = sortSeqlevels(new.segs)
    new.segs = sort(new.segs)

    return(new.segs)
}

#' @name grab_cncp
#' @title grab_cncp
#' 
#' @description
#'
#' get copy number changepoints from simplifed CN segmentation
#'
#' @param gr (GRanges) with metadata field $cn
#' @param include.na (logical) include segments with NA CN as CN changepoint
#' @param verbose (logical)
#'
#' @return GRanges each with width 1 giving CN change points
grab_cncp = function(gr = NULL,
                     include.na = TRUE,
                     verbose = FALSE) {

    if (is.null(gr$cn)) {
        stop("gr missing $cn metadata")
    }

    if (verbose) {
        message("Sorting gr...")
    }

    gr = sortSeqlevels(gr)
    gr = sort(gr)

    dt = as.data.table(gr)

    dt[, shifted.cn := data.table::shift(cn)]
    dt[, shifted.seqnames := data.table::shift(seqnames)]

    if (include.na) {
        cncp.dt = dt[(shifted.cn != cn | is.na(shifted.cn) | is.na(cn)) & shifted.seqnames == seqnames,]
    } else {
        cncp.dt = dt[shifted.cn != cn & !is.na(cn) & !is.na(shifted.cn) & shifted.seqnames == seqnames,]
    }
    
    cncp.gr = GRanges(seqnames = cncp.dt[, seqnames],
                      ranges = IRanges(start = cncp.dt[, start], width = 1),
                      cn = cncp.dt[, cn])

    return(cncp.gr)
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
    tiles = gr.fix(tiles, sl)
    
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
compare_cncp = function(gr = NULL, gs = NULL, max.dist = 1e3, verbose = FALSE) {

    if (verbose) {
        message("Computing pairwise distances between CN change points")
    }

    ov = gr.dist(gr1 = gr, gr2 = gs, ignore.strand = TRUE)
    ov.dt = data.table(which(ov<max.dist, arr.ind=T))
    setnames(ov.dt, c("row", "col"), c("gr.id", "gs.id"))
    ov.dt[, dist := ov[cbind(gr.id, ov.dt$gs.id)]]

    if (verbose) {
        message("computing match stats")
    }

    res = data.table(n.cncp = length(gr),
                     n.cncp.gs = length(gs),
                     tp = length(unique(ov.dt[, gr.id])))

    res[, fp := n.cncp - tp]
    res[, fn := n.cncp.gs - length(unique(ov.dt[, gs.id]))]
    res[, precision.cncp := tp / (tp + fp)]
    res[, recall.cncp := tp / (tp + fn)]
    res[, f1.cncp := 2 * precision.cncp * recall.cncp / (precision.cncp + recall.cncp)]

    ## prepare query GRanges and distance to closest gold standard GRanges
    cncp.gr = gr[, c()]
    closest.dt = ov.dt[, .(gs.id = .SD$gs.id[which.min(.SD$dist)],
                           dist = min(.SD$dist, na.rm = TRUE)),
                       by = gr.id]
    cncp.gr$closest.gs = closest.dt$gs.id[match(1:length(cncp.gr), closest.dt$gr.id)]
    cncp.gr$dist = closest.dt$dist[match(1:length(cncp.gr), closest.dt$gr.id)]

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

    gr = gr %Q% (width(gr) >= min.width)

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

