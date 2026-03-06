## Simulation functions

#Load to other scripts with:
#library(kinship2) 
#library(optiSel)
#source("SimFunctions.R")

## Parameters
# - Simulate offspring
#   - max # of progeny per cross [set max.progeny and then it will sim 1:max]
# - Simulate choosing pair crosses
#   - # pair crosses per generation
#   - Random mating, but with:
#     - maximum kinship allowed (e.g., 0.625 to avoid mating cousins)
#     - do not pair up multiple fish from the same 2 families [to avoid double cousins 
#       which are genetically the same as half-sibs]
#     - do not cross more than 4 fish from any family 
#     - each fish can only mate once
#     - fish availability scheme
#       - # of spawning sessions (not all fish available from beginning of season)
#       - proportion of fish to make available first week
#       - proportion of fish lost after each session (sim. fish dying & females never ripe again)
#     - # of tries to randomly pair fish per session, given other parameters (because it is 
#       random, it may run out of useable pairs even if a different combo of fish would have
#       been possible, so let it try assigning pairs multiple times)

## Conventions
# - New genIDs are coded as gen_fam_ind (e.g., first progeny for sim offspring
#   after latest pedigree PCs (e.g., 190011x190012) will be 20_001_01, 20_001_02, etc.)


##### Functions ####

simulate_progeny <- function(ped.sim, pairs, max.progeny=10){
  
  # Create progeny df
  progeny <- NULL
  # Loop through PCs to generate progeny with unique IDs
  for (i in 1:nrow(pairs)) {
    #Randomize # of progeny [1-n.progeny] 
    #LATER: change to 0-n.progeny, then add if statement for rest of function to bypass if 0
    n.progeny <- sample(1:max.progeny,1)
    
    #Identify generation for unique IDs
    if (!grepl("_", pairs$Sire[i])){
      g <- as.integer(substr(pairs$Sire[i],1,2)) + 1
    } else {
      g <- as.integer(sub("_.*$","", pairs$Sire[i])) + 1
    }
    
    #ID can have any gen #, PC # up to 3 digits, progeny # up to 3 digits
    x <- data.frame(ID=paste(g, formatC(pairs$PC[i],width=3,flag="0"), 
                             #formatC(1:10,width=3,flag="0"),sep="_"),
                             formatC(1:n.progeny,width=3,flag="0"),sep="_"),
                    Sire=pairs$Sire[i],Dam=pairs$Dam[i],
                    PC=pairs$PC[i],
                    yr=pairs$yr[i]+1)
    progeny <- rbind(progeny,x)
  }
  
  # Combine with original pedigree
  ped.sim <- rbind(ped.sim,progeny)
  return(ped.sim)
}
#Run function
#ped.sim <- simulate_progeny(ped.sim,pairs,max.progeny=10)


choose_pairs_avail <- function(spawners, # list of available spawners
                               Ks, # kinship matrix
                               target_pairs=312,
                               kin_thresh=0.0625, # avoid cousins
                               max.fam=4, # number of spawners per family
                               weeks=2, # number of spawning weeks
                               loss=0, # proportion of fish lost each week
                               early_frac=NULL, # fraction available in week 1; NULL = even
                               tries=1) {
  
  ids_all <- as.character(spawners$ID)
  fam <- setNames(as.character(spawners$PC), ids_all)
  n_spawners <- length(ids_all)
  
  ## crosses per week
  crosses_perweek <- floor(target_pairs / weeks)
  crosses_lastweek <- target_pairs - crosses_perweek * (weeks - 1)
  
  ## divide spawners among weeks
  spawner_groups <- vector("list", weeks)
  
  if(is.null(early_frac)) {
    
    ## EVEN split across weeks
    spawners_perweek <- floor(n_spawners / weeks)
    idx <- 1
    for(w in 1:(weeks - 1)) {
      spawner_groups[[w]] <- ids_all[idx:(idx + spawners_perweek - 1)]
      idx <- idx + spawners_perweek
    }
    spawner_groups[[weeks]] <- ids_all[idx:n_spawners]
    
  } else {
    
    ## FRONT-LOADED split
    if(early_frac <= 0 || early_frac >= 1)
      stop("early_frac must be between 0 and 1")
    
    n_early <- ceiling(n_spawners * early_frac)
    spawner_groups[[1]] <- ids_all[1:n_early]
    
    remaining <- ids_all[(n_early + 1):n_spawners]
    
    #Linear decay split among remaining weeks
    n_weeks_remaining <- weeks - 1
    n_remaining <- length(remaining)
    
    if(n_weeks_remaining > 0 && n_remaining > 0) {
      # linear decay: assign decreasing counts each week
      weights <- rev(seq_len(n_weeks_remaining))
      counts <- floor(weights / sum(weights) * n_remaining)
      
      # adjust last week for rounding
      counts[n_weeks_remaining] <- n_remaining - sum(counts[-n_weeks_remaining])
      
      idx <- 1
      for(w in 2:weeks) {
        spawner_groups[[w]] <- remaining[idx:(idx + counts[w-1] - 1)]
        idx <- idx + counts[w-1]
      }
    }
    
    # #Even split among remaining weeks
    # if(weeks > 1 && length(remaining) > 0) {
    #   if(weeks == 2) {
    #     spawner_groups[[2]] <- remaining
    #   } else {
    #     splits <- split(
    #       remaining,
    #       cut(seq_along(remaining), weeks - 1, labels = FALSE)
    #     )
    #     for(w in 2:weeks) {
    #       spawner_groups[[w]] <- splits[[w - 1]]
    #     }
    #   }
    # }
  }
  
  ## global state across weeks
  used_id <- setNames(rep(FALSE, length(ids_all)), ids_all) #exclude previously selected fish
  fam_count <- integer(0) #count # of times family has been selected
  used_fam <- new.env(parent=emptyenv()) #exclude families selected 4x
  
  sel_a <- character(0)
  sel_b <- character(0)
  available_ids <- character(0)
  
  for(w in 1:weeks) {
    
    ## add newly available spawners
    available_ids <- c(available_ids, spawner_groups[[w]])
    
    ## remove already used fish
    available_ids <- available_ids[!used_id[available_ids]]
    
    ## apply loss AFTER week 1
    if(loss > 0 && w > 1 && length(available_ids) > 0) {
      keep <- runif(length(available_ids)) > loss
      available_ids <- available_ids[keep]
    }
    
    if(length(available_ids) < 2) next
    
    weekly_target <- if(w < weeks) crosses_perweek else crosses_lastweek
    
    combos <- t(combn(available_ids, 2))
    kinvals <- mapply(function(a, b) Ks[a, b], combos[,1], combos[,2])
    df <- data.frame(a=combos[,1], b=combos[,2], kin=kinvals,
                     stringsAsFactors=FALSE)
    df <- df[df$kin < kin_thresh, ]
    if(nrow(df) == 0) next
    
    # Save state for this week
    best_sel_n <- 0
    best_sel_a <- character(0)
    best_sel_b <- character(0)
    best_used_id <- NULL
    best_fam_count <- NULL
    best_used_fam <- NULL
    
    for(t in 1:tries) {
      # make local copies of the state
      used_id_local <- used_id
      fam_count_local <- fam_count
      used_fam_local <- new.env(parent=emptyenv())
      list2env(as.list(as.environment(used_fam)), envir=used_fam_local)
      
      sel_a_local <- character(0)
      sel_b_local <- character(0)
      
      ord <- sample(nrow(df))
      
      for(i in ord) {
        if(length(sel_a_local) >= weekly_target) break
        
        x <- df$a[i]; y <- df$b[i]
        if(used_id_local[x] || used_id_local[y]) next
        
        fx <- fam[x]; fy <- fam[y]
        nx <- ifelse(is.na(fam_count_local[fx]),0,fam_count_local[fx])
        ny <- ifelse(is.na(fam_count_local[fy]),0,fam_count_local[fy])
        if(nx + 1 > max.fam || ny + 1 > max.fam) next
        
        key <- ifelse(fx < fy, paste(fx,fy,sep="|"), paste(fy,fx,sep="|"))
        if(exists(key, envir=used_fam_local, inherits=FALSE)) next
        
        # accept pair
        used_id_local[x] <- used_id_local[y] <- TRUE
        assign(key, TRUE, envir=used_fam_local)
        fam_count_local[fx] <- nx + 1
        fam_count_local[fy] <- ny + 1
        
        sel_a_local <- c(sel_a_local, x)
        sel_b_local <- c(sel_b_local, y)
      }
      
      # keep the best try
      if(length(sel_a_local) > best_sel_n) {
        best_sel_n <- length(sel_a_local)
        best_sel_a <- sel_a_local
        best_sel_b <- sel_b_local
        best_used_id <- used_id_local
        best_fam_count <- fam_count_local
        best_used_fam <- used_fam_local
      }
      if(length(sel_a_local) >= weekly_target) break
    }
    
    # commit best pairs for this week
    sel_a <- c(sel_a, best_sel_a)
    sel_b <- c(sel_b, best_sel_b)
    used_id <- best_used_id
    fam_count <- best_fam_count
    used_fam <- best_used_fam
    
  }
  
  if(length(sel_a) == 0) stop("no pairs found")
  
  best <- data.frame(
    Sire=sel_a,
    Dam=sel_b,
    kin=Ks[cbind(sel_a, sel_b)],
    stringsAsFactors=FALSE
  )
  
  if(nrow(best) < target_pairs)
    warning("could only find ", nrow(best),
            " pairs (target ", target_pairs, ")")
  
  best$PC <- seq_len(nrow(best))
  best$yr <- max(spawners$yr)
  
  best[1:min(nrow(best), target_pairs), , drop=FALSE]
}
#Run function
# pairs_new <- choose_pairs_avail(spawners, #list of available spawners
#                                 Ks, #kinship matrix
#                                 target_pairs=312,
#                                 kin_thresh=0.0625, #avoid cousins
#                                 max.fam=4, #number of spawners per family (Mandi limits to 4 in hatchery)
#                                 weeks=4, #number of subsets of spawners (AKA # weeks of spawning)
#                                 early_frac=0.5, #front load spawners in 1st week (to match reality better)
#                                 loss=0.25, #proportion of fish lost each week
#                                 tries=20)


find_ancestors <- function(ped, ids) {
  all_anc <- character(0)
  new_ids <- ids
  while(length(new_ids) > 0) {
    parents <- unique(c(ped$Sire[ped$ID %in% new_ids], ped$Dam[ped$ID %in% new_ids]))
    parents <- parents[!is.na(parents) & !(parents %in% all_anc)]
    all_anc <- c(all_anc, parents)
    new_ids <- parents
  }
  all_anc
}
#Run function
# ancestors <- find_ancestors(ped.sim, rep_ids)



run_forward_sim <- function(input_pedigree,
                            pairs,
                            n.generations=10,
                            target_pairs=312,
                            max.progeny=10,
                            kin_seq=c(0.03125, 0.0625),
                            max.fam=10,
                            weeks=2,
                            loss=0,
                            early_frac=NULL,
                            tries=1) {
  
  ped.sim <- input_pedigree
  ped.new <- input_pedigree
  
  # track which kin_thresh was used per generation
  kin_used <- numeric(0)
  
  for(g in 1:n.generations) {
    message("Generation ", g)
    
    # 1. Simulate progeny
    ped.sim <- simulate_progeny(
      ped.sim,
      pairs=pairs,
      max.progeny=max.progeny
    )
    
    # Identify new spawners
    spawners <- ped.sim[ped.sim$yr == max(ped.sim$yr), ]
    if(nrow(spawners) == 0) {
      message("Stopped: no spawners in generation ", g)
      break
    }
    
    # pick one representative progeny per family (PC)
    rep_ids <- tapply(spawners$ID, spawners$PC, `[`, 1)
    
    # 2. Kinship for representative progeny only
    ancestors <- find_ancestors(ped.sim, rep_ids)
    ped.sub <- ped.sim[ped.sim$ID %in% c(rep_ids, ancestors), ]
    
    K_rep <- kinship(
      id=ped.sub$ID,
      dadid=ped.sub$Sire,
      momid=ped.sub$Dam
    )
    K_rep <- K_rep[rep_ids, rep_ids]
    
    # map each spawner to its representative
    rep_for_spawner <- rep_ids[as.character(spawners$PC)]
    
    # expand kinship back to individual spawners
    Ks <- K_rep[rep_for_spawner, rep_for_spawner]
    rownames(Ks) <- spawners$ID
    colnames(Ks) <- spawners$ID
    
    # try progressively relaxed kinship thresholds
    pairs_new <- NULL
    kin_used_g <- NA
    
    for(k in kin_seq) {
      message("  Trying kin_thresh = ", k)
      
      tmp <- try(
        choose_pairs_avail(
          spawners,
          Ks,
          target_pairs=target_pairs,
          kin_thresh=k,
          max.fam=max.fam,
          weeks=weeks,
          early_frac=early_frac,
          loss=loss,
          tries=tries
        ),
        silent=TRUE
      )
      
      if(!inherits(tmp, "try-error") && nrow(tmp) >= target_pairs) {
        pairs_new <- tmp
        kin_used_g <- k
        break
      }
    }
    
    # if no kin_thresh worked, stop simulation
    if(is.null(pairs_new)) {
      message(
        "Stopped: could not find enough pairs in generation ",
        g,
        " even after relaxing kin_thresh"
      )
      break
    }
    
    kin_used[g] <- kin_used_g
    
    # 4. Update ped.new with selected spawners
    # ped.new <- rbind(
    #   ped.new,
    #   spawners[
    #     spawners$ID %in% pairs_new$Sire |
    #       spawners$ID %in% pairs_new$Dam,
    #   ]
    # )
    
    ## extract selected spawners
    sel <- spawners[
      spawners$ID %in% pairs_new$Sire |
        spawners$ID %in% pairs_new$Dam,
    ]
    
    ## map new PC onto selected spawners
    pc_map <- c(
      setNames(pairs_new$PC, pairs_new$Sire),
      setNames(pairs_new$PC, pairs_new$Dam)
    )
    sel$PC <- pc_map[sel$ID]
    sel <- sel[order(sel$PC),]
    
    ## add to pedigree
    ped.new <- rbind(ped.new, sel)
    
    
    # 5. Update pairs for next generation
    pairs <- pairs_new
    
    # 6. Reset ped.sim
    ped.sim <- ped.new
  }
  
  # Final F calculation
  ped.new[order(ped.new$yr,ped.new$PC),] #order by PC within each year
  ped.new$F_R <- pedInbreeding(ped.new)[,2]
  
  # Attach metadata
  attr(ped.new,"n.generations") <- n.generations
  attr(ped.new,"target_pairs") <- target_pairs
  attr(ped.new,"max.progeny") <- max.progeny
  attr(ped.new,"max.fam") <- max.fam
  attr(ped.new,"weeks") <- weeks
  attr(ped.new,"loss") <- loss
  attr(ped.new,"early_frac") <- early_frac
  attr(ped.new,"tries") <- tries
  attr(ped.new, "kin_used_by_generation") <- kin_used
  
  return(ped.new)
}
#Run function
# run_forward_sim(input_pedigree,pairs,
#                 n.generations=10,
#                 target_pairs=300,
#                 kin_seq=c(0.03125, 0.0625),
#                 n.progeny=10,
#                 n.fam=4,
#                 tries=1)


# #### Manual sim - keep for reference ####
# 
# ## Simulate progeny
# ped.sim <- simulate_progeny(ped.sim,pairs,n.progeny=10)
# 
# ## Identify available spawners
# spawners <- ped.sim[ped.sim$yr==max(ped.sim$yr),]
# 
# ## Compute kinship matrix for spawners
# K <- kinship(id=ped.sim$ID,dadid=ped.sim$Sire,momid=ped.sim$Dam) # full kinship matrix
# Ks <- K[spawners$ID, spawners$ID]   # restrict to spawners only
# 
# ## Simulate pair crosses
# pairs <- choose_pairs(spawners, Ks, target_pairs=300, kin_thresh=0.0625,tries=1)
# 
# ## Update pedigree to only include newly selected spawners
# ped.new <- rbind(ped.new,ped.sim[ped.sim$ID %in% pairs$Sire | ped.sim$ID %in% pairs$Dam,])
# 
# ## Calc F values for updated pedigree
# ped.new$F_R <- pedInbreeding(ped.new)[,2]
